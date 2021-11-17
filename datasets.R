##############
# functions: #
##############
create_contrasts = function(k){
  c<-contr.treatment(k)
  my.coding<-matrix(rep(1/k, k*(k-1)), ncol=(k-1))
  my.simple<-c-my.coding
  return(my.simple)
}

add_contrasts = function(data,cols){
  for (i in 1:length(cols)){
    # print(cols[i])
    data[,cols[i]] = as.factor(data[,cols[i]])
    contrasts(data[,cols[i]]) = create_contrasts(length(table(data[,cols[i]])))
  }
  return(data)
}

nb.glm.theta.data = function(theta,data,curr_model_names){
  # print(theta)
  model.nb = glm(curr_model_names,
                 family = negative.binomial(theta),data = data,maxit = 1000)
  logLik_model.nb = logLik(model.nb)
  return(logLik_model.nb)
}

find_model.nb.data = function(curr_data,curr_model_names){
  
  theta <- tryCatch({
    a = optimize(nb.glm.theta.data, data=curr_data,curr_model_names = curr_model_names,
                 interval=c(0.05, 100), maximum=TRUE)
    theta = as.numeric(a[1])
  }, warning = function(war) {
    a = optimize(nb.glm.theta.data, data=curr_data,curr_model_names = curr_model_names,
                 interval=c(0.05, 100), maximum=TRUE)
    theta = as.numeric(a[1])
    return(theta)
  }, error = function(err) {
    print(paste("MY_ERROR find_model.nb.data:  ",err))
    theta = -1
    return(theta)
  }, finally = {
  }) # END tryCatch
  
  # a = optimize(nb.glm.theta.data, data=curr_data, interval=c(0.05, 100), maximum=TRUE)
  # theta = as.numeric(a[1])
  if(theta!=-1){
    model.nb = glm(curr_model_names,
                   family = negative.binomial(theta),data = curr_data,maxit = 1000)
  }else{model.nb = 0}
  
  return(model.nb)
}

find_model.P.data = function(data,curr_model_names){
  model.p <- tryCatch({
    model.p = glm(curr_model_names,data = data,family = "poisson")
    # GLR = 2*(logLik(model.nb)[1]-logLik(model.p)[1])
  }, warning = function(war) {
    model.p = glm(curr_model_names,data = data,family = "poisson")
    # GLR = 2*(logLik(model.nb)[1]-logLik(model.p)[1])
    return(model.p)
  }, error = function(err) {
    print(paste("MY_ERROR find_model.P.data:  ",err))
    model.p = -1
    return(model.p)
  }, finally = {
  }) # END tryCatch
  
  return(model.p)
}

remove_linear_dependant_cols = function(curr_data,log_plcs){
  #recursive!
  rankifremoved <- sapply(1:ncol(curr_data), function (x) qr(curr_data[,-x])$rank)
  if (length(unique(rankifremoved))>1){
    col_to_remove = which(rankifremoved == max(rankifremoved))
    col_to_remove = col_to_remove[1]   #removes only 1 col each time and makes sure it's not the output
    if(col_to_remove>=(ncol(curr_data)-length(log_plcs)) ){
      return(curr_data) # doesn't remove output and exposure cols
    }else{
      curr_data = curr_data[,-col_to_remove]
      remove_linear_dependant_cols(curr_data,log_plcs)
    }
  }else{
    return(curr_data) 
  }
}

prepare_data_rm_cols = function(curr_data){
  if (colnames(curr_data)[dim(curr_data)[2]]=="syn" & length(unique(curr_data$syn_exposure))>1 ){
    curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
                                        ".-exposure-syn_exposure+offset(log(exposure)+log(syn_exposure))"))
    # curr_data = curr_data[,-which(colnames(curr_data)=="non_syn_exposure")]
    log_plcs = c(which(colnames(curr_data)=="syn_exposure"), which(colnames(curr_data)=="exposure"))
  }else if (colnames(curr_data)[dim(curr_data)[2]]=="non_syn" & length(unique(curr_data$non_syn_exposure))>1){
    curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
                                        ".-exposure-non_syn_exposure+offset(log(exposure)+log(non_syn_exposure))"))
    # curr_data = curr_data[,-which(colnames(curr_data)=="syn_exposure")]
    log_plcs = c(which(colnames(curr_data)=="non_syn_exposure"), which(colnames(curr_data)=="exposure"))
  }else{
    # curr_data = curr_data[,-c(which(colnames(curr_data)=="syn_exposure"),
    # which(colnames(curr_data)=="non_syn_exposure"))]
    curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
                                        " . -exposure + offset(log(exposure))"))
    log_plcs = which(colnames(curr_data)=="exposure")
  }
  if(dim(curr_data)[2]>4){
    curr_data_input = curr_data[,1:(ncol(curr_data)-4),drop=FALSE ] # remove exposure, syn_exposure, non_syn_exposure and output
    factor_cols = which(sapply(curr_data_input, class)=="factor")
    for (i in factor_cols){
      curr_data_input[,i] = as.numeric(factor(curr_data_input[,i]))
    }
    curr_data_factored = cbind(curr_data_input,log(curr_data[,log_plcs]),curr_data[,ncol(curr_data)])
  }else{
    curr_data_factored = cbind(log(curr_data[,log_plcs]),curr_data[,ncol(curr_data)])
  }
  
  colnames(curr_data_factored)[(ncol(curr_data)-3):ncol(curr_data_factored)] = colnames(curr_data)[c(log_plcs,ncol(curr_data))]
  # if((ncol(curr_data)-4)==1){
  #   colnames(curr_data_factored)[1] = colnames(curr_data)[1]
  # }
  curr_data_factored = remove_linear_dependant_cols(curr_data_factored,log_plcs)
  curr_data = curr_data[,colnames(curr_data_factored)]
  
  tmp = list(curr_data,curr_model_names)
  return(tmp)
}

find_df = function(model){
  df = as.numeric(0.5*model$aic+logLik(model))
  return(df)
}

need_to_divide.GLR = function(model_nums,details,num_missing_models){
  sep_logLik = 0
  sep_df = num_missing_models
  model_nums = as.numeric(model_nums)
  for (r in model_nums[1:(length(model_nums)-1)]){
    sep_logLik = sep_logLik + as.numeric(details[r,"logLik"])
    sep_df = sep_df + as.numeric(details[r,"df"])
  }
  GLR = 2*(sep_logLik-as.numeric(details[model_nums[length(model_nums)],"logLik"]))
  df_diff = sep_df-as.numeric(details[model_nums[length(model_nums)],"df"])
  
  GOF = (1-pchisq(GLR, df=df_diff))
  return(GOF)
}

need_to_divide.AIC = function(model_nums,details,num_missing_models){
  # if models are missing (because y was the same for all the subset), logLik=0 and df=1 for
  # every missing model.
  sep_logLik = 0
  sep_df = num_missing_models
  model_nums = as.numeric(model_nums)
  for (r in model_nums[1:(length(model_nums)-1)]){
    sep_logLik = sep_logLik + as.numeric(details[r,"logLik"])
    sep_df = sep_df + as.numeric(details[r,"df"])
  }
  AIC.sep = 2*sep_df-2*sep_logLik
  AIC.all = as.numeric(details[model_nums[length(model_nums)],"AIC"])
  need_to_divide = AIC.sep<AIC.all
  # print(paste("AIC.sep = ",AIC.sep," AIC.all = ",AIC.all))
  return(need_to_divide)
}

get_theta = function(model.nb){
  tmp = strsplit(model.nb$family[[1]],"")
  tmp2 = tmp[[1]][19:(length(tmp[[1]])-1)]
  theta = as.numeric(paste(tmp2,collapse = ""))
  return(theta)
}

binary_matrix = function(data,input_iterate){
  tmp1 = NULL
  for (i in 1:length(input_iterate)){
    curr_col = input_iterate[i]
    vals = sort(unique(data[,curr_col]))
    # print(curr_col)
    # print(vals)
    tmp = matrix(0,dim(data)[1],(length(vals)+1))
    colnames(tmp) = paste(curr_col,1:(length(vals)+1),sep = ".")
    for (j in 1:length(vals)){
      tmp[,j] = data[,curr_col]==vals[j]
    }
    tmp[,(j+1)] = 1
    tmp1 = cbind(tmp1,tmp)
  }
  syn_exposure. = data[,"syn_exposure"]>0
  non_syn_exposure. = data[,"non_syn_exposure"]>0
  tmp1 = cbind(tmp1,syn_exposure.,non_syn_exposure.)
  return(as.matrix(tmp1))
}

process_data = function(curr_data,selected_input,selected_output,curr_details,count,output_num){
  curr_details[,colnames(selected_input)] = selected_input
  curr_details[,"output"] = output_num
  curr_details[,"rows_num"] = dim(curr_data)[1]
  curr_details[,"num_non_0_rows_1"] = (length(which(curr_data[,dim(curr_data)[2]]!=0))+1)
  curr_details[,"output_sum"] = sum(curr_data[,dim(curr_data)[2]])
  if (curr_details[,"output_sum"]==0 | dim(curr_data)[1]==1){
    curr_details[,c("logLik","P.logLik")] = 0
    curr_details[,c("df","P.df")] = c(2,1)
    curr_details[,c("AIC","P.AIC")] = c(4,2)
    curr_details[,c("theta","NB.converged","P.converged","P.test","P.GLR")] = NA
    curr_details[,"cols_num"] = dim(curr_data)[2]
  }else{   
    tmp = prepare_data_rm_cols(curr_data)
    curr_data = tmp[[1]]; curr_model_names = tmp[[2]]
    model.nb = find_model.nb.data(curr_data,curr_model_names)
    model.p = find_model.P.data(curr_data, curr_model_names)
    curr_details[,"cols_num"] = dim(curr_data)[2]
    if (is.double(model.nb)==FALSE){
      curr_details[,"logLik"] = logLik(model.nb)
      curr_details[,"df"] = find_df(model.nb) + 1 #+1 because theta is found outside of the glm
      curr_details[,"AIC"] = AIC(model.nb) + 2 #+2 because theta is found outside of the glm
      curr_details[,"theta"] = get_theta(model.nb)
      curr_details[,"NB.converged"] = model.nb$converged
    }else{
      # print(selected_input)
      curr_details[,"logLik"] = NA
      curr_details[,"df"] = NA
      curr_details[,"AIC"] = NA
      curr_details[,"theta"] = NA
      curr_details[,"NB.converged"] = NA
    }
    if (is.double(model.p)==FALSE){
      if (is.double(model.nb)==FALSE){
        GLR.model.p = 2*(logLik(model.nb)[1]-logLik(model.p)[1])
        test.model.p = (1-pchisq(GLR.model.p, df=1))
      }else{
        GLR.model.p = NA
        test.model.p = NA
      }
      logLik.model.p = logLik(model.p)
      df.model.p = find_df(model.p)
      AIC.model.p = AIC(model.p)
      converged.p = model.p$converged 
    }else{
      GLR.model.p = NA
      test.model.p = NA
      logLik.model.p = NA 
      df.model.p = NA
      AIC.model.p = NA
      converged.p = NA
    }
    curr_details[,"P.GLR"] = GLR.model.p
    curr_details[,"P.test"] = test.model.p
    curr_details[,"P.logLik"] = logLik.model.p
    curr_details[,"P.df"] = df.model.p
    curr_details[,"P.AIC"] = AIC.model.p
    curr_details[,"P.converged"] = converged.p
  }
  return(curr_details)
}

is_neighbor_same_as_codon = function(curr_data,rm_cols,selected_input){
  if (is.na(selected_input[1,"codon.pos"])==FALSE){
    if (selected_input[1,"codon.pos"]==1){
      rm_cols = c(rm_cols,which(colnames(curr_data)=="R.neighbor"))
    }
    if (selected_input[1,"codon.pos"]==2){
      rm_cols = c(rm_cols, 
                  which(colnames(curr_data)=="R.neighbor"),
                  which(colnames(curr_data)=="L.neighbor"))
    }
    if (selected_input[1,"codon.pos"]==3){
      rm_cols = c(rm_cols, which(colnames(curr_data)=="L.neighbor"))
    }
  }
  rm_cols = sort(unique(rm_cols))
  return(rm_cols)
}

loop_process = function(curr_data,curr_binary_data,selected_input,iterate_col,x){
  if (iterate_col!="syn_exposure" & iterate_col!="non_syn_exposure"){
    selected_input[,iterate_col] = x
  }
  curr_col = paste(iterate_col,x,sep = ".")
  if (is.null(dim(curr_binary_data)[1])==FALSE){
    plcs = which(curr_binary_data[,curr_col]==1)
    if (length(plcs)>1){
      curr_data = curr_data[plcs,]
      curr_binary_data = curr_binary_data[plcs,]
    }else{
      curr_data = NA
      curr_binary_data = NA
    }
  }else{
    curr_data = NA
    curr_binary_data = NA
  }
  tmp = list(selected_input,curr_binary_data,curr_data)
  return(tmp)
}

datasets_preprocess = function(codons_table){
  codons_table$mat_peptide[which(codons_table$mat_peptide!="0")] = 1
  codons_table$mat_peptide = as.numeric(codons_table$mat_peptide)
  tmp = c("L.neighbor","R.neighbor","codon","codon.pos","gene",            
          "stem_loop","amino_acid","base")
  codons_table = add_contrasts(codons_table,tmp)
  
  # output = c("syn","non_syn","transitions","transversions","y","A","C","G","T")
  output = c("syn","non_syn")
  
  # data:
  input = c("L.neighbor","R.neighbor","codon","site",
            "codon.pos","gene","mat_peptide","stem_loop","CG",             
            "amino_acid","base","exposure", "syn_exposure","non_syn_exposure")
  
  input_iterate = c("L.neighbor","R.neighbor","codon","codon.pos","gene",
                    "mat_peptide","stem_loop","CG","amino_acid","base") 
  data = codons_table
  data_input = data[,input]
  data_output = data[,output]
  binary_data = binary_matrix(data,input_iterate)
  
  colnames_curr_details = c(input_iterate,"output","ID","logLik","df","AIC","theta","NB.converged","rows_num","cols_num",
                            "P.GLR","P.test","P.logLik","P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum")
  # details = matrix(0,500000,length(colnames_curr_details))
  selected_input = matrix(0,1,length(input_iterate))
  colnames(selected_input) = c(input_iterate)
  curr_details = matrix(0,1,length(colnames_curr_details))
  # colnames(details) = colnames(curr_details)
  colnames(curr_details) = colnames_curr_details
  iterate_vals = matrix(0,1,length(input_iterate))
  colnames(iterate_vals) = input_iterate
  for (i in 1:length(input_iterate)){
    curr_col = input_iterate[i]
    iterate_vals[,curr_col] = length(unique(data[,curr_col]))
  }
  iterate_vals = iterate_vals+2
  iterate_vals = data.frame(iterate_vals)
  
  data_rows = matrix(0,dim(data)[1],1)
  rownames(data_rows) = rownames(data)
  
  setwd(paste0(begin_path,"vars/"))
  save(iterate_vals,file = "iterate_vals")
  save(data_rows,file = "data_rows")
  save(output,file = "output")
  save(selected_input,file = "selected_input")
  save(curr_details,file = "curr_details")
  save(data_input,file = "data_input")
  save(data_output,file = "data_output")
  save(binary_data,file = "binary_data")
  save(input_iterate,file = "input_iterate")
}

#############
# DataSets: #
#############
datasets = function(codons_table,q,w,b,j,parallel_num,iterate_vals,data_rows,output,selected_input,
                    curr_details,data_input,data_output,binary_data,input_iterate){
  colnames_short_details =c(input_iterate,"output","reg_name2","count",
                            "u","t","w","q","a","p","k","j","b")
  short_details = NULL
  other_outputs_counts = NULL
  
  count = 0
  for (i in 1:length(output)){
    selected_output = output[i]
    selected_input[1:length(selected_input)]=0
    curr_details[1:length(curr_details)]=0
    curr_data = cbind(data_input,data_output[,i])
    curr_binary_data = binary_data
    i_plcs = NULL
    if (i==6){
      i_plcs = which(curr_data$base=="A")
    }else if (i==7){
      i_plcs = which(curr_data$base=="C") 
    }else if (i==8){
      i_plcs = which(curr_data$base=="G")
    }else if (i==9){
      i_plcs = which(curr_data$base=="T")
    }
    if (length(i_plcs)>0){
      curr_data = curr_data[-i_plcs,]
      curr_binary_data = curr_binary_data[-i_plcs,]
    }
    colnames(curr_data)[dim(curr_data)[2]] = selected_output
    rm_cols = NULL
    if (i<3){
      if (i==1){tmp = loop_process(curr_data,binary_data,selected_input,"syn_exposure","")}
      else{ #i=2
        tmp = loop_process(curr_data,binary_data,selected_input,"non_syn_exposure","")
      }
      selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
    }else{
      # curr_binary_data = binary_data 
    }
    for (j in j){
      curr_loop_data.j = curr_data; curr_loop_binary_data.j = curr_binary_data; rm_cols.j = rm_cols
      if (j==iterate_vals[1,"codon.pos"]){
        rm_cols = c(rm_cols,which(colnames(curr_data)=="codon.pos"))
        selected_input[,"codon.pos"] = NA
      }else{
        tmp = loop_process(curr_data,curr_binary_data,selected_input,"codon.pos",j)
        selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
        if(is.na(tmp[2])){
          curr_data = curr_loop_data.j; curr_binary_data = curr_loop_binary_data.j; rm_cols = rm_cols.j 
          next
        }
      }
      for (b in b){
        curr_loop_data.b = curr_data; curr_loop_binary_data.b = curr_binary_data; rm_cols.b = rm_cols
        if (b==iterate_vals[1,"base"]){
          rm_cols = c(rm_cols,which(colnames(curr_data)=="base"))
          selected_input[,"base"] = NA
        }else{
          tmp = loop_process(curr_data,curr_binary_data,selected_input,"base",b)
          selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
          if(is.na(tmp[2])){
            curr_data = curr_loop_data.b; curr_binary_data = curr_loop_binary_data.b; rm_cols = rm_cols.b 
            next
          }
        }
        for (k in 1:iterate_vals[1,"mat_peptide"]){
          curr_loop_data.k = curr_data; curr_loop_binary_data.k = curr_binary_data; rm_cols.k = rm_cols
          if (k==iterate_vals[1,"mat_peptide"]){
            rm_cols = c(rm_cols,which(colnames(curr_data)=="mat_peptide"))
            selected_input[,"mat_peptide"] = NA
          }else{
            tmp = loop_process(curr_data,curr_binary_data,selected_input,"mat_peptide",k)
            selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
            if(is.na(tmp[2])){
              curr_data = curr_loop_data.k; curr_binary_data = curr_loop_binary_data.k; rm_cols = rm_cols.k 
              next
            }
          }
          for (p in 1:iterate_vals[1,"CG"]){
            curr_loop_data.p = curr_data; curr_loop_binary_data.p = curr_binary_data; rm_cols.p = rm_cols
            if (p==iterate_vals[1,"CG"]){
              rm_cols = c(rm_cols,which(colnames(curr_data)=="CG"))
              selected_input[,"CG"] = NA
            }else{
              tmp = loop_process(curr_data,curr_binary_data,selected_input,"CG",p)
              selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
              if(is.na(tmp[2])){
                curr_data = curr_loop_data.p; curr_binary_data = curr_loop_binary_data.p; rm_cols = rm_cols.p
                next
              }
            }
            for (a in 1:iterate_vals[1,"stem_loop"]){
              curr_loop_data.a = curr_data; curr_loop_binary_data.a = curr_binary_data; rm_cols.a = rm_cols
              if (a==iterate_vals[1,"stem_loop"]){
                rm_cols = c(rm_cols,which(colnames(curr_data)=="stem_loop"))
                selected_input[,"stem_loop"] = NA
              }else{
                tmp = loop_process(curr_data,curr_binary_data,selected_input,"stem_loop",a)
                selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
                if(is.na(tmp[2])){
                  curr_data = curr_loop_data.a; curr_binary_data = curr_loop_binary_data.a; rm_cols = rm_cols.a
                  next
                }
              }
              for (q in q){#in 1:iterate_vals[1,"codon"]){
                curr_loop_data.q = curr_data; curr_loop_binary_data.q = curr_binary_data; rm_cols.q = rm_cols
                if (q==iterate_vals[1,"codon"]){
                  rm_cols = c(rm_cols,which(colnames(curr_data)=="codon"))
                  selected_input[,"codon"] = NA
                  # use amino_acid instead of codons
                  w_max = 1
                }else{
                  tmp = loop_process(curr_data,curr_binary_data,selected_input,"codon",q)
                  selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
                  w_max = iterate_vals[1,"amino_acid"]
                  if(is.na(tmp[2])){
                    curr_data = curr_loop_data.q; curr_binary_data = curr_loop_binary_data.q; rm_cols = rm_cols.q
                    next
                  }
                }
                for (w in w){# w_max:iterate_vals[1,"amino_acid"]){
                  w = max(w,w_max)
                  curr_loop_data.w = curr_data; curr_loop_binary_data.w = curr_binary_data; rm_cols.w = rm_cols
                  if (w==iterate_vals[1,"amino_acid"]){
                    rm_cols = c(rm_cols,which(colnames(curr_data)=="amino_acid"))
                    selected_input[,"amino_acid"] = NA
                  }else{
                    tmp = loop_process(curr_data,curr_binary_data,selected_input,"amino_acid",w)
                    selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
                    if(is.na(tmp[2])){
                      curr_data = curr_loop_data.w; curr_binary_data = curr_loop_binary_data.w; rm_cols = rm_cols.w
                      next
                    }
                  }
                  if (j==1 |j==2){
                    t_max=iterate_vals[1,"R.neighbor"]
                  }else{t_max = 1}
                  for (t in t_max:iterate_vals[1,"R.neighbor"]){
                    curr_loop_data.t = curr_data; curr_loop_binary_data.t = curr_binary_data; rm_cols.t = rm_cols
                    if (t == iterate_vals[1,"R.neighbor"]){
                      rm_cols = c(rm_cols,which(colnames(curr_data)=="R.neighbor"))
                      selected_input[,"R.neighbor"] = NA
                    }else{
                      tmp = loop_process(curr_data,curr_binary_data,selected_input,"R.neighbor",t)
                      selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
                      if(is.na(tmp[2])){
                        curr_data = curr_loop_data.t; curr_binary_data = curr_loop_binary_data.t; rm_cols = rm_cols.t
                        next
                      }
                    }
                    if (j==2 |j==3){
                      u_max=iterate_vals[1,"L.neighbor"]
                    }else{u_max = 1}
                    for (u in u_max:iterate_vals[1,"L.neighbor"]){
                      curr_loop_data.u = curr_data; curr_loop_binary_data.u = curr_binary_data; rm_cols.u = rm_cols
                      if (u == iterate_vals[1,"L.neighbor"]){
                        rm_cols = c(rm_cols,which(colnames(curr_data)=="L.neighbor"))
                        selected_input[,"L.neighbor"] = NA
                      }else{
                        tmp = loop_process(curr_data,curr_binary_data,selected_input,"L.neighbor",u)
                        selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
                        if(is.na(tmp[2])){
                          curr_data = curr_loop_data.u; curr_binary_data = curr_loop_binary_data.u; rm_cols = rm_cols.u
                          next
                        }
                      }
                      for (m in 1:iterate_vals[1,"gene"]){
                        curr_loop_data.m = curr_data; curr_loop_binary_data.m = curr_binary_data; rm_cols.m = rm_cols
                        if (m==iterate_vals[1,"gene"]){
                          rm_cols = c(rm_cols,which(colnames(curr_data)=="gene"))
                          selected_input[,"gene"] = NA
                          # use regions2 instead of regions
                          n_max = 1
                        }else{
                          tmp = loop_process(curr_data,curr_binary_data,selected_input,"gene",m)
                          selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]; 
                          #n_max = iterate_vals[1,"regions2"]
                          if(is.na(tmp[2])){
                            curr_data = curr_loop_data.m; curr_binary_data = curr_loop_binary_data.m; rm_cols = rm_cols.m
                            next
                          }
                        }
                        if (dim(curr_data)[1]!=0){
                          # if(dim(curr_data)[1]>2000){
                          #   tmp = 1
                          # }
                          rm_cols = c(rm_cols,which(colnames(curr_data)=="site"))
                          rel_cols = c(1:(dim(curr_data)[2]-4))
                          for (v in rel_cols){ # so that exposure/syn_exposure/non_syn_exposure and output will not be included!
                            if (length(unique(curr_data[,v]))==1){
                              rm_cols = c(rm_cols,v)
                            }
                          }
                          rm_cols = is_neighbor_same_as_codon(curr_data,rm_cols,selected_input)
                          if (is.na(rm_cols)[1]==FALSE){
                            curr_data = curr_data[,-rm_cols]
                          }
                          if (is.na(dim(curr_data)[1])==FALSE){
                            count = count +1
                            # print(count)
                            data_rows[1:length(data_rows)] = 0
                            data_rows[rownames(curr_data)] = 1
                            reg_name = paste0(c(paste0(data_rows,collapse = ""),
                                                paste0(colnames(curr_data),collapse = ".")),
                                              collapse = ".")
                            reg_name2 = digest(reg_name,algo = "md5")
                            curr_short_details = cbind(selected_input,i,reg_name2,count,
                                                       u,t,w,q,a,p,k,j,b)
                            # if (sum(data$A[data_rows])>0 | sum(data$C[data_rows])>0 | sum(data$G[data_rows])>0){
                            #   other_outputs_counts = c(other_outputs_counts,count)
                            # }
                            colnames(curr_short_details) = colnames_short_details
                            # print(curr_short_details)
                            short_details = rbind(short_details,curr_short_details)
                          }
                          # }
                          # curr_data = curr_loop_data.n; curr_binary_data = curr_loop_binary_data.n; rm_cols = rm_cols.n
                        }
                        curr_data = curr_loop_data.m; curr_binary_data = curr_loop_binary_data.m; rm_cols = rm_cols.m
                      }
                      curr_data = curr_loop_data.u; curr_binary_data = curr_loop_binary_data.u; rm_cols = rm_cols.u
                    }
                    curr_data = curr_loop_data.t; curr_binary_data = curr_loop_binary_data.t; rm_cols = rm_cols.t
                  }
                  curr_data = curr_loop_data.w; curr_binary_data = curr_loop_binary_data.w; rm_cols = rm_cols.w
                }
                curr_data = curr_loop_data.q; curr_binary_data = curr_loop_binary_data.q; rm_cols = rm_cols.q
              }
              curr_data = curr_loop_data.a; curr_binary_data = curr_loop_binary_data.a; rm_cols = rm_cols.a
            }
            curr_data = curr_loop_data.p; curr_binary_data = curr_loop_binary_data.p; rm_cols = rm_cols.p
          }
          curr_data = curr_loop_data.k; curr_binary_data = curr_loop_binary_data.k; rm_cols = rm_cols.k
        }   
        curr_data = curr_loop_data.b; curr_binary_data = curr_loop_binary_data.b; rm_cols = rm_cols.b
      }   
      curr_data = curr_loop_data.j; curr_binary_data = curr_loop_binary_data.j; rm_cols = rm_cols.j
    }
  }
  
  file_name = paste0("short_details.parallel_num.",parallel_num,"_q.",q,"_w.",w,"_b",b,"_j",j,collapse = "")
  setwd(paste0(begin_path,"vars/datasets/"))
  
  if(count>0){
    short_details = short_details[1:(count),]
    colnames(short_details) = c(input_iterate,"output","reg_name2","count",
                                "u","t","w","q","a","p","k","j","b")
    short_details = data.frame(short_details)
    # if(length(other_outputs_counts)>0){
    #   bug_name = paste0("other_outputs_counts.parallel_num.",parallel_num,"_q.",q,"_w.",w,"_b",b,"_j",j,collapse = "")
    #   save(other_outputs_counts,file = bug_name)
    # }
  }else{
    short_details = NULL
  }
  save(short_details,file = file_name)
}