##############
# Functions: #
##############
change_cols_to_numeric = function(numeric_cols,data){
  for (i in 1:length(numeric_cols)){
    data[,numeric_cols[i]] = as.numeric(data[,numeric_cols[i]])
  }
  return(data)
}

create_sites_data = function(curr_codons_table,curr_mat){
  all_sites = unique(curr_codons_table$ref_site)
  sites_data = data.frame(matrix(0,length(all_sites),dim(curr_mat)[2]))
  colnames(sites_data) = colnames(curr_mat)
  for (i in 1:length(all_sites)){
    curr_site = all_sites[i]
    plcs = which(curr_codons_table$ref_site==all_sites[i])
    if (length(plcs)>1){
      sites_data[i,] = apply(curr_mat[plcs,],2,sum)
    }else{
      sites_data[i,] = curr_mat[plcs,]
    }
  }
  ref_site = all_sites
  sites_data = cbind(sites_data,ref_site)
  return(sites_data)
}

true_hits = function(output,prediction){
  data = data.frame(cbind(output,prediction))
  data = data[sample(1:dim(data)[1]),]
  data = data[order(data[,"prediction"],decreasing = TRUE),]
  num_sites = dim(data)[1]
  ratios = matrix(0,num_sites,1)
  for (i in 1:num_sites){
    ratios[i] = length(which(data[1:i,"output"]>0))/i
  }
  # print(paste(c(ratios[c(300,600,900,1200,1500)])))
  return(ratios)
}

get_FP_TP = function(ratios,curr_output){
  true_hits_nums = ratios*(1:length(ratios))
  TP = true_hits_nums/length(which(curr_output>0))
  FP =(c(1:length(ratios))-true_hits_nums)/length(which(curr_output==0))
  return(list(TP,FP))
}

get_plot_vars = function(spike_flag,curr_codons_table,curr_predictions,
                         testing_plcs){
  # curr_codons_table = codons_table_amino_diff_output
  # curr_predictions = amino_predictions_diff_output
  if (spike_flag){
    plcs = which(curr_codons_table$gene=="S")
  }else{
    plcs = c(1:dim(curr_codons_table)[1])
  }
  curr_codons_table = curr_codons_table[plcs,]
  for (j in testing_plcs){
    curr_predictions[[j]] = curr_predictions[[j]][plcs,]
  }
  syn_plcs =which(curr_predictions[[1]][,1]>0)
  non_syn_plcs = which(curr_predictions[[4]][,1]>0)
  return(list(syn_plcs,non_syn_plcs,curr_codons_table,curr_predictions))
}

divide_codons_table_to_outputs = function(curr_codons_table,syn_ti_const,syn_tv_const,
                                          non_syn_ti_const,non_syn_tv_const){
  input_cols = c("L.neighbor","R.neighbor","base","codon","site","ref_site","codon.pos","gene","mat_peptide",
                 "stem_loop","ref_seq","ref_codon","exposure","CG","amino_acid","diff_seqs_num","amino_num")
  output_cols = c("y","syn","non_syn","syn_exposure","non_syn_exposure","transitions","transversions")
  bases = c("A","C","G","T")
  num_states = dim(curr_codons_table)[1]
  codons_table_divided = matrix(0,3*num_states,26)
  cat("\n")
  pb <- txtProgressBar(min = 0, max = num_states, style = 3, width = 50, char = "=")
  for (j in 1:num_states){
    # print(j)
    curr_row = curr_codons_table[j,]
    curr_row_partial = curr_row[,input_cols]
    output_base = bases[is.element(bases,curr_row$base)==FALSE]
    output_amino = curr_row[c(paste0(output_base[1],".amino"),paste0(output_base[2],".amino"),paste0(output_base[3],".amino"))]
    tmp = data.frame(matrix(curr_row_partial,3,length(curr_row_partial),byrow = TRUE))
    colnames(tmp) = colnames(curr_row_partial)
    tmp2 = matrix(0,3,length(output_cols))
    colnames(tmp2) = output_cols
    tmp = cbind(tmp,output_base,t(output_amino),tmp2)
    colnames(tmp)[length(curr_row_partial)+2] = "output_amino"
    for (i in 1:3){
      tmp$y[i] = curr_row[,tmp$output_base[i]]
      AG_sum = sum(is.element(c(tmp$base[i],tmp$output_base[i]),c("A","G")))
      if (AG_sum==0 | AG_sum==2){ # This is a transition
        tmp$transitions[i] = tmp$y[i]
        ti_flag = 1
      }else{ # This is a transversion
        tmp$transversions[i] = tmp$y[i]
        ti_flag = 0
      }
      if (curr_row$amino_acid==tmp$output_amino[i][[1]]){
        # Syn. sub:
        if (ti_flag){
          curr_exposure = syn_ti_const
        }else{
          curr_exposure = syn_tv_const
        }
        tmp$syn[i] = tmp$y[i]
        tmp$syn_exposure[i] = curr_exposure/curr_row$syn_exposure
      }else{
        # Non-syn sub:
        if (ti_flag){
          # curr_exposure = 1
          curr_exposure = non_syn_ti_const
        }else{
          curr_exposure = non_syn_tv_const
        }
        tmp$non_syn[i] = tmp$y[i]
        tmp$non_syn_exposure[i] = curr_exposure/curr_row$non_syn_exposure
      }
    }
    codons_table_divided[((3*j-2):(3*j)),] = unlist(tmp)
    setTxtProgressBar(pb, j)
  }
  colnames(codons_table_divided) = colnames(tmp)
  codons_table_divided = data.frame(codons_table_divided)
  
  codons_table_divided = change_cols_to_numeric(numeric_cols = c("exposure","diff_seqs_num",
                                                                 "y","syn","non_syn","syn_exposure",
                                                                 "non_syn_exposure","transitions","transversions"),
                                                data = codons_table_divided)
  return(codons_table_divided)
}

get_only_rel_states = function(diff_seqs_num_flag,codons_table_testing,codons_table_testing_divided){
  if (diff_seqs_num_flag){
    codons_table_testing = codons_table_testing[which(codons_table_testing$diff_seqs_num==0 |
                                                        codons_table_testing$diff_seqs_num>1),]
    codons_table_testing_divided = 
      codons_table_testing_divided[which(codons_table_testing_divided$diff_seqs_num==0 |
                                           codons_table_testing_divided$diff_seqs_num>1),]
  }else{
    codons_table_testing = codons_table_testing[which(codons_table_testing$diff_seqs_num==0 |
                                                        codons_table_testing$diff_seqs_num==1),]
    codons_table_testing_divided = 
      codons_table_testing_divided[which(codons_table_testing_divided$diff_seqs_num==0 |
                                           codons_table_testing_divided$diff_seqs_num==1),]
  }
  return(list(codons_table_testing,codons_table_testing_divided))
}

handle_predictions_na = function(testing_plcs,predictions,internal_predictions,num_models){
  for (i in testing_plcs){
    predictions[[i]] = predictions[[i]][,1:num_models]
    internal_predictions[[i]] = internal_predictions[[i]][,1:num_models]
    for (j in 1:dim(predictions[[i]])[2]){
      predictions[[i]][which(is.na(predictions[[i]][,j])),j] = 0
      internal_predictions[[i]][which(is.na(internal_predictions[[i]][,j])),j] = 0
    }
  }
  return(list(predictions,internal_predictions))
}

get_divided_predictions = function(testing_plcs,predictions,codons_table_testing_divided){
  divided_predictions = list()
  for (i in testing_plcs){
    tmp = t(cbind(predictions[[i]],predictions[[i]],predictions[[i]]))
    tmp = matrix(tmp,3*dim(predictions[[i]])[1],dim(predictions[[i]])[2],byrow = TRUE)
    if (i==1 | i==2){
      curr_exposure = codons_table_testing_divided$syn_exposure
    }else if (i==4 |i==5){
      curr_exposure = codons_table_testing_divided$non_syn_exposure
    }
    for (j in 1:dim(tmp)[2]){
      tmp[,j] = tmp[,j]*curr_exposure
    }
    divided_predictions[[i]] = tmp
  }
  return(divided_predictions)
}

get_codons_table_amino = function(testing_plcs,codons_table_testing_divided,divided_predictions,
                                  num_models){
  cat("\n")
  print("Creating codons_table_amino")
  amino_predictions = list()
  states = apply(codons_table_testing_divided[,c("gene","amino_num","amino_acid","output_amino")],1,paste0,collapse = ".")
  unique_states = unique(states)
  colnames_codons_table_amino = c("gene","amino_num","amino","output_amino","mat_peptide","stem_loop",
                                  "ref_codon","exposure","y","syn","non_syn","transitions","transversions")
  codons_table_amino = matrix(0,length(unique_states),length(colnames_codons_table_amino))
  colnames(codons_table_amino) = colnames_codons_table_amino
  tmp = matrix(0,length(unique_states),num_models)
  row.names(tmp) = row.names(codons_table_amino) = unique_states
  amino_predictions[[1]] = amino_predictions[[2]] = amino_predictions[[4]] = amino_predictions[[5]] = tmp
  pb <- txtProgressBar(min = 0, max = length(unique_states), style = 3, width = 50, char = "=")
  for (i in 1:length(unique_states)){
    # print(i)
    plcs = which(states==unique_states[i])
    codons_table_amino[i,c("gene","amino_num","amino","output_amino")] = strsplit(unique_states[i],split = "\\.")[[1]] 
    codons_table_amino[i,c("mat_peptide","stem_loop","ref_codon")] = 
      unlist(codons_table_testing_divided[plcs[1],c("mat_peptide","stem_loop","ref_codon")])
    codons_table_amino[i,c("exposure","y","syn","non_syn","transitions","transversions")] = 
      apply(codons_table_testing_divided[plcs,c("exposure","y","syn","non_syn","transitions","transversions")],2,sum)
    for (j in testing_plcs){
      tmp = divided_predictions[[j]][plcs,]
      if (length(plcs)>1){
        tmp = apply(tmp,2,sum)
      }
      amino_predictions[[j]][i,] = tmp
    }
    setTxtProgressBar(pb, i)
  }
  codons_table_amino = data.frame(codons_table_amino)
  codons_table_amino = change_cols_to_numeric(numeric_cols = c("exposure","y","syn","non_syn","transitions","transversions"),
                                              data = codons_table_amino)
  return(codons_table_amino)
}

get_codons_table_amino_diff_output = function(testing_plcs,codons_table_testing_divided,
                                              divided_predictions,num_models){
  cat("\n")
  print("Creating codons_table_amino_diff_output")
  amino_predictions_diff_output = list()
  states = apply(codons_table_testing_divided[,c("gene","amino_num","output_amino")],1,paste0,collapse = ".")
  unique_states = unique(states)
  colnames_codons_table_amino_diff_output = c("gene","amino_num","output_amino","mat_peptide","stem_loop","ref_seq",
                                              "ref_codon","exposure","y","syn","non_syn","transitions","transversions")
  codons_table_amino_diff_output = matrix(0,length(unique_states),length(colnames_codons_table_amino_diff_output))
  colnames(codons_table_amino_diff_output) = colnames_codons_table_amino_diff_output
  tmp = matrix(0,length(unique_states),num_models)
  row.names(tmp) = row.names(codons_table_amino_diff_output)= unique_states
  amino_predictions_diff_output[[1]] = amino_predictions_diff_output[[2]] =
    amino_predictions_diff_output[[4]] = amino_predictions_diff_output[[5]] = tmp
  pb <- txtProgressBar(min = 0, max = length(unique_states), style = 3, width = 50, char = "=")
  for (i in 1:length(unique_states)){
    # print(i)
    plcs = which(states==unique_states[i])
    codons_table_amino_diff_output[i,c("gene","amino_num","output_amino")] = strsplit(unique_states[i],split = "\\.")[[1]] 
    codons_table_amino_diff_output[i,c("mat_peptide","stem_loop","ref_seq","ref_codon")] = 
      unlist(codons_table_testing_divided[plcs[1],c("mat_peptide","stem_loop","ref_seq","ref_codon")])
    codons_table_amino_diff_output[i,c("exposure","y","syn","non_syn","transitions","transversions")] = 
      apply(codons_table_testing_divided[plcs,c("exposure","y","syn","non_syn","transitions","transversions")],2,sum)
    for (j in testing_plcs){
      tmp = divided_predictions[[j]][plcs,]
      if (length(plcs)>1){
        tmp = apply(tmp,2,sum)
      }
      amino_predictions_diff_output[[j]][i,] = tmp
    }
    setTxtProgressBar(pb, i)
  }
  
  codons_table_amino_diff_output = data.frame(codons_table_amino_diff_output)
  codons_table_amino_diff_output = change_cols_to_numeric(numeric_cols = c("exposure","y","syn","non_syn","transitions","transversions"),
                                                          data = codons_table_amino_diff_output)
  return(list(codons_table_amino_diff_output,amino_predictions_diff_output))
}

get_ratios = function(testing_plcs,curr_codons_table,syn_plcs,non_syn_plcs,
                      curr_predictions,spike_flag,num_models){
  pb <- txtProgressBar(min = 0, max = length(testing_plcs)*num_models, style = 3, width = 50, char = "=")
  ratios = list()
  count = 0
  for (j in testing_plcs){
    if (j<4){
      curr_output = as.numeric(curr_codons_table$syn[syn_plcs])
      plcs = syn_plcs
      if (j==1){
        curr_name = "syn, NB"
      }else{
        curr_name = "syn, Poisson"
      }
    }else{
      curr_output = as.numeric(curr_codons_table$non_syn[non_syn_plcs])
      plcs = non_syn_plcs
      if (j==4){
        curr_name = "non-syn, NB"
      }else{
        curr_name = "non-syn, Poisson"
      }
    }
    ratios[[j]] = matrix(0,length(plcs),num_models)
    delta = round(length(plcs)/100)
    plot_plcs = delta*c(1:99)
    for (i in 1:num_models){
      # print(paste("model #",i," ",curr_name,"spike:",spike_flag))
      tmp = true_hits(output = curr_output,
                      prediction = as.numeric(curr_predictions[[j]][plcs,i]))
      ratios[[j]][,i] = tmp
      count = count+1
      setTxtProgressBar(pb, count)
    }
  }
  return(ratios)
}

get_results = function(testing_plcs,curr_codons_table,syn_plcs,non_syn_plcs,ratios,num_models){
  results = NULL
  for(j in testing_plcs){
    if (j<4){
      curr_output = as.numeric(curr_codons_table$syn[syn_plcs])
      plcs = syn_plcs
      if (j==1){
        curr_name = "syn,NB"
      }else{
        curr_name = "syn,Poisson"
      }
    }else{
      curr_output = as.numeric(curr_codons_table$non_syn[non_syn_plcs])
      plcs = non_syn_plcs
      if (j==4){
        curr_name = "non-syn,NB"
      }else{
        curr_name = "non-syn,Poisson"
      }
    }
    delta = round(length(plcs)/100)
    plot_plcs = delta*c(1:98)
    for (i in c(1:num_models)){
      tmp = ratios[[j]][,i]
      max_lift = 1/(length(which(curr_output>0))/length(curr_output))
      
      # plot(plot_plcs,tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)),
      #      xlab = "Index",ylab = "Lift")
      # points(plot_plcs,matrix(1,length(plot_plcs),1),type = "l",lty = 2)
      # points(plot_plcs,tmp[plot_plcs]/ratios[[j]][plot_plcs,1],col = "blue")
      # title(paste("model num =",i," ",curr_name))
      # legend(plot_plcs[55],max_lift/2,legend = c("constant null","all exposures null","time exposure null"),
      #        col = c("black","blue","red"), lty=1, cex=0.8)
      
      tmp2 = get_FP_TP(ratios = tmp,curr_output = curr_output)
      TP = tmp2[[1]]; FP = tmp2[[2]]
      tmp_all_exposures = get_FP_TP(ratios = ratios[[j]][,1],curr_output = curr_output)
      TP_all_exposures = tmp_all_exposures[[1]]; FP_all_exposures = tmp_all_exposures[[2]]
      tmp_time_exposure = get_FP_TP(ratios = ratios[[j]][,2],curr_output = curr_output)
      TP_time_exposure = tmp_time_exposure[[1]]; FP_time_exposure = tmp_time_exposure[[2]]
      
      AUC = sum(diff(c(0,FP))*TP)
      AUC_all_exposures = sum(diff(c(0,FP_all_exposures))*TP_all_exposures)
      AUC_time_exposure = sum(diff(c(0,FP_time_exposure))*TP_time_exposure)
      
      curr_results = c(curr_name,AUC,
                       (tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)))[3],
                       (tmp[plot_plcs]/ratios[[j]][plot_plcs,1])[3])
      
      results = rbind(results, curr_results)
    }
  }
  colnames(results) = c("type","AUC","Lift[3]","Lift_base[3]")
  return(results)
}

get_plots = function(i,j,spike_flag,ratios,syn_plcs,non_syn_plcs,curr_codons_table,curr_predictions){
  if (j<3){
    plcs = syn_plcs
    curr_output = as.numeric(curr_codons_table$syn[syn_plcs])
  }else{
    plcs = non_syn_plcs
    curr_output = as.numeric(curr_codons_table$non_syn[non_syn_plcs])
  }
  
  delta = round(length(plcs)/100)
  plot_plcs = delta*c(1:98)
  
  if (spike_flag==1){
    title_name = "Spike gene:"
  }else{
    title_name = "All genes:"
  }
  
  tmp = ratios[[j]][,i]
  plot_data = data.frame(cbind(0.01*c(1:98),tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)),
                               matrix(1,length(plot_plcs),1),tmp[plot_plcs]/ratios[[j]][plot_plcs,1]))
  colnames(plot_data) = c("Index","model","ones","null_lift")
  lift_plot = ggplot(plot_data,aes(x=Index,y=Lift))+
    geom_point(aes(y=model,col = "model"),size = 2)+
    geom_point(aes(y=null_lift,col = "null_lift"),size = 2)+
    geom_line(aes(y=ones),linetype = "longdash",size = 1.5,color = "black")+
    scale_y_continuous(breaks=c(0:8),limits = c(0,8))+
    # theme_classic()+ggtitle(paste(title_name,"Lift curve"))+xlab("Proportion of states")+ylab("Lift")+
    theme_classic()+xlab("Proportion of states")+ylab("Lift")+
    theme(plot.title = element_text(hjust = 0.5, size = 20,face = "bold"),axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"), legend.title = element_blank(),
          legend.position = c(.6, .8),legend.text=element_text(size=12))+
    scale_color_hue(labels=c("Winning model vs. random model", "Winning model vs. base model"))
  
  tmp2 = get_FP_TP(ratios = tmp,curr_output = curr_output)
  TP = tmp2[[1]]; FP = tmp2[[2]]
  tmp_all_exposures = get_FP_TP(ratios = ratios[[j]][,1],curr_output = curr_output)
  TP_all_exposures = tmp_all_exposures[[1]]; FP_all_exposures = tmp_all_exposures[[2]]
  tmp_time_exposure = get_FP_TP(ratios = ratios[[j]][,2],curr_output = curr_output)
  TP_time_exposure = tmp_time_exposure[[1]]; FP_time_exposure = tmp_time_exposure[[2]]
  
  AUC = sum(diff(c(0,FP))*TP)
  AUC_all_exposures = sum(diff(c(0,FP_all_exposures))*TP_all_exposures)
  AUC_time_exposure = sum(diff(c(0,FP_time_exposure))*TP_time_exposure)
  
  AUC_data = data.frame(cbind(FP,TP,FP_all_exposures,TP_all_exposures))
  roc_plot = ggplot(AUC_data,aes(x=FP,y=AUC,color = "variable"))+
    geom_line(aes(x = FP_all_exposures,y=TP_all_exposures,col = "TP_all_exposures"),size = 1.5)+
    geom_line(aes(x = FP,y=TP,col = "TP"),size = 1.8)+
    theme_classic()+ggtitle(paste(title_name,"ROC curve"))+xlab("False positive rate")+ylab("True positive rate")+
    theme(plot.title = element_text(hjust = 0.5, size = 20,face = "bold"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"), legend.title = element_blank(),
          legend.position = c(.7, .6),legend.text=element_text(size=12))+
    scale_color_hue(labels=c(paste("Winning model AUC:",signif(AUC,3)),
                             paste("Base model AUC:",signif(AUC_all_exposures,3))))
  
  return(list(lift_plot,roc_plot,
              max(max(plot_data$model[which(plot_data$model<Inf)]), 
                  max(plot_data$null_lift[which(plot_data$null_lift<Inf)]))))
}
#########################
analyse_predcitions_preprocess = function(begin_path){
  library("ggplot2")
  library("gtable")
  library("grid")
  library("gridExtra")
  
  setwd(paste0(begin_path,"/testing/vars/"))
  load("predictions")
  load("internal_predictions")
  load("debug.nb.converged")
  load("debug.p.converged")
  setwd(paste0(begin_path,"vars/"))
  load("codons_table_testing")
  load("syn_ti_const")
  load("syn_tv_const")
  load("non_syn_ti_const")
  load("non_syn_tv_const")
  load("codons_table")
  load("site_details")
  
  diff_seqs_num_flag = 1
  testing_plcs = c(1,2,4,5)
  num_models = dim(predictions[[1]])[2]
  
  codons_table_testing_divided = divide_codons_table_to_outputs(codons_table_testing,
                                                                syn_ti_const,syn_tv_const,
                                                                non_syn_ti_const,non_syn_tv_const)
  tmp = get_only_rel_states(diff_seqs_num_flag,codons_table_testing,codons_table_testing_divided)
  codons_table_testing = tmp[[1]]; codons_table_testing_divided = tmp[[2]]; rm(tmp)
  
  tmp = handle_predictions_na(testing_plcs,predictions,internal_predictions,num_models)
  predictions = tmp[[1]]; internal_predictions = tmp[[2]]; rm(tmp)
  divided_predictions = get_divided_predictions(testing_plcs,predictions,codons_table_testing_divided)
  codons_table_amino = get_codons_table_amino(testing_plcs,codons_table_testing_divided,
                                              divided_predictions,num_models)
  tmp = get_codons_table_amino_diff_output(testing_plcs,codons_table_testing_divided,
                                           divided_predictions,num_models)
  codons_table_amino_diff_output = tmp[[1]]; amino_predictions_diff_output = tmp[[2]]
  prediction_names = names(predictions)
  
  setwd(paste0(begin_path,"/vars/"))
  save(codons_table_amino_diff_output,file = "codons_table_amino_diff_output")
  save(amino_predictions_diff_output,file = "amino_predictions_diff_output")
  save(testing_plcs,file = "testing_plcs")
  save(prediction_names,file = "prediction_names")
  save(num_models,file = "num_models")
}

AUC_and_Lift_plot = function(spike_flag,codons_table_amino_diff_output,
                             amino_predictions_diff_output,testing_plcs,model_type,
                             prediction_names,examined_model_num,num_models){
  tmp = get_plot_vars(spike_flag,codons_table_amino_diff_output,amino_predictions_diff_output,
                      testing_plcs)
  syn_plcs = tmp[[1]]; non_syn_plcs = tmp[[2]];
  ratios = get_ratios(testing_plcs,codons_table_amino_diff_output,syn_plcs,non_syn_plcs,
                      amino_predictions_diff_output,spike_flag,num_models)
  results = get_results(testing_plcs,codons_table_amino_diff_output,syn_plcs,non_syn_plcs,
                        ratios,num_models)
  
  j = which(prediction_names==model_type)
  
  if (spike_flag){
    tmp = get_plots(i=(examined_model_num+2),j,spike_flag=1,ratios,syn_plcs,non_syn_plcs,
                    codons_table_amino_diff_output,amino_predictions_diff_output)
    lift_plot = tmp[[1]]; roc_plot = tmp[[2]]; y_lim = tmp[[3]]
  }else{
    tmp = get_plots(i=(examined_model_num+2),j,spike_flag=0,ratios,syn_plcs,non_syn_plcs,
                    codons_table_amino_diff_output,amino_predictions_diff_output)
    lift_plot = tmp[[1]]; roc_plot = tmp[[2]]; y_lim = tmp[[3]]
    
  }
  plot(lift_plot)
  return(list(results,lift_plot))
}
# grid.arrange(all_lift_plot + ylim(0,max(all_y_lim,spike_y_lim)),spike_lift_plot + 
#                ylim(0,max(all_y_lim,spike_y_lim)),ncol = 2)
# grid.arrange(all_roc_plot,spike_roc_plot,ncol = 2)
