################33
# Functions:
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
  print(paste(c(ratios[c(300,600,900,1200,1500)])))
  return(ratios)
}

get_FP_TP = function(ratios,curr_output){
  true_hits_nums = ratios*(1:length(ratios))
  TP = true_hits_nums/length(which(curr_output>0))
  FP =(c(1:length(ratios))-true_hits_nums)/length(which(curr_output==0))
  return(list(TP,FP))
}

get_plot_vars = function(y_flag,sites_flag,divided_flag,amino_flag,spike_flag,states_flag,
                         sites_codons_table_prediction_divided,sites_divided_predictions,
                         codons_table_prediction_sites,sites_predictions,
                         codons_table_prediction_divided,
                         divided_predictions,
                         codons_table_amino,
                         amino_predictions,
                         codons_table_amino_diff_output,
                         amino_predictions_diff_output,
                         codons_table_prediction,predictions){
  if (sites_flag){
    if (divided_flag){
      curr_codons_table = sites_codons_table_prediction_divided
      curr_predictions = sites_divided_predictions
    }else{
      curr_codons_table = codons_table_prediction_sites
      curr_predictions = sites_predictions
    }
  }else{
    curr_codons_table = codons_table_prediction_divided
    curr_predictions = divided_predictions
  }
  if (amino_flag){
    if (divided_flag){
      curr_codons_table = codons_table_amino_diff_output
      curr_predictions = amino_predictions_diff_output
    }else{
      curr_codons_table = codons_table_amino
      curr_predictions = amino_predictions
    }
  }
  if (states_flag){
    curr_codons_table = codons_table_prediction
    curr_predictions = predictions
  }
  if (spike_flag){
    plcs = which(curr_codons_table$gene=="S")
  }else{
    plcs = c(1:dim(curr_codons_table)[1])
  }
  curr_codons_table = curr_codons_table[plcs,]
  for (j in prediction_plcs){
    curr_predictions[[j]] = curr_predictions[[j]][plcs,]
  }
  syn_plcs =which(curr_predictions[[1]][,1]>0)
  non_syn_plcs = which(curr_predictions[[4]][,1]>0)
  return(list(syn_plcs,non_syn_plcs,curr_codons_table,curr_predictions))
}

##################3
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
# begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/prediction//"))
diff_seqs_num_flag = 1
if (diff_seqs_num_flag==1){
  load("no_sites_normalization.predictions_syn_non_syn.more_diff_seqs")
  load("no_sites_normalization.internal_predictions_syn_non_syn.more_diff_seqs")
  load("no_sites_normalization.debug.nb.converged.syn_non_syn.more_diff_seqs")
  # load("non_syn.no_sites_normalization.predictions_syn_non_syn.more_diff_seqs")
  # load("non_syn.no_sites_normalization.internal_predictions_syn_non_syn.more_diff_seqs")
  # load("non_syn.no_sites_normalization.debug.nb.converged.syn_non_syn.more_diff_seqs")
  diff_seqs_name = "more_diff_seqs"
}else if (diff_seqs_num_flag==0){
  load("predictions_syn_non_syn.1_diff_seq")
  load("internal_predictions_syn_non_syn.1_diff_seq")
  load("debug.nb.syn_non_syn.converged.1_diff_seq")
  diff_seqs_name = "1_diff_seq"
}else{
  load("predictions_syn_non_syn.all_seqs")
  load("internal_predictions_syn_non_syn.all_diff_seqs")
  load("debug.nb.converged.syn_non_syn.all_diff_seqs")
  load("debug.p.converged.syn_non_syn.all_diff_seqs")
  diff_seqs_name = "all_diff_seqs"
}


setwd(paste0(begin_path,"vars/"))
load("codons_table")
load("codons_table_prediction")
load("codons_table_prediction_divided")
load("site_details")

if (diff_seqs_num_flag){
  codons_table_prediction = codons_table_prediction[which(codons_table_prediction$diff_seqs_num==0 |
                                                            codons_table_prediction$diff_seqs_num>1),]
  codons_table_prediction_divided = 
    codons_table_prediction_divided[which(codons_table_prediction_divided$diff_seqs_num==0 |
                                            codons_table_prediction_divided$diff_seqs_num>1),]
}else{
  codons_table_prediction = codons_table_prediction[which(codons_table_prediction$diff_seqs_num==0 |
                                                            codons_table_prediction$diff_seqs_num==1),]
  codons_table_prediction_divided = 
    codons_table_prediction_divided[which(codons_table_prediction_divided$diff_seqs_num==0 |
                                            codons_table_prediction_divided$diff_seqs_num==1),]
}

prediction_plcs = c(1,2,4,5)
# prediction_plcs = c(1,2) #y
max_j = dim(predictions[[1]])[2]
for (i in prediction_plcs){
  predictions[[i]] = predictions[[i]][,1:max_j]
  internal_predictions[[i]] = internal_predictions[[i]][,1:max_j]
  for (j in 1:dim(predictions[[i]])[2]){
    predictions[[i]][which(is.na(predictions[[i]][,j])),j] = 0
    internal_predictions[[i]][which(is.na(internal_predictions[[i]][,j])),j] = 0
  }
}

for(i in 1:max_j){
  print(i)
  print(length(which(is.na(predictions[["check.syn"]][,i]) & is.na(predictions[["check.non_syn"]][,i]))))
}

############3
sites_predictions = list()
sites_internal_predictions = list()
for(i in prediction_plcs){
  print(i)
  sites_predictions[[i]] = create_sites_data(curr_codons_table = codons_table_prediction,
                                             curr_mat = predictions[[i]])
  sites_internal_predictions[[i]] = create_sites_data(curr_codons_table = codons_table,
                                                      curr_mat = internal_predictions[[i]])
}

codons_table_sites = create_sites_data(curr_codons_table = codons_table,
                                       curr_mat = codons_table
                                       [,which(is.element(colnames(codons_table),
                                                          c("exposure","syn_exposure","non_syn_exposure",
                                                            "y","syn","non_syn","A","C","G","T")))])
codons_table_prediction_sites = 
  create_sites_data(curr_codons_table = codons_table_prediction,
                    curr_mat = codons_table_prediction
                    [,which(is.element(colnames(codons_table),
                                       c("exposure","syn_exposure","non_syn_exposure",
                                         "y","syn","non_syn","A","C","G","T")))])

#################3
divided_predictions = list()
for (i in prediction_plcs){
  tmp = t(cbind(predictions[[i]],predictions[[i]],predictions[[i]]))
  tmp = matrix(tmp,3*dim(predictions[[i]])[1],dim(predictions[[i]])[2],byrow = TRUE)
  if (i==1 | i==2){
    curr_exposure = codons_table_prediction_divided$syn_exposure
  }else if (i==4 |i==5){
    curr_exposure = codons_table_prediction_divided$non_syn_exposure
  }
  for (j in 1:dim(tmp)[2]){
    tmp[,j] = tmp[,j]*curr_exposure
  }
  divided_predictions[[i]] = tmp
}

############
amino_predictions = list()
states = apply(codons_table_prediction_divided[,c("gene","amino_num","amino_acid","output_amino")],1,paste0,collapse = ".")
unique_states = unique(states)
colnames_codons_table_amino = c("gene","amino_num","amino","output_amino","mat_peptide","stem_loop",
                                "ref_codon","exposure","y","syn","non_syn","transitions","transversions")
codons_table_amino = matrix(0,length(unique_states),length(colnames_codons_table_amino))
colnames(codons_table_amino) = colnames_codons_table_amino
tmp = matrix(0,length(unique_states),max_j)
row.names(tmp) = row.names(codons_table_amino) = unique_states
amino_predictions[[1]] = amino_predictions[[2]] = amino_predictions[[4]] = amino_predictions[[5]] = tmp
for (i in 1:length(unique_states)){
  print(i)
  plcs = which(states==unique_states[i])
  codons_table_amino[i,c("gene","amino_num","amino","output_amino")] = strsplit(unique_states[i],split = "\\.")[[1]] 
  codons_table_amino[i,c("mat_peptide","stem_loop","ref_codon")] = 
    unlist(codons_table_prediction_divided[plcs[1],c("mat_peptide","stem_loop","ref_codon")])
  codons_table_amino[i,c("exposure","y","syn","non_syn","transitions","transversions")] = 
    apply(codons_table_prediction_divided[plcs,c("exposure","y","syn","non_syn","transitions","transversions")],2,sum)
  for (j in prediction_plcs){
    tmp = divided_predictions[[j]][plcs,]
    if (length(plcs)>1){
      tmp = apply(tmp,2,sum)
    }
    amino_predictions[[j]][i,] = tmp
  }
}
codons_table_amino = data.frame(codons_table_amino)
codons_table_amino = change_cols_to_numeric(numeric_cols = c("exposure","y","syn","non_syn","transitions","transversions"),
                                            data = codons_table_amino)

############
amino_predictions_diff_output = list()
# tmp = unique(apply(codons_table_prediction_divided[,c("gene","amino_num")],1,paste0,collapse = "."))
# a = matrix(unlist(strsplit(tmp,split = "\\.")),length(tmp),2,byrow = TRUE)
# amino_acids = unique(codons_table_prediction_divided$output_amino)
# # Create a matrix,all_states.amino_divided, with all amino acids for each state in tmp
# all_states.amino_divided = matrix(0,length(tmp)*21,3)
# for (i in 1:length(amino_acids)){
#   all_states.amino_divided[((i-1)*length(tmp)+1):(i*length(tmp)),1:2] = a
#   all_states.amino_divided[((i-1)*length(tmp)+1):(i*length(tmp)),3] = amino_acids[i]
# }
# all_states = apply(all_states.amino_divided,1,paste0,collapse = ".")

states = apply(codons_table_prediction_divided[,c("gene","amino_num","output_amino")],1,paste0,collapse = ".")
unique_states = unique(states)
colnames_codons_table_amino_diff_output = c("gene","amino_num","output_amino","mat_peptide","stem_loop","ref_seq",
                                            "ref_codon","exposure","y","syn","non_syn","transitions","transversions")
codons_table_amino_diff_output = matrix(0,length(unique_states),length(colnames_codons_table_amino_diff_output))
colnames(codons_table_amino_diff_output) = colnames_codons_table_amino_diff_output
tmp = matrix(0,length(unique_states),max_j)
row.names(tmp) = row.names(codons_table_amino_diff_output)= unique_states
amino_predictions_diff_output[[1]] = amino_predictions_diff_output[[2]] =
  amino_predictions_diff_output[[4]] = amino_predictions_diff_output[[5]] = tmp
for (i in 1:length(unique_states)){
  print(i)
  plcs = which(states==unique_states[i])
  codons_table_amino_diff_output[i,c("gene","amino_num","output_amino")] = strsplit(unique_states[i],split = "\\.")[[1]] 
  codons_table_amino_diff_output[i,c("mat_peptide","stem_loop","ref_seq","ref_codon")] = 
    unlist(codons_table_prediction_divided[plcs[1],c("mat_peptide","stem_loop","ref_seq","ref_codon")])
  codons_table_amino_diff_output[i,c("exposure","y","syn","non_syn","transitions","transversions")] = 
    apply(codons_table_prediction_divided[plcs,c("exposure","y","syn","non_syn","transitions","transversions")],2,sum)
  for (j in prediction_plcs){
    tmp = divided_predictions[[j]][plcs,]
    if (length(plcs)>1){
      tmp = apply(tmp,2,sum)
    }
    amino_predictions_diff_output[[j]][i,] = tmp
  }
}

# missing_states = all_states[which(is.element(all_states,unique_states)==FALSE)]
# tmp = cbind(matrix(unlist(strsplit(missing_states,split = "\\.")),
#                    length(missing_states),3,byrow = TRUE), 
#             matrix(0,length(missing_states),10))
# colnames(tmp) = colnames(codons_table_amino_diff_output)
# codons_table_amino_diff_output = rbind(codons_table_amino_diff_output,tmp)
codons_table_amino_diff_output = data.frame(codons_table_amino_diff_output)
codons_table_amino_diff_output = change_cols_to_numeric(numeric_cols = c("exposure","y","syn","non_syn","transitions","transversions"),
                                                        data = codons_table_amino_diff_output)

# tmp = matrix(0,length(missing_states),dim(amino_predictions_diff_output[[1]])[2])
# row.names(tmp) = missing_states
# for (j in prediction_plcs){
#   amino_predictions_diff_output[[j]] = rbind(amino_predictions_diff_output[[j]],tmp)
# }
#####################
# Create divided_sites_predictions and codons_table:
sites_codons_table_prediction_divided = NULL
sites_divided_predictions = list()
bases = c("A","C","G","T")
for (j in prediction_plcs){
  tmp = NULL
  print(j)
  for (i in 1:length(bases)){
    curr_base = bases[i]
    plcs = which(codons_table_prediction_divided$output_base==curr_base)
    if (j==1){
      sites_codons_table_prediction_divided = rbind(sites_codons_table_prediction_divided,
                                                    create_sites_data(curr_codons_table = codons_table_prediction_divided[plcs,],
                                                                      curr_mat = 
                                                                        codons_table_prediction_divided[plcs,c("exposure","syn_exposure","non_syn_exposure",
                                                                                                               "y","syn","non_syn")]))
    }
    tmp = rbind(tmp,create_sites_data(curr_codons_table = codons_table_prediction_divided[plcs,],
                                      curr_mat = divided_predictions[[j]][plcs,]))
  }
  sites_divided_predictions[[j]] = tmp
}

############

ratios = list()
y_flag = 0 ; sites_flag = 0; divided_flag = 1; amino_flag = 1; spike_flag = 0; states_flag = 0;
tmp = get_plot_vars(y_flag,sites_flag,divided_flag,amino_flag,spike_flag,states_flag,
                    sites_codons_table_prediction_divided,sites_divided_predictions,
                    codons_table_prediction_sites,sites_predictions,
                    codons_table_prediction_divided,
                    divided_predictions,
                    codons_table_amino,
                    amino_predictions,
                    codons_table_amino_diff_output,
                    amino_predictions_diff_output,
                    codons_table_prediction,predictions)
syn_plcs = tmp[[1]]; non_syn_plcs = tmp[[2]]; curr_codons_table = tmp[[3]]; curr_predictions = tmp[[4]]

for (j in prediction_plcs){
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
  if (y_flag){
    if (j>2){
      break
    }
    curr_output = as.numeric(curr_codons_table$y)
    plcs = 1:dim(curr_codons_table)[1]
    if (j==1){
      curr_name = "y, NB"
    }else{
      curr_name = "y, Poisson"
    }
  }
  ratios[[j]] = matrix(0,length(plcs),max_j)
  delta = round(length(plcs)/100)
  plot_plcs = delta*c(1:99)
  for (i in 1:max_j){
    print(paste("model #",i," ",curr_name,"amino:",amino_flag,"spike:",spike_flag,
                "sites:",sites_flag,"divided:",divided_flag,"states:",states_flag))
    if (y_flag){
      tmp = true_hits(output = curr_output,
                      prediction = (as.numeric(curr_predictions[[j]][,i])+as.numeric(curr_predictions[[(j+3)]][,i])))
    }else{
      tmp = true_hits(output = curr_output,
                      prediction = as.numeric(curr_predictions[[j]][plcs,i]))
    }
    ratios[[j]][,i] = tmp
  }
}
setwd(paste0(begin_path,"vars/prediction/syn_non_syn/"))
save(ratios,file = paste0(c("ratios",diff_seqs_name,"spike",spike_flag,"amino",amino_flag,
                            "sites",sites_flag,"states",states_flag,
                            "y",y_flag,"divided",divided_flag),collapse = "."))


# y_flag = 0 ; sites_flag = 0; divided_flag = 0; amino_flag = 1; spike_flag = 1;
# y_flag = 0 ; sites_flag = 0; divided_flag = 1; amino_flag = 1; spike_flag = 0; 
# y_flag = 0 ; sites_flag = 1; divided_flag = 1; amino_flag = 0; spike_flag = 0;
# y_flag = 0 ; sites_flag = 1; divided_flag = 0; amino_flag = 0; spike_flag = 0;
# y_flag = 1 ; sites_flag = 1; divided_flag = 0; amino_flag = 0; spike_flag = 0;
# y_flag = 0 ; sites_flag = 0; divided_flag = 1; amino_flag = 0; spike_flag = 0;
# y_flag = 0 ; sites_flag = 0; divided_flag = 0; amino_flag = 0; spike_flag = 0; states_flag = 1
y_flag = 0 ; sites_flag = 0; divided_flag = 1; amino_flag = 1; spike_flag = 0; states_flag = 0;
y_flag = 0 ; sites_flag = 0; divided_flag = 1; amino_flag = 1; spike_flag = 1; states_flag = 0;

load(paste0(c("ratios",diff_seqs_name,"spike",spike_flag,"amino",amino_flag,
              "sites",sites_flag,"states",states_flag,
              "y",y_flag,"divided",divided_flag),collapse = "."))
tmp = get_plot_vars(y_flag,sites_flag,divided_flag,amino_flag,spike_flag,states_flag,
                    sites_codons_table_prediction_divided,sites_divided_predictions,
                    codons_table_prediction_sites,sites_predictions,
                    codons_table_prediction_divided,
                    divided_predictions,
                    codons_table_amino,
                    amino_predictions,
                    codons_table_amino_diff_output,
                    amino_predictions_diff_output,
                    codons_table_prediction,predictions)
syn_plcs = tmp[[1]]; non_syn_plcs = tmp[[2]]; curr_codons_table = tmp[[3]]; curr_predictions = tmp[[4]]

results = results2 = results3 = NULL
for(j in prediction_plcs){
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
  if (y_flag){
    if (j>2){
      break
    }
    curr_output = as.numeric(curr_codons_table$y)
    plcs = 1:dim(curr_codons_table)[1]
    if (j==1){
      curr_name = "y,NB"
    }else{
      curr_name = "y,Poisson"
    }
  }
  delta = round(length(plcs)/100)
  plot_plcs = delta*c(1:98)
  for (i in c(1:max_j)){
    tmp = ratios[[j]][,i]
    max_lift = 1/(length(which(curr_output>0))/length(curr_output))
    
    plot(plot_plcs,tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)),
         xlab = "Index",ylab = "Lift")
    points(plot_plcs,matrix(1,length(plot_plcs),1),type = "l",lty = 2)
    points(plot_plcs,tmp[plot_plcs]/ratios[[j]][plot_plcs,1],col = "blue")
    # points(plot_plcs,tmp[plot_plcs]/ratios[[j]][plot_plcs,2],col = "red")
    title(paste("model num =",i," ",curr_name))
    legend(plot_plcs[55],max_lift/2,legend = c("constant null","all exposures null","time exposure null"),
           col = c("black","blue","red"), lty=1, cex=0.8)
    
    tmp2 = get_FP_TP(ratios = tmp,curr_output = curr_output)
    TP = tmp2[[1]]; FP = tmp2[[2]]
    tmp_all_exposures = get_FP_TP(ratios = ratios[[j]][,1],curr_output = curr_output)
    TP_all_exposures = tmp_all_exposures[[1]]; FP_all_exposures = tmp_all_exposures[[2]]
    tmp_time_exposure = get_FP_TP(ratios = ratios[[j]][,2],curr_output = curr_output)
    TP_time_exposure = tmp_time_exposure[[1]]; FP_time_exposure = tmp_time_exposure[[2]]
    
    AUC = sum(diff(c(0,FP))*TP)
    AUC_all_exposures = sum(diff(c(0,FP_all_exposures))*TP_all_exposures)
    AUC_time_exposure = sum(diff(c(0,FP_time_exposure))*TP_time_exposure)
    
    # plot(FP_all_exposures,TP_all_exposures,col = "blue",xlab = "FP",ylab = "TP")
    # points(FP_time_exposure,TP_time_exposure,col = "red")
    # points(FP,TP)
    # legend(0.4,0.4,legend = c(paste("model AUC = ",signif(AUC,5)),
    #                           paste("null all exposures AUC = ",signif(AUC_all_exposures,5)),
    #                           paste("null time exposure AUC = ",signif(AUC_time_exposure,5))),
    #        col = c("black","blue","red"), lty=1, cex=0.8)
    # title(paste("model #",i," ",curr_name,"amino:",amino_flag,"spike:",spike_flag,
    #             "sites:",sites_flag,"divided:",divided_flag,"states:",states_flag))
    
    curr_results = c(i,curr_name,AUC,AUC_all_exposures,AUC_time_exposure,
                     (tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)))[1:7],
                     diff(tmp[plot_plcs]*plot_plcs)[1:7])
    
    curr_results2 = c(i,curr_name,AUC,AUC_all_exposures,AUC_time_exposure,
                     (tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)))[c(1,2,3)],
                     (tmp[plot_plcs]/ratios[[j]][plot_plcs,1])[c(1,2,3)])
    curr_results3 = c(i,curr_name,AUC,
                      (tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)))[2],
                      (tmp[plot_plcs]/ratios[[j]][plot_plcs,1])[2])
    
    print(paste(curr_results))
    print(paste(curr_results2))
    results = rbind(results, curr_results)
    results2 = rbind(results2, curr_results2)
    results3 = rbind(results3, curr_results3)
    # print(AUC)
  }
}
colnames(results) = c("model_num","type","AUC","AUC_null_all","AUC_null_time","Lift[1]","Lift[2]","Lift[3]","Lift[4]","Lift[5]"
                      ,"Lift[6]","Lift[7]","Diffs[1]","Diffs[2]","Diffs[3]","Diffs[4]","Diffs[5]"
                      ,"Diffs[6]","Diffs[7]")
colnames(results2) = c("model_num","type","AUC","AUC_null_all","AUC_null_time","Lift[1]","Lift[2]","Lift[3]","Lift_base[1]",
                      "Lift_base[2]","Lift_base[3]")
colnames(results3) = c("model_num","type","AUC","Lift[5]","Lift_base[5]")
tmp = results2[,c("type","AUC","Lift[3]","Lift_base[3]")]
setwd("/Users/keren/Dropbox/covid/new2/figures/")
write.csv(tmp,file = "results_all_ncbi.csv")


# Debug:
data = cbind(curr_output,curr_predictions[[5]][non_syn_plcs,3],curr_predictions[[6]][non_syn_plcs,3],curr_codons_table[non_syn_plcs,])
data = cbind(curr_output,curr_predictions[[4]][non_syn_plcs,3],curr_predictions[[6]][non_syn_plcs,3],curr_codons_table[non_syn_plcs,])
colnames(data)[1:3] = c("output","prediction","sub_model_num")
data = data[order(data$prediction,decreasing = TRUE),]
setwd(paste0(begin_path,"vars/prediction/syn_non_syn/"))
save(data,file = "non_syn.poisson.model_3")


###############3

get_plots = function(i,j,y_flag,sites_flag,divided_flag,amino_flag,spike_flag,states_flag){
  setwd(paste0(begin_path,"vars/prediction/syn_non_syn"))
  load(paste0(c("ratios",diff_seqs_name,"spike",spike_flag,"amino",amino_flag,
                "sites",sites_flag,"states",states_flag,
                "y",y_flag,"divided",divided_flag),collapse = "."))
  tmp = get_plot_vars(y_flag,sites_flag,divided_flag,amino_flag,spike_flag,states_flag,
                      sites_codons_table_prediction_divided,sites_divided_predictions,
                      codons_table_prediction_sites,sites_predictions,
                      codons_table_prediction_divided,
                      divided_predictions,
                      codons_table_amino,
                      amino_predictions,
                      codons_table_amino_diff_output,
                      amino_predictions_diff_output,
                      codons_table_prediction,predictions)
  syn_plcs = tmp[[1]]; non_syn_plcs = tmp[[2]]; curr_codons_table = tmp[[3]]; curr_predictions = tmp[[4]]
  
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

library("ggplot2")
library("gtable")
library("grid")
library("gridExtra")
tmp = get_plots(i=5,j=5,y_flag=0,sites_flag=0,divided_flag=1,amino_flag=1,spike_flag=1,states_flag=0)
spike_lift_plot = tmp[[1]]; spike_roc_plot = tmp[[2]]; spike_y_lim = tmp[[3]]
tmp = get_plots(i=5,j=5,y_flag=0,sites_flag=0,divided_flag=1,amino_flag=1,spike_flag=0,states_flag=0)
all_lift_plot = tmp[[1]]; all_roc_plot = tmp[[2]]; all_y_lim = tmp[[3]]
grid.arrange(all_lift_plot + ylim(0,max(all_y_lim,spike_y_lim)),spike_lift_plot + 
               ylim(0,max(all_y_lim,spike_y_lim)),ncol = 2)
grid.arrange(all_roc_plot,spike_roc_plot,ncol = 2)

setwd(paste0(begin_path,"figures/"))

# #########
# # Divide to different outputs:
# # Add if 
# begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
# setwd(paste0(begin_path,"vars/"))
# load("codons_table_prediction")
# load("syn_ti_const2")
# load("syn_tv_const2")
# load("non_syn_ti_const2")
# load("non_syn_tv_const2")
# 
# # normalize exposure for all sites:
# input_cols = c("L.neighbor","R.neighbor","base","codon","site","ref_site","codon.pos","gene","mat_peptide",
#                "stem_loop","ref_seq","ref_codon","exposure","CG","amino_acid","diff_seqs_num","amino_num")
# output_cols = c("y","syn","non_syn","syn_exposure","non_syn_exposure","transitions","transversions")
# bases = c("A","C","G","T")
# num_states = dim(codons_table_prediction)[1]
# codons_table_prediction_divided = matrix(0,3*num_states,26)
# for (j in 1:num_states){
#   print(j)
#   curr_row = codons_table_prediction[j,]
#   curr_row_partial = curr_row[,input_cols]
#   output_base = bases[is.element(bases,curr_row$base)==FALSE]
#   output_amino = curr_row[c(paste0(output_base[1],".amino"),paste0(output_base[2],".amino"),paste0(output_base[3],".amino"))]
#   tmp = data.frame(matrix(curr_row_partial,3,length(curr_row_partial),byrow = TRUE))
#   colnames(tmp) = colnames(curr_row_partial)
#   tmp2 = matrix(0,3,length(output_cols))
#   colnames(tmp2) = output_cols
#   tmp = cbind(tmp,output_base,t(output_amino),tmp2)
#   colnames(tmp)[length(curr_row_partial)+2] = "output_amino"
#   for (i in 1:3){
#     tmp$y[i] = curr_row[,tmp$output_base[i]]
#     AG_sum = sum(is.element(c(tmp$base[i],tmp$output_base[i]),c("A","G")))
#     if (AG_sum==0 | AG_sum==2){ # This is a transition
#       tmp$transitions[i] = tmp$y[i]
#       ti_flag = 1
#     }else{ # This is a transversion
#       tmp$transversions[i] = tmp$y[i]
#       ti_flag = 0
#     }
#     if (curr_row$amino_acid==tmp$output_amino[i][[1]]){
#       # Syn. sub:
#       if (ti_flag){
#         curr_exposure = syn_ti_const
#       }else{
#         curr_exposure = syn_tv_const
#       }
#       tmp$syn[i] = tmp$y[i]
#       tmp$syn_exposure[i] = curr_exposure/curr_row$syn_exposure
#     }else{
#       # Non-syn sub:
#       if (ti_flag){
#         # curr_exposure = 1
#         curr_exposure = non_syn_ti_const
#       }else{
#         curr_exposure = non_syn_tv_const
#       }
#       tmp$non_syn[i] = tmp$y[i]
#       tmp$non_syn_exposure[i] = curr_exposure/curr_row$non_syn_exposure
#     }
#   }
#   codons_table_prediction_divided[((3*j-2):(3*j)),] = unlist(tmp)
# }
# colnames(codons_table_prediction_divided) = colnames(tmp)
# codons_table_prediction_divided = data.frame(codons_table_prediction_divided)
# 
# codons_table_prediction_divided = change_cols_to_numeric(numeric_cols = c("exposure","diff_seqs_num",
#                                                                           "y","syn","non_syn","syn_exposure",
#                                                                           "non_syn_exposure","transitions","transversions"),
#                                                          data = codons_table_prediction_divided)
# save(codons_table_prediction_divided,file = "codons_table_prediction_divided")
