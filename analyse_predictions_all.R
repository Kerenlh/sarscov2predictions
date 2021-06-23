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
  # if (amino_flag){
  #   # syn_plcs = which(curr_codons_table$amino==curr_codons_table$output_amino)
  #   # non_syn_plcs = which(curr_codons_table$amino!=curr_codons_table$output_amino)
  #   # syn_plcs = which(curr_predictions[[1]][,1]>0)
  #   # non_syn_plcs = which(curr_predictions[[4]][,1]>0)
  #   syn_plcs = non_syn_plcs = c(1:dim(curr_codons_table)[1]) 
  # }else{
  #   
  # }
  syn_plcs =which(curr_predictions[[1]][,1]>0)
  non_syn_plcs = which(curr_predictions[[4]][,1]>0)
  return(list(syn_plcs,non_syn_plcs,curr_codons_table,curr_predictions))
}

##################3
begin_path = "/Users/keren/Dropbox/covid/new2/"
# begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/prediction/all/"))
load("no_sites_normalization.predictions_syn_non_syn")
load("no_sites_normalization.internal_predictions_syn_non_syn")
load("no_sites_normalization.debug.nb.converged.syn_non_syn")
load("no_sites_normalization.debug.p.converged.syn_non_syn")

setwd(paste0(begin_path,"vars/"))
load("codons_table")
load("codons_table_prediction_all")
load("codons_table_prediction_all_divided")
load("site_details")

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
y_flag = 0 ; sites_flag = 0; divided_flag = 1; amino_flag = 1; spike_flag = 1; states_flag = 0;
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

##########################3
# Supplementary TableS3:
predicted_non_syn_subs_rate = data.frame(cbind(curr_codons_table$amino_num[non_syn_plcs],
                                    curr_codons_table$output_amino[non_syn_plcs],
                                    curr_predictions[[5]][non_syn_plcs,5]))
colnames(predicted_non_syn_subs_rate) = c("amino_num","amino_acid_output","predicted_rate")
predicted_non_syn_subs_rate$predicted_rate = as.numeric(predicted_non_syn_subs_rate$predicted_rate)
predicted_non_syn_subs_rate = predicted_non_syn_subs_rate[order(predicted_non_syn_subs_rate$predicted_rate,
                                                                decreasing = TRUE),]
codons_table_states = apply(codons_table[which(codons_table$gene=="S"),
                                         c("gene","amino_num","amino_acid")],1,paste0,collapse = ".")
codons_table_states = gsub(" ","",codons_table_states)
appears_in_training = is.element(row.names(predicted_non_syn_subs_rate),codons_table_states)

codons_table_prediction_states = apply(codons_table_prediction[which(codons_table_prediction$gene=="S"),
                                         c("gene","amino_num","amino_acid")],1,paste0,collapse = ".")
codons_table_prediction_states = gsub(" ","",codons_table_prediction_states)
appears_in_test = is.element(row.names(predicted_non_syn_subs_rate),codons_table_prediction_states)

predicted_non_syn_subs_rate = cbind(predicted_non_syn_subs_rate,appears_in_training,appears_in_test)
setwd("/Users/keren/Dropbox/covid/new2/figures/")
write.csv(predicted_non_syn_subs_rate,"sub_rates_amino_acids.csv")
####################

non_syn_codons_table = curr_codons_table[non_syn_plcs,]
non_syn_predictions = list()
for (j in prediction_plcs){
  non_syn_predictions[[j]] = curr_predictions[[j]][non_syn_plcs,]
}
sorted_codons_table = non_syn_codons_table[order(non_syn_predictions[[5]][,5],decreasing = TRUE),]

spike_plcs = which(codons_table$gene=="S")
codons_table_spike_variants = NULL
for (i in spike_plcs){
  if (codons_table$non_syn[i]>0){
    for (j in c("A","C","G","T")){
      if (codons_table[i,j]>0){
        codons_table_spike_variants = 
          c(codons_table_spike_variants,paste0(codons_table$gene[i], ".",
                                               codons_table$amino_num[i],".",
                                               codons_table[i,paste0(j,".amino")],collapse = ""))
      }
    }
  }
}
codons_table_spike_variants = unique(codons_table_spike_variants)

setwd(paste0(begin_path,"vars/"))
library(stringr)
gisaid_variants = read.csv("variantswatchlist.csv")
gisaid_variants = strsplit(gisaid_variants$Variant,split = "_")
regexp = "[[:digit:]]+"
gisaid_variants_states = NULL
for (i in 1:length(gisaid_variants)){
  # print(i)
  for (j in 1:length(gisaid_variants[[i]])){
    curr_site = str_extract(gisaid_variants[[i]][j], regexp)
    length_site = length(strsplit(curr_site,split = "")[[1]])
    tmp = strsplit(gisaid_variants[[i]][j],split = "")[[1]]
    curr_amino = paste0(tmp[(length_site+1):length(tmp)],collapse = "")
    if (curr_amino=="X"){
      print(i)
      print(j)
    }
    if (curr_amino!="del" & curr_amino!="X"){
      gisaid_variants_states = c(gisaid_variants_states,paste0("S.",curr_site,".",curr_amino))
    }
  }
}
gisaid_variants_states = unique(gisaid_variants_states)
new_gisaid_states = gisaid_variants_states[which(is.element(gisaid_variants_states,codons_table_spike_variants)==FALSE)]

variants_of_interest = c("S.67.A.V","S.484.E.K","S.614.D.G","S.677.Q.H","S.888.F.L",
                         "S.95.T.I","S.253.D.G","S.477.S.N","S.701.A.V","S.80.D.G",
                         "S.157.F.S","S.452.L.R","S.791.T.I","S.859.T.N","S.950.D.H",
                         "S.484.E.Q","S.142.G.D","S.154.E.K","S.681.P.R","S.1071.Q.H",
                         "S.158.R.G","S.478.T.K","S.950.D.N","S.19.T.R","S.565.F.L",
                         "S.1176.V.F")
variants_of_interest = c("S.67.V","S.484.K","S.614.G","S.677.H","S.888.L",
                         "S.95.I","S.253.G","S.477.N","S.701.V","S.80.G",
                         "S.157.S","S.452.R","S.791.I","S.859.N","S.950.H",
                         "S.484.Q","S.142.D","S.154.K","S.681.R","S.1071.H",
                         "S.158.G","S.478.K","S.950.N","S.19.R","S.565.L",
                         "S.1176.F")
variants_of_concern = c("S.484.E.K","S.494.S.P","S.501.N.Y","S.570.A.D","S.614.D.G",
                        "S.681.P.H","S.716.T.I","S.982.S.A","S.1118.D.H","S.1191.K.N",
                        "S.80.D.A","S.215.D.G","S.417.K.N","S.701.A.V","S.452.L.R","S.13.S.I",
                        "S.152.W.C","S.18.L.F","S.20.T.N","S.26.P.S","S.138.D.Y","S.190.R.S",
                        "S.417.K.T","S.655.H.Y","S.1027.T.I")

variants_from_infectivity_paper = c("S.495.Y.N","S.489.Y.H","S.423.Y.C","S.489.Y.D","S.505.Y.H",
                                    "S.400.F.C","S.421.Y.S","S.489.Y.F","S.423.Y.F","S.501.N.S",
                                    "S.500.T.N","S.511.V.E","S.489.Y.C","S.423.Y.S","S.400.F.V",
                                    "S.338.F.C","S.512.V.L","S.489.Y.N","S.452.L.P","S.460.N.T")

find_variants_plcs = function(sorted_codons_table,variants){
  variants_plcs = data.frame(matrix(0,length(variants),2))
  variants_plcs[,1] = variants
  colnames(variants_plcs) = c("variant","plc")
  for (i in 1:length(variants)){
    curr_plc = which(row.names(sorted_codons_table)==variants[i])
    if (length(curr_plc)>0){
      variants_plcs[i,2] = curr_plc
    }
  }
  variants_plcs = variants_plcs[order(variants_plcs$plc),]
  return(variants_plcs)
}

variants_from_infectivity_paper_plcs = find_variants_plcs(sorted_codons_table,
                                                          variants_from_infectivity_paper)
variants_of_interest_plcs = find_variants_plcs(sorted_codons_table,
                                               variants_of_interest)
variants_of_concern_plcs = find_variants_plcs(sorted_codons_table,
                                              variants_of_concern)

library("ggplot2")
library("gtable")
library("grid")
library("gridExtra")
library("tidyquant")
library("dplyr")
library("tidyverse")

df = data.frame(matrix(0,dim(sorted_codons_table)[1],2))
colnames(df) = c("state_num","variant_indicator") 
df$state_num = c(1:dim(sorted_codons_table)[1])
df$variant_indicator[variants_of_interest_plcs$plc] = 1
ggplot(df,aes(x = state_num, y = variant_indicator)) +
  # geom_line() + 
  geom_ma(ma_fun = SMA, n = 2000)


gisaid_variants_plcs = find_variants_plcs(sorted_codons_table,gisaid_variants_states)
new_gisaid_variants_plcs = find_variants_plcs(sorted_codons_table,new_gisaid_states)

amino_pairs = apply(sorted_codons_table[,c("amino","output_amino")],1,paste0,collapse = ".")
sort(table(amino_pairs[1:1000]))
sort(table(sorted_codons_table$amino[1:1000]))

amino_acids = unique(codons_table$amino_acid)
amino_acids_exposures = matrix(0,length(amino_acids),1)
names(amino_acids_exposures) = amino_acids
for (i in 1:length(amino_acids)){
  plcs = which(sorted_codons_table$amino==amino_acids[i])
  amino_acids_exposures[i] = sum(sorted_codons_table$exposure[plcs])
}
sort(amino_acids_exposures)


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
setwd(paste0(begin_path,"vars/prediction/all/"))
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

results = results2 = NULL
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
  plot_plcs = delta*c(1:99)
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
                      (tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)))[c(1,5,10)],
                      (tmp[plot_plcs]/ratios[[j]][plot_plcs,1])[c(1,5,10)])
    
    print(paste(curr_results))
    print(paste(curr_results2))
    results = rbind(results, curr_results)
    results2 = rbind(results2, curr_results2)
    # print(AUC)
  }
}
colnames(results) = c("model_num","type","AUC","AUC_null_all","AUC_null_time","Lift[1]","Lift[2]","Lift[3]","Lift[4]","Lift[5]"
                      ,"Lift[6]","Lift[7]","Diffs[1]","Diffs[2]","Diffs[3]","Diffs[4]","Diffs[5]"
                      ,"Diffs[6]","Diffs[7]")
colnames(results2) = c("model_num","type","AUC","AUC_null_all","AUC_null_time","Lift[1]","Lift[5]","Lift[10]","Lift_base[1]",
                       "Lift_base[5]","Lift_base[10]")
write.csv(results2,file = "results_spike.csv")


# Debug:
data = cbind(curr_output,curr_predictions[[5]][non_syn_plcs,3],curr_predictions[[6]][non_syn_plcs,3],curr_codons_table[non_syn_plcs,])
data = cbind(curr_output,curr_predictions[[4]][non_syn_plcs,3],curr_predictions[[6]][non_syn_plcs,3],curr_codons_table[non_syn_plcs,])
colnames(data)[1:3] = c("output","prediction","sub_model_num")
data = data[order(data$prediction,decreasing = TRUE),]
setwd(paste0(begin_path,"vars/prediction/syn_non_syn/"))
save(data,file = "non_syn.poisson.model_3")


###############3

get_plots = function(i,j,y_flag,sites_flag,divided_flag,amino_flag,spike_flag,states_flag){
  setwd(paste0(begin_path,"vars/prediction/syn_non_syn/"))
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
  plot_plcs = delta*c(1:99)
  
  if (spike_flag==1){
    title_name = "Spike gene:"
  }else{
    title_name = "All genes:"
  }
  
  tmp = ratios[[j]][,i]
  plot_data = data.frame(cbind(0.01*c(1:99),tmp[plot_plcs]/(length(which(curr_output>0))/length(curr_output)),
                               matrix(1,length(plot_plcs),1),tmp[plot_plcs]/ratios[[j]][plot_plcs,1]))
  colnames(plot_data) = c("Index","model","ones","null_lift")
  lift_plot = ggplot(plot_data,aes(x=Index,y=Lift))+
    geom_point(aes(y=model,col = "model"),size = 2)+
    geom_point(aes(y=null_lift,col = "null_lift"),size = 2)+
    geom_line(aes(y=ones),linetype = "longdash",size = 1.5,color = "black")+
    theme_classic()+ggtitle(paste(title_name,"Lift curve"))+xlab("Proportion of states")+ylab("Lift")+
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
tmp = get_plots(i=3,j=5,y_flag=0,sites_flag=0,divided_flag=1,amino_flag=1,spike_flag=1,states_flag=0)
spike_lift_plot = tmp[[1]]; spike_roc_plot = tmp[[2]]; spike_y_lim = tmp[[3]]
tmp = get_plots(i=3,j=5,y_flag=0,sites_flag=0,divided_flag=1,amino_flag=1,spike_flag=0,states_flag=0)
all_lift_plot = tmp[[1]]; all_roc_plot = tmp[[2]]; all_y_lim = tmp[[3]]
grid.arrange(all_lift_plot + ylim(0,max(all_y_lim,spike_y_lim)),spike_lift_plot + 
               ylim(0,max(all_y_lim,spike_y_lim)),ncol = 2)
grid.arrange(all_roc_plot,spike_roc_plot,ncol = 2)

setwd(paste0(begin_path,"figures/"))

#########
# Divide to different outputs:
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("codons_table_prediction_all")
load("syn_ti_const2")
load("syn_tv_const2")
load("non_syn_ti_const2")
load("non_syn_tv_const2")
all_predictions_flag=1
if (all_predictions_flag){
  input_cols = c("L.neighbor","R.neighbor","base","codon","site","ref_site","codon.pos","gene","mat_peptide",
                 "stem_loop","ref_seq","ref_codon","exposure","CG","amino_acid","amino_num")
}else{
  input_cols = c("L.neighbor","R.neighbor","base","codon","site","ref_site","codon.pos","gene","mat_peptide",
                 "stem_loop","ref_seq","ref_codon","exposure","CG","amino_acid","diff_seqs_num","amino_num")
}
output_cols = c("y","syn","non_syn","syn_exposure","non_syn_exposure","transitions","transversions")
bases = c("A","C","G","T")
num_states = dim(codons_table_prediction)[1]
codons_table_prediction_divided = matrix(0,3*num_states,25)
for (j in 1:num_states){
  print(j)
  curr_row = codons_table_prediction[j,]
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
  codons_table_prediction_divided[((3*j-2):(3*j)),] = unlist(tmp)
}
colnames(codons_table_prediction_divided) = colnames(tmp)
codons_table_prediction_divided = data.frame(codons_table_prediction_divided)

codons_table_prediction_divided = change_cols_to_numeric(numeric_cols = c("exposure",#"diff_seqs_num",
                                                                          "y","syn","non_syn","syn_exposure",
                                                                          "non_syn_exposure","transitions","transversions"),
                                                         data = codons_table_prediction_divided)
save(codons_table_prediction_divided,file = "codons_table_prediction_all_divided")
