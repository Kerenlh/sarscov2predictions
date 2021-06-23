begin_path = "/Users/keren/Dropbox/covid/new2/"
# begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
# begin_path = "C:/Users/Keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/details/"))

# all_details = NULL
# for (i in 1:2){
#   load(paste0("details.",i))
#   print(dim(details))
#   all_details = rbind(all_details,details)
# }
# details = all_details
# save(details,file = "details.12")
# rm(all_details)

# load("details.6789")
# load("details.12")
load("details.2")
setwd(paste0(begin_path,"vars/"))
load("iterate_vals")
load("codons_table")
data = codons_table
data = data[which(data$non_syn_exposure>0),]
# syn_non_syn_rows_num = length(which(data$syn_exposure>0))+length(which(data$non_syn_exposure>0))
syn_non_syn_rows_num = length(which(data$non_syn_exposure>0))
codon.pos_pos = which(colnames(iterate_vals)=="codon.pos")
L.N.pos = which(colnames(iterate_vals)=="L.neighbor")
R.N.pos = which(colnames(iterate_vals)=="R.neighbor")
# source(paste0(begin_path,"debug_functions_apple.R"))

##############
# Functions: #
##############
get_model_ids = function(details,iterate_vals){
  colnames_model_ids = c(colnames(iterate_vals),"output","ID")
  model_ids = matrix(1,dim(details)[1],length(colnames_model_ids))
  colnames(model_ids) = colnames_model_ids
  
  model_ids[which(details$codon.pos==4),"codon.pos"] = 2
  model_ids[which(is.na(details$codon.pos)),"codon.pos"] = 3
  model_ids[,c("R.neighbor","L.neighbor")] = 0
  # 1 = divide according to neighbors with respect to codon.pos
  model_ids[which(details$codon.pos==1 & details$L.neighbor<5 & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 1 
  model_ids[which(details$codon.pos==3 & is.na(details$L.neighbor) & details$R.neighbor<5),c("R.neighbor","L.neighbor")] = 1
  # 2 = have neighbors as an explaining variable outside the codon 
  model_ids[which(details$codon.pos==1 & details$L.neighbor==5 & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 2 
  model_ids[which(details$codon.pos==3 & is.na(details$L.neighbor) & details$R.neighbor==5),c("R.neighbor","L.neighbor")] = 2
  # 3 = don't include neighbors
  model_ids[which(details$codon.pos==1 & is.na(details$L.neighbor) & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 3 
  model_ids[which(details$codon.pos==3 & is.na(details$L.neighbor) & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 3
  # codon.pos==2 is the same for all options and remaining codon.pos are foreign for different options
  model_ids[which(details$codon.pos==2 & is.na(details$L.neighbor) & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 5
  
  # codon.pos==4:
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$L.neighbor<5),c("L.neighbor")] = 1
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$L.neighbor==5),c("L.neighbor")] = 2
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & is.na(details$L.neighbor)),c("L.neighbor")] = 3
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$R.neighbor<5),c("R.neighbor")] = 1
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$R.neighbor==5),c("R.neighbor")] = 2
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & is.na(details$R.neighbor)),c("R.neighbor")] = 3
  
  
  model_ids[which(details$mat_peptide==3),"mat_peptide"] = 2
  model_ids[which(is.na(details$mat_peptide)),"mat_peptide"] = 3
  model_ids[which(details$CG==4),"CG"] = 2
  model_ids[which(is.na(details$CG)),"CG"] = 3
  model_ids[which(details$codon==65),"codon"] = 2
  model_ids[which(is.na(details$codon)),"codon"] = 3
  model_ids[which(details$amino_acid==22),"amino_acid"] = 2
  model_ids[which(is.na(details$amino_acid)),"amino_acid"] = 3
  model_ids[which(details$gene==12),"gene"] = 2
  model_ids[which(is.na(details$gene)),"gene"] = 3
  model_ids[which(details$stem_loop==4),"stem_loop"] = 2
  model_ids[which(is.na(details$stem_loop)),"stem_loop"] = 3
  # model_ids[which(details$output==2),"output"] = 1 # Take only non-syn models
  # model_ids[which(details$output==1),"output"] = 5
  model_ids[which(details$output==1 | details$output==2),"output"] = 1
  model_ids[which(details$output==3 | details$output==4),"output"] = 2
  model_ids[which(details$output==5),"output"] = 3
  model_ids[which(details$output==6 | details$output==7 | 
                    details$output==8 | details$output==9),"output"] = 4
  
  
  model_ids[which(details$base==5),"base"] = 2
  model_ids[is.na(details$base),"base"] = 3
  
  model_ids[,"ID"] = as.numeric(apply(model_ids[,1:(dim(model_ids)[2]-1)],1,paste0,collapse = ""))
  model_ids = data.frame(model_ids)
  details_ids = apply(details[,1:(dim(model_ids)[2]-1)],1,paste0,collapse = "")
  tmp = list(model_ids,details_ids)
  return(tmp)
}

find_id_plcs = function(model_ids,unique_models){
  output_pos = dim(model_ids)[2]-1
  ID_pos = dim(model_ids)[2]
  
  for (i in 1:dim(unique_models)[1]){
    print(i)
    curr_id = unique_models$ID[i]
    plcs = which(model_ids$ID==curr_id)
    if(substring(curr_id,first = codon.pos_pos,last = codon.pos_pos)=="1" &
       substring(curr_id,first = L.N.pos,last = R.N.pos)=="55"){next}
    if(substring(curr_id,first = codon.pos_pos,last = codon.pos_pos)=="1" &
       substring(curr_id,first = L.N.pos,last = R.N.pos)!="55"){
      # if the position of L.N,R.N changes change this accordingly
      curr_id2 = as.numeric(paste0("55",
                                   substring(curr_id,first = 3,last=output_pos),collapse = ""))
      plcs = c(plcs, which(model_ids$ID==curr_id2))
    }
    tmp = apply(details[plcs,],1,paste0,collapse=".")
    plcs = plcs[duplicated(tmp)==FALSE]
    
    unique_models[i,colnames(model_ids)] = model_ids[plcs[1],]
    unique_models$num_models[i] = length(plcs)
    unique_models$sum_rows[i] = sum(details[plcs,"rows_num"])
    unique_models$logLik.NB[i] = sum(details[plcs,"logLik"])
    unique_models$logLik.P[i] = sum(details[plcs,"P.logLik"])
    
    if (substring(curr_id,first = output_pos,last = output_pos)==3){
      unique_models$num_missing_rows[i] = dim(data)[1]-unique_models$sum_rows[i]
    }
    if (substring(curr_id,first = output_pos,last = output_pos)==2){
      unique_models$num_missing_rows[i] = 2*dim(data)[1]-unique_models$sum_rows[i]
    }
    if (substring(curr_id,first = output_pos,last = output_pos)==1){ # syn/non_syn together have 27144 rows
      unique_models$num_missing_rows[i] = syn_non_syn_rows_num-unique_models$sum_rows[i]
    }
    if (substring(curr_id,first = output_pos,last = output_pos)==4){ # ACGT
      #!!!!!!!!!!! change back to 3*dim(data)[1]!!!!!!!!!!!!!!!!!!!1
      unique_models$num_missing_rows[i] = 3*dim(data)[1]-unique_models$sum_rows[i]
    }
    
    unique_models$total_num_models[i] = unique_models$num_models[i]+unique_models$num_missing_rows[i]
    
    unique_models$num_missing_models[i] =unique_models$num_missing_rows[i] # assumes all missing models include only 1 row!
    
    unique_models$df.NB[i] = sum(details$df[plcs])+unique_models$num_missing_models[i]
    unique_models$df.P[i] = sum(details$P.df[plcs])+unique_models$num_missing_models[i]
    
    unique_models$AIC.NB[i] = 2*(unique_models$df.NB[i]-unique_models$logLik.NB[i])
    unique_models$AIC.P[i] = 2*(unique_models$df.P[i]-unique_models$logLik.P[i])
    unique_models$NA.NB[i] = sum(details$added_log_lik[plcs])
    unique_models$NA.P[i] = sum(details$added_log_lik.P[plcs])
    
    if (length(which(details$added_log_lik[plcs]==1))>0){
      good_plcs.NB = plcs[-which(details$added_log_lik[plcs]==1)]
    }else{
      good_plcs.NB = plcs
    }
    unique_models[i,26:31] = summary(details[good_plcs.NB,"theta"])
    unique_models[i,32] = weighted.mean(details[good_plcs.NB,"theta"],
                                        details[good_plcs.NB,"rows_num"]/sum(details[good_plcs.NB,"rows_num"]))
    unique_models[i,33] = sum(details$output_0_or_1_row[plcs])
  }
  return(unique_models)
}

########

tmp = get_model_ids(details,iterate_vals)
model_ids = tmp[[1]]; model_details_ids = tmp[[2]]
tmp = unique(model_ids[,"ID"])

colnames_unique_models = c(colnames(model_ids),"num_models","logLik.NB","df.NB",
                           "num_missing_models","AIC.NB","NA.NB","logLik.P","df.P",
                           "AIC.P","NA.P","sum_rows","total_num_models","num_missing_rows",
                           "theta_min","theta_1st Qu.","theta_Median","theta_Mean","theta3rd Qu.",
                           "theta_Max","theta_weighted_mean","output_0_or_1_row_models")
unique_models = matrix(0,length(tmp),length(colnames_unique_models))
colnames(unique_models) = colnames_unique_models
unique_models = data.frame(unique_models)
unique_models$ID = tmp


unique_models = find_id_plcs(model_ids,unique_models)
unique_models = unique_models[-which(unique_models$codon.pos==0),]
setwd(paste0(begin_path,"vars/models/"))
unique_models = unique_models[order(apply(unique_models[,c("AIC.NB","AIC.P")],1,min)),]
save(unique_models,file = "unique_models.non_syn")
save(model_ids,file = "model_ids.non_syn")
