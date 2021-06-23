begin_path = "/Users/keren/Dropbox/covid/new2/"
# begin_path = "C:/Users/Keren/Dropbox/covid/new2/"
# begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
begin_path_ncbi = paste0(begin_path,"ncbi_tree/")

diff_seqs_num_flag = 1
ncbi_flag =0
if (ncbi_flag){
  curr_begin_path = begin_path_ncbi
}else{
  curr_begin_path = begin_path
}

setwd(paste0(curr_begin_path,"vars/"))
load("codons_table")
load("codons_table_prediction_all")

# if (diff_seqs_num_flag==1){
#   codons_table_prediction = codons_table_prediction[which(codons_table_prediction$diff_seqs_num==0 |
#                                                             codons_table_prediction$diff_seqs_num>1),]
# }else if (diff_seqs_num_flag==0){
#   codons_table_prediction = codons_table_prediction[which(codons_table_prediction$diff_seqs_num==0 |
#                                                             codons_table_prediction$diff_seqs_num==1),]
# }
# # if diff_seq_num_flag==2 predict for all rows of codons_table_prediction


# Insert syn/non_syn exposure to the exposure:
multiply_exposures = function(curr_codons_table,col_name_to_multiply,col_name_to_make_1){
  new_codons_table = curr_codons_table
  plcs = which(curr_codons_table[,col_name_to_multiply]>0)
  new_codons_table[plcs,"exposure"] = curr_codons_table[plcs,"exposure"]*
    curr_codons_table[plcs,col_name_to_multiply]
  new_codons_table[plcs,col_name_to_make_1] = 1
  return(new_codons_table)
}
codons_table_prediction_syn = multiply_exposures(curr_codons_table = codons_table_prediction,
                                                 col_name_to_multiply = "syn_exposure",
                                                 col_name_to_make_1 = "syn_exposure")
codons_table_prediction_non_syn = multiply_exposures(curr_codons_table = codons_table_prediction,
                                                     col_name_to_multiply = "non_syn_exposure",
                                                     col_name_to_make_1 = "non_syn_exposure")

codons_table_syn = multiply_exposures(curr_codons_table = codons_table,
                                      col_name_to_multiply = "syn_exposure",
                                      col_name_to_make_1 = "syn_exposure")
codons_table_non_syn = multiply_exposures(curr_codons_table = codons_table,
                                          col_name_to_multiply = "non_syn_exposure",
                                          col_name_to_make_1 = "non_syn_exposure")

load("iterate_vals")
setwd(paste0(begin_path,"vars/models/"))
# load("unique_models.5")
# load("unique_models.34")
if (ncbi_flag){
  load("ncbi_models")
  unique_models = unique_models_ncbi
}else{
  load("Lanfear_models")
  unique_models = unique_models_Lanfear
}
setwd(paste0(curr_begin_path,"vars/models/"))
load("model_ids.12")
# load("model_ids.5")
# load("model_ids.34")
setwd(paste0(curr_begin_path,"vars/details/"))
load("details.12")
# load("details.5")
# load("details.34")

# null_model_num = which(unique_models$ID=="33333333333")
null_model_num = which(unique_models$ID=="33333333331")
# null_model_num = which(unique_models$ID=="33333333332")
# null_model_num = which(unique_models$ID=="22232322321")
null_time = null_model_num
model_nums = c(null_model_num,null_time,c(1:10))
# model_nums = c(null_model_num,null_time,c(454,6,12069,12181,12396,12,632,1325,347,933)) # Poisson winning models
num_models = length(model_nums)
tmp = data.frame(matrix(NA,dim(codons_table_prediction)[1],(num_models)))
row.names(tmp) = row.names(codons_table_prediction)
colnames(tmp)[1:2] = c("null_all_exposures","null_time_exposure")
colnames(tmp)[3:dim(tmp)[2]] = model_nums[3:length(model_nums)]
internal_tmp = data.frame(matrix(NA,dim(codons_table)[1],(num_models)))
row.names(internal_tmp) = row.names(codons_table)
predictions = internal_predictions = list()
predictions_names = c("nb.syn","p.syn","check.syn","nb.non_syn","p.non_syn","check.non_syn")
# predictions_names = c("nb.y","p.y","check.y")
for (i in 1:length(predictions_names)){
  predictions[[i]] = tmp
  internal_predictions[[i]] = internal_tmp
}
names(predictions) = predictions_names
names(internal_predictions) = predictions_names

codon.pos_pos = which(colnames(iterate_vals)=="codon.pos")
L.N.pos = which(colnames(iterate_vals)=="L.neighbor")
R.N.pos = which(colnames(iterate_vals)=="R.neighbor")
output_pos = dim(model_ids)[2]-1

setwd(curr_begin_path)
debug.nb.converged = debug.p.converged = NULL

curr_model_names = as.formula(paste("y~offset(log(exposure))"))
data = codons_table
source(paste0(begin_path,"debug_functions_apple.R"))
a = optimize(nb.glm.theta.data, data=codons_table,curr_model_names = curr_model_names,
             interval=c(0.05, 100), maximum=TRUE)
theta = as.numeric(a[1])

for (j in 1:length(model_nums)){
  print(j)
  curr_id = unique_models$ID[model_nums[j]] 
  unique_models[which(unique_models$ID==curr_id),]
  plcs = which(model_ids$ID==curr_id)
  if(substring(curr_id,first = codon.pos_pos,last = codon.pos_pos)=="1" &
     substring(curr_id,first = L.N.pos,last = R.N.pos)!="55"){
    # if the position of L.N,R.N changes change this accordingly
    curr_id2 = as.numeric(paste0("55",
                                 substring(curr_id,first = 3,last=output_pos),collapse = ""))
    plcs = c(plcs, which(model_ids$ID==curr_id2))
  }
  tmp = apply(details[plcs,],1,paste0,collapse=".")
  plcs = plcs[duplicated(tmp)==FALSE]
  details[plcs,]
  print(length(plcs))
  
  for (i in 1:length(plcs)){ 
    print(i)
    output_type = details$output[plcs[i]]
    if (output_type==1 | output_type==3){
      prediction_plcs = c(1:3)
      if (output_type==1){
        curr_codons_table = codons_table_syn
        curr_codons_table_prediction = codons_table_prediction_syn
        if (j==2){
          curr_codons_table = codons_table
          curr_codons_table_prediction = codons_table_prediction
          curr_model_names_null = "syn ~ offset(log(exposure))"
        }
      }
    }else if (output_type==2| output_type==4){
      prediction_plcs = c(4:6)
      if (output_type==2){
        curr_codons_table = codons_table_non_syn
        curr_codons_table_prediction = codons_table_prediction_non_syn
        if (j==2){
          curr_codons_table = codons_table
          curr_codons_table_prediction = codons_table_prediction
          curr_model_names_null = "non_syn ~ offset(log(exposure))"
        }
      }
    }else if (output_type==5){
      prediction_plcs = c(1:3)
    }
    tmp.nb = tmp.p = model.nb = model.p = NA
    # get model from original data
    data = curr_codons_table
    source(paste0(begin_path,"debug_functions_apple.R"))
    tmp = get_data(details[plcs[i],],iterate_vals)
    curr_data = tmp[[1]]; curr_details = tmp[[2]]
    tmp = prepare_data_rm_cols(curr_data)
    curr_data = tmp[[1]]; curr_model_names = tmp[[2]]
    rel_colnames = colnames(curr_data)
    if (j==2){
      curr_model_names = curr_model_names_null
    }
    model.nb = find_model.nb.data(curr_data,curr_model_names)
    model.p = find_model.P.data(curr_data, curr_model_names)
    if(is.list(model.nb)==FALSE & is.list(model.p)){
      model.nb = model.p
      debug.nb.converged = c(debug.nb.converged,paste0(c("model.",j,".plc.",i,"changed_to_poisson"),collapse = ""))
    }else if(is.list(model.p)==FALSE & is.list(model.nb)){
      model.nb = model.p
      debug.p.converged = c(debug.p.converged,paste0(c("model.",j,".plc.",i,"changed_to_NB"),collapse = ""))
    }
    if (is.list(model.nb) & is.list(model.p)){
      tmp.nb = predict(model.nb,curr_data)
      tmp.p = predict(model.p,curr_data)
      
      if(details[plcs[i],"output_sum"]==0){
        internal_predictions[[prediction_plcs[1]]][names(tmp.nb),j] = 0
        internal_predictions[[prediction_plcs[2]]][names(tmp.p),j] = 0
      }else{
        internal_predictions[[prediction_plcs[1]]][names(tmp.nb),j] = exp(tmp.nb)
        internal_predictions[[prediction_plcs[2]]][names(tmp.p),j] = exp(tmp.p)
      }
      if (sum(is.na(internal_predictions[[prediction_plcs[3]]][names(tmp.nb),j]))==length(tmp.nb)){
        internal_predictions[[prediction_plcs[3]]][names(tmp.nb),j] = i
      }else{
        print("2 models for the same row!")
        internal_predictions[[prediction_plcs[3]]][names(tmp.nb),j] = -1000*i
        # break
      }
      
      # predict for leaves
      data = curr_codons_table_prediction
      source(paste0(begin_path,"debug_functions_apple.R"))
      tmp = get_data(details[plcs[i],],iterate_vals)
      curr_data.predict = tmp[[1]]; curr_details = tmp[[2]]
      curr_data.predict = data[row.names(curr_data.predict),rel_colnames]
      if (sum(is.element(rel_colnames,colnames(curr_data.predict)))!=length(rel_colnames)){
        print("Bug! there are missing cols in the predict data")
        break
      }
      check_prediction = tryCatch({predict(model.nb,curr_data.predict)},
                                  error=function(cond){})
      if (is.null(check_prediction)){
        print("new categorical vars")
        tmp.nb = tmp.p = 0
        # Handle missing levels in the original data: find the column with missing levels,
        # first predict for all rows except for the missing levels rows. Then create a new model
        # without this explaining variable and predict for the missing levels rows.
        ks = missing_val_rows =  NULL
        for (k in 1:(dim(curr_data.predict)[2]-2)){
          unique_vals_predict = unique(curr_data.predict[,k])
          unique_vals_orig = unique(curr_data[,k])
          if (sum(is.element(unique_vals_predict,unique_vals_orig))!=length(unique_vals_predict)){
            ks = c(ks,k)
            missing_val = unique_vals_predict[which(is.element(unique_vals_predict,unique_vals_orig)==FALSE)]
            missing_val_rows = c(missing_val_rows,
                                 which(is.element(curr_data.predict[,k],missing_val)))
          }
        }
        missing_val_rows = unique(missing_val_rows)
        missing_val_data = curr_data.predict[missing_val_rows,-ks]
        missing_prediction.nb = predict(model.nb,curr_data.predict[-missing_val_rows,]) 
        missing_prediction.p = predict(model.p,curr_data.predict[-missing_val_rows,])
        curr_data_missing_col = curr_data[,-ks]
        rel_colnames = colnames(curr_data_missing_col)
        model.nb = find_model.nb.data(curr_data_missing_col,curr_model_names)
        model.p = find_model.P.data(curr_data_missing_col, curr_model_names)
        tmp.nb = c(missing_prediction.nb,predict(model.nb,missing_val_data))
        tmp.p = c(missing_prediction.p,predict(model.p,missing_val_data))
      }else{
        tmp.nb = predict(model.nb,curr_data.predict)
        tmp.p = predict(model.p,curr_data.predict)
      }
      
      if (length(tmp.nb)>0){
        if(details[plcs[i],"output_sum"]==0){
          predictions[[prediction_plcs[1]]][names(tmp.nb),j] = 0
          predictions[[prediction_plcs[2]]][names(tmp.p),j] = 0 
        }else{
          predictions[[prediction_plcs[1]]][names(tmp.nb),j] = exp(tmp.nb)
          predictions[[prediction_plcs[2]]][names(tmp.p),j] = exp(tmp.p)
        }
        if (is.na(predictions[[prediction_plcs[3]]][names(tmp.nb),j])){
          predictions[[prediction_plcs[3]]][names(tmp.nb),j] = i
        }else{
          print("2 models for the same row!")
          predictions[[prediction_plcs[3]]][names(tmp.nb),j] = -1000*i
          # break
        }
      }
      if (dim(predictions[[1]])[1]>dim(codons_table_prediction)[1]){
        print(dim(predictions[[1]]))
        break
      }
    }else{
      if (model.nb==0){
        debug.nb.converged = c(debug.nb.converged,paste0(c("model.",j,".plc.",i),collapse = ""))
      }else{
        debug.p.converged = c(debug.p.converged,paste0(c("model.",j,".plc.",i),collapse = ""))
      } 
    }
  }
}

######################
setwd(paste0(curr_begin_path,"vars/prediction/all/"))
table(codons_table_prediction$diff_seqs_num)
unique_models[1,]
save(predictions,file = "no_sites_normalization.predictions_syn_non_syn")
save(internal_predictions,file = "no_sites_normalization.internal_predictions_syn_non_syn")
save(debug.nb.converged,file = "no_sites_normalization.debug.nb.converged.syn_non_syn")
save(debug.p.converged,file = "no_sites_normalization.debug.p.converged.syn_non_syn")




