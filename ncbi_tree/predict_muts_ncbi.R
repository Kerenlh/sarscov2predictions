begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
# begin_path = "C:/Users/Keren/Dropbox/covid/new2/ncbi_tree/"
# begin_path = "C:/Users/ifog/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("codons_table_prediction")

load("codons_table")
load("iterate_vals")
setwd(paste0(begin_path,"vars/models/"))
load("Lanfear_models")
unique_models = unique_models_Lanfear
load("model_ids.12")
# load("model_ids.5")
setwd(paste0(begin_path,"vars/details/"))
load("details.12")
# load("details.5")

predictions.nb.syn = predictions.nb.non_syn = 
  predictions.p.syn =predictions.p.non_syn = data.frame(matrix(NA,dim(codons_table_prediction)[1],100))
internal_predictions.nb.syn = internal_predictions.nb.non_syn = 
  internal_predictions.p.syn = internal_predictions.p.non_syn = data.frame(matrix(NA,dim(codons_table_prediction)[1],100))
row.names(predictions.nb.syn) = row.names(predictions.nb.non_syn) = 
  row.names(predictions.p.syn) = row.names(predictions.p.non_syn) = row.names(codons_table_prediction)
row.names(internal_predictions.nb.syn) = row.names(internal_predictions.nb.non_syn) = 
  row.names(internal_predictions.p.syn) = row.names(internal_predictions.p.non_syn) = row.names(codons_table_prediction)

for (j in 1:100){
  print(j)
  curr_ID = unique_models$ID[j] #19
  # curr_ID = "33333333333"
  unique_models[which(unique_models$ID==curr_ID),]
  plcs = which(model_ids$ID==curr_ID)
  details[plcs,]
  # details$output[plcs] = "1"
  setwd(begin_path)
  for (i in 1:length(plcs)){  
    # print(i)
    # get model from original data
    data = codons_table
    source(paste0(begin_path,"debug_functions_apple.R"))
    tmp = get_data(details[plcs[i],],iterate_vals)
    curr_data = tmp[[1]]; curr_details = tmp[[2]]
    # print(curr_data)
    # print(curr_details)
    # print(details[plcs[i],])
    tmp = prepare_data_rm_cols(curr_data)
    curr_data = tmp[[1]]; curr_model_names = tmp[[2]]
    print(dim(curr_data))
    rel_colnames = colnames(curr_data)
    model.nb = find_model.nb.data(curr_data,curr_model_names)
    model.p = find_model.P.data(curr_data, curr_model_names)
    tmp = predict(model.nb,curr_data)
    # print(tmp)
    if (details$output[plcs[i]]==1){
      internal_predictions.nb.syn[names(tmp),j] = tmp
    }else if ((details$output[plcs[i]]==2)){
      internal_predictions.nb.non_syn[names(tmp),j] = tmp
    }
    tmp = predict(model.p,curr_data)
    if (details$output[plcs[i]]==1){
      internal_predictions.p.syn[names(tmp),j] = tmp
    }else if ((details$output[plcs[i]]==2)){
      internal_predictions.p.non_syn[names(tmp),j] = tmp
    }
    
    # predict for leaves
    data = codons_table_prediction
    source(paste0(begin_path,"debug_functions_apple.R"))
    tmp = get_data(details[plcs[i],],iterate_vals)
    curr_data = tmp[[1]]; curr_details = tmp[[2]]
    curr_data = data[row.names(curr_data),rel_colnames]
    # tmp = prepare_data_rm_cols(curr_data)
    # curr_data = tmp[[1]]; curr_model_names = tmp[[2]]
    print(dim(curr_data))
    tmp = predict(model.nb,curr_data)
    if (details$output[plcs[i]]==1){
      predictions.nb.syn[names(tmp),j] = tmp
    }else if ((details$output[plcs[i]]==2)){
      predictions.nb.non_syn[names(tmp),j] = tmp
    }
    tmp = predict(model.p,curr_data)
    if (details$output[plcs[i]]==1){
      predictions.p.syn[names(tmp),j] = tmp
    }else if ((details$output[plcs[i]]==2)){
      predictions.p.non_syn[names(tmp),j] = tmp
    } 
  }
}
cor(exp(internal_predictions.nb[,"syn"]),codons_table$syn)
cor(exp(internal_predictions.nb[,"non_syn"]),codons_table$non_syn)
plot(exp(internal_predictions.nb[,"non_syn"]),codons_table$non_syn)
cor((exp(internal_predictions.nb[,"syn"])+exp(internal_predictions.nb[,"non_syn"])),codons_table$y)
plot((exp(internal_predictions.nb[,"syn"])+exp(internal_predictions.nb[,"non_syn"])),codons_table$y)

cor(exp(predictions.nb[,"syn"]),codons_table_prediction$syn)
cor(exp(predictions.nb[,"non_syn"]),codons_table_prediction$non_syn)
plot(codons_table_prediction$syn,exp(predictions.nb[,"syn"]))
plot(exp(predictions.nb[,"non_syn"]),codons_table_prediction$non_syn)
cor((exp(predictions.nb[,"syn"])+exp(predictions.nb[,"non_syn"])),codons_table_prediction$y)
plot((exp(predictions.nb[,"syn"])+exp(predictions.nb[,"non_syn"])),codons_table_prediction$y)

syn_sites = unique(codons_table_prediction$ref_site[which(codons_table_prediction$syn>0)])
non_syn_sites = unique(codons_table_prediction$ref_site[which(codons_table_prediction$non_syn>0)])
tmp = unique(codons_table_prediction$ref_site[which(codons_table_prediction$syn==0 & codons_table_prediction$non_syn==0)])

predictions.nb = cbind(predictions.nb,codons_table_prediction$ref_site)
tmp = predictions.nb[order(predictions.nb$syn,decreasing = TRUE),]
table(is.element(tmp$`codons_table_prediction$ref_site`[1:2000],syn_sites))
tmp = predictions.nb[order(predictions.nb$non_syn,decreasing = TRUE),]
table(is.element(tmp$`codons_table_prediction$ref_site`[1:2000],non_syn_sites))
