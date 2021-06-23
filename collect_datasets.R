begin_path = "C:/Users/ifog/Dropbox/covid/new2/ncbi_tree/"
# begin_path = "/Users/keren/Dropbox/covid/new2/"
# begin_path = "C:/Users/Keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/datasets/datasets_12"))


##################################3
# Get datasets, divide to syn/non_syn, transition/transversion, y, A/C/G/T
file_names = dir()
count = 0
count2 = 0
output_num = 1
dataset = data.frame(matrix("-",3543468, 22))
start_time = Sys.time()
for (i in 1:length(file_names)){
  print(i)
  if (paste0(strsplit(file_names[i],split = "")[[1]][1:14]
             ,collapse = "")=="short_details."){
    load(file_names[i])
    count2 = count2+1
    print(count2)
    if (is.null(short_details)==FALSE){
      output_plcs = which(short_details$output==output_num)
      if (length(output_plcs)>0){
        rel_short_details = short_details[output_plcs,]
        dataset[(count+1):(count+dim(rel_short_details)[1]),] = rel_short_details
        count = count + dim(rel_short_details)[1]
        rm(short_details)
        rm(rel_short_details)
      }
    }
  }
}
print(Sys.time()-start_time)
print(dim(dataset))
setwd(paste0(begin_path,"vars/datasets/datasets_12"))
load("short_details.parallel_num.1006_q.66_w.18_b3_j5")
colnames(dataset) = colnames(short_details)
save(dataset,file = paste0("dataset.",output_num))
dataset_reg_name = dataset$reg_name2
save(dataset_reg_name, file = paste0("dataset_reg_name.",output_num))
plcs = which(duplicated(dataset$reg_name2)==FALSE)
dataset = dataset[plcs,]
save(dataset,file = paste0("unique_dataset.",output_num))

########
# Check:
load("dataset_syn_non_syn")
tmp = apply(dataset[,1:11],1,paste0,collapse = ".")
tmp2 = which(duplicated(tmp)==FALSE)

load("dataset_ti_tv")
tmp = apply(dataset[,1:11],1,paste0,collapse = ".")
tmp2 = which(duplicated(tmp)==FALSE)


load("dataset_ACGT")
tmp = apply(dataset[,1:11],1,paste0,collapse = ".")
tmp2 = which(duplicated(tmp)==FALSE)

#
load("unique_dataset_1")
all_dataset = dataset
load("unique_dataset_2")
all_dataset = rbind(all_dataset,dataset)
load("unique_dataset_3")
all_dataset = rbind(all_dataset,dataset)
load("unique_dataset_4")
all_dataset = rbind(all_dataset,dataset)
load("unique_dataset_5")
all_dataset = rbind(all_dataset,dataset)
load("unique_dataset_6")
all_dataset = rbind(all_dataset,dataset)
dataset = all_dataset
# tmp = apply(dataset,1,paste0,collapse = "")
# tmp2 = which(duplicated(tmp)==FALSE)
plcs = which(duplicated(dataset$reg_name2)==FALSE)
dataset = dataset[plcs,]
save(dataset,file = "unique_dataset_all")

dataset_reg_name_all = NULL
for (i in 1:6){
  print(i)
  load(paste0("dataset_reg_name_",i,collapse = ""))
  dataset_reg_name_all = c(dataset_reg_name_all,dataset_reg_name)
}
save(dataset_reg_name_all,file = "dataset_reg_name_all")

######3
# Check:
tmp_all = NULL
for(i in 1:6){
  print(i)
  load(paste0("dataset_",i,collapse = ""))
  tmp = apply(dataset[,-12],1,paste0,collapse = "")
  tmp2 = which(duplicated(tmp)==FALSE)
  # tmp2 = which(duplicated(tmp)==TRUE)
  if (length(tmp)!=length(tmp2)){
    print(paste("dataset ",i,"contains duplicates!"))
  }
  rm(dataset)
  tmp_all = c(tmp_all,tmp)
}

tmp2 = which(duplicated(tmp_all)==TRUE)
load("dataset_reg_name_all")
count = 0
for (i in 30112138:30114675){#length(dataset_reg_name)){
  if (dataset_reg_name[[i]]!="-"){
    print(i)
    count = count+1
  }
}

# ###############
# j=9
# setwd("/Users/keren/Dropbox/covid/vars/datasets/base_all/")
# file_names = dir()
# count = 0
# dataset = NULL
# count2 = 0
# for (i in 1:length(file_names)){
#   print(i)
#   if (paste0(strsplit(file_names[i],split = "")[[1]][1:14]
#              ,collapse = "")=="short_details."){
#     load(file_names[i])
#     if (is.null(short_details)==FALSE){
#       count2 = count2+1
#       print(count2)
#       count = count + dim(short_details)[1]
#       plcs = which(short_details$output==j)
#       if (length(plcs)>0){
#         dataset = rbind(dataset,short_details[plcs,])
#       }
#       rm(short_details)
#     }
#   }
# }
# print(Sys.time()-start_time)
# print(dim(dataset))
# setwd("/Users/keren/Dropbox/covid/vars/datasets/")
# save(dataset,file = paste0("dataset_base_all_i.",j,collapse = ""))
# dataset_reg_name = dataset$reg_name2
# save(dataset_reg_name, file = 
#        paste0("dataset_reg_name_base_all_i.",j,collsape = ""))
# plcs = which(duplicated(dataset$reg_name2)==FALSE)
# dataset = dataset[plcs,]
# save(dataset,file = paste0("unique_dataset_base_all_i.",j,collapse = ""))
# 
# 
# ###########3
# # Merge all unique_datasets into one file so that regressions will be found
# # for all together
# setwd("/Users/keren/Dropbox/covid/vars/datasets/")
# dataset2 = NULL
# for (i in 3:9){
#   load(paste0("unique_dataset_base_all_i.",i,collapse = ""))
#   dataset2 = rbind(dataset2,dataset)
# }
# plcs = which(duplicated(dataset2$reg_name2)==FALSE)
# dataset2 = dataset2[plcs,]
# dataset = dataset2
# save (dataset,file = "unique_dataset_base_all_i_3_9")