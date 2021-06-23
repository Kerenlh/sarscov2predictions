cluster = 1
if (cluster){
  # begin_path = "/a/home/cc/math/kerenlh/covid_new/"
  # begin_path = "/Users/keren/Dropbox/covid/new2/"
  begin_path = "C:/Users/ifog/Dropbox/covid/new2/ncbi_tree/"
  # begin_path = "C:/Users/Keren/Dropbox/covid/new2/"
}else{
  begin_path = "/a/home/cc/math/saharon/keren/"
}

wait = runif(1, 0, 15)
Sys.sleep(wait)
tmp = as.numeric(Sys.time())
set.seed(tmp)
setwd(paste0(begin_path,"vars/datasets/datasets_12"))
library(hash)
load("dataset_reg_name.1")
reg_hash_names = unique(dataset_reg_name)
print(length(reg_hash_names))
iter_num = 500
total_num_files = round(length(reg_hash_names)/iter_num)+1

setwd(paste0(begin_path,"vars/reg_hash_plcs/dataset_1/"))
if (file.exists("in_progress_files")){
  load("in_progress_files")
}else{
  in_progress_files = NULL
  save(in_progress_files,file = "in_progress_files")
}
file_names = dir()
nums = NULL
for (i in 1:length(file_names)){
  if (strsplit(file_names[i],split = "[.]")[[1]][1]=="reg_hashes_plcs"){
    nums = c(nums,strsplit(file_names[i],split = "[.]")[[1]][2])
  }
}
nums = as.numeric(nums)
nums = c(nums,in_progress_files)
missing_nums = which(is.element(1:total_num_files,nums)==FALSE)
# print("missing_files:")
# print(missing_nums)
print("in_progress:")
print(in_progress_files)
print(paste("length(missing_nums):",length(missing_nums)))
reg_hash_num = missing_nums[sample(1:length(missing_nums),1)]

print(reg_hash_num)
in_progress_files = c(in_progress_files,reg_hash_num)
save(in_progress_files,file = "in_progress_files")

for (k in 1:length(missing_nums)){
  reg_hash_num = missing_nums[k]
  print(k)
  print(reg_hash_num)
  reg_hashes_plcs = hash()
  for (j in ((reg_hash_num-1)*iter_num+1):min((reg_hash_num*iter_num),length(reg_hash_names))){
    # print(j)
    plcs = which(dataset_reg_name==reg_hash_names[j])
    reg_hashes_plcs[[ reg_hash_names[j] ]] = plcs
  }
  file_name = paste0("reg_hashes_plcs.",reg_hash_num,collapse = "")
  save(reg_hashes_plcs,file = file_name)
}
load("in_progress_files")
in_progress_files = in_progress_files[-which(in_progress_files==reg_hash_num)]
save(in_progress_files,file = "in_progress_files")
count = 0
while(file.exists(file_name)==FALSE){
  count = count+1
  if(count>10){
    print("endless loop")
    print(reg_hash_num)
    break
  }
  print ("ERROR!!! didn't save file!")
  print(reg_hash_num)
  print(paste("count  ",count))
  save(reg_hashes_plcs,file = file_name)
  if (file.exists(file_name)){
    print("Good, saved!")
  }else{
    print("Not saved!")
  }
}

wait = 15
Sys.sleep(wait)
stop("so that the workspace is not saved")



