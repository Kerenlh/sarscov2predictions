begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/site_patterns/all2/"))
load("all_site_patterns_plcs")
load("all_no_subs_plcs")
load("all_site_patterns")
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
length_seqs = length(aligned_seqs[1][[1]])
num_seqs = length(aligned_seqs)
init_seq_no_subs_plcs = strsplit(as.character(aligned_seqs[[1]]),split ="")[[1]][no_subs_plcs]
rm(aligned_seqs)


sort_file_names = function(file_names){
  tmp = unlist(strsplit(file_names,split = "\\."))
  nums = as.numeric(tmp[seq(2,length(tmp),2)])
  file_names = file_names[order(nums)]
  return(file_names)
}

setwd("/Users/keren/Desktop/covid_files/new2/Lanfear trees/")
file_names = dir()
load(file_names[1])
nums = NULL
for (i in 1:length(file_names)){
  if (strsplit(file_names[i],split = "[.]")[[1]][1]=="tree_list"){
    nums = c(nums,strsplit(file_names[i],split = "[.]")[[1]][2])
  }
}
nums = as.numeric(nums)
for (j in 1:30){
  setwd("/Users/keren/Desktop/covid_files/new2/Lanfear trees/")
  start_time = Sys.time()
  if (j==30){
    part_seqs = data.frame(matrix("M",length(bases),(length_seqs-29000)))
  }else{
    part_seqs = data.frame(matrix("M",length(bases),1000))
  }
  rownames(part_seqs) = names(bases)
  rel_plcs = c(((j-1)*1000+1):(j*1000))
  rel_file_plcs = which(is.element(nums,rel_plcs))
  for (i in 1:length(rel_file_plcs)){
    print(i)
    load(file_names[rel_file_plcs[i]])
    curr_num = as.numeric(strsplit(file_names[rel_file_plcs[i]],split = "[.]")[[1]][2])
    part_seqs[names(bases),(curr_num-((j-1)*1000))] = bases
  }
  rel_no_subs_plcs = which(is.element(no_subs_plcs,rel_plcs))
  if (length(rel_no_subs_plcs)>0){
    for(i in 1:length(rel_no_subs_plcs)){
      print(i)
      part_seqs[,(no_subs_plcs[rel_no_subs_plcs[i]]-((j-1)*1000))] = 
        init_seq_no_subs_plcs[rel_no_subs_plcs[i]]
    }
  }
  setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  save(part_seqs,file = paste0("part_seqs.",j))
  print(Sys.time()-start_time)
}

#####
# Add duplicated site patterns:
# setwd(paste0(begin_path,"vars/site_patterns/all2/"))
# load("all_site_patterns")
# load("all_site_patterns_plcs")
# unique_site_patterns = site_patterns[which(duplicated(site_patterns)==FALSE)]
# site_patterns_numbers = site_patterns_plcs[which(duplicated(site_patterns)==FALSE)]
# site_patterns_plcs = cbind(site_patterns_plcs,matrix(0,length(site_patterns_plcs),1))
# colnames(site_patterns_plcs) = c("site","tree_list_num")
# site_patterns_plcs = data.frame(site_patterns_plcs)
# for (i in 1:length(unique_site_patterns)){
#   print(i)
#   site_patterns_plcs$tree_list_num[which(site_patterns==unique_site_patterns[i])] =
#     site_patterns_numbers[i]
# }
# save(site_patterns_plcs,file = "site_patterns_plcs_matrix")
setwd(paste0(begin_path,"vars/site_patterns/all2/"))
load("site_patterns_plcs_matrix")

setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
file_names = dir()
file_names = file_names[-1] #rm "all_seqs"
file_names = sort_file_names(file_names)
print(file_names) #make sure it's sorted!
for (i in 1:15){
  setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  print(file_names[i])
  start_time = Sys.time()
  load(file_names[i])
  missing_plcs_local = which(part_seqs[1,]=="M")
  missing_plcs = missing_plcs_local +(i-1)*1000
  setwd("/Users/keren/Desktop/covid_files/new2/Lanfear trees/")
  print(length(missing_plcs))
  for (j in 1:length(missing_plcs)){
    rel_tree_list_num = 
      site_patterns_plcs$tree_list_num[which(site_patterns_plcs$site==missing_plcs[j])]
    tree_name = paste0("tree_list.",rel_tree_list_num)
    if(file.exists(tree_name)){
      load(tree_name)
    }else{
      load(paste0(tree_name,".filepart"))
    }
    part_seqs[names(bases),missing_plcs_local[j]] = bases
  }
  print(Sys.time()-start_time)
  print(length(which(part_seqs[1,]=="M")))
  setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  save(part_seqs,file = paste0("P",file_names[i]))
}

setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
file_names = dir()
file_names = file_names[32:length(file_names)]
file_names = sort_file_names(file_names)
all_seqs = NULL
for (i in 1:length(file_names)){
  print(i)
  load(file_names[i])
  tmp = apply(part_seqs,1,paste0,collapse = "")
  rm(part_seqs)
  all_seqs = apply(cbind(all_seqs,tmp),1,paste0,collapse = "")
}
save(all_seqs,file = "all_seqs")

