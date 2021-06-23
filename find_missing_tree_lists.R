# Check existing and missing tree_list files:
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/trees/"))
setwd("/Users/keren/Desktop/covid_files/new2/Lanfear trees/")
file_names = dir()
for (i in 1:length(file_names)){
  tryCatch(load(file_names[i]),error=function(cond){
    print(file_names[i])
    file.remove(file_names[i])})
}
file_names = dir()
tmp = strsplit(file_names,"\\.")
nums = NULL
for (i in 1:length(tmp)){
  nums = c(nums,as.numeric(tmp[[i]][2]))
}
cluster_nums = unique(nums)
save(cluster_nums,file = "cluster_nums")


begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("site_to_alignment")
setwd(paste0(begin_path,"vars/muts/"))
file_names = dir()
for (i in 1:length(file_names)){
  tryCatch(load(file_names[i]),error=function(cond){file.remove(file_names[i])})
}
file_names = dir()
for (i in 1:length(file_names)){
  tryCatch(load(file_names[i]),error=function(cond){print(i)})
}

tmp = strsplit(file_names,"\\.")
num_muts = NULL
for (i in 1:length(tmp)){
  print(i)
  num_muts = c(num_muts,as.numeric(tmp[[i]][3]))
}
nums = NULL
insertions = deletions = subs = matrix(0,sum(num_muts),6)
insertions_count = deletions_count = subs_count = 0

update_matrix = function(mat,plcs,count,muts){
  mat[((count+1):(count+length(plcs))),] = as.matrix(muts[plcs,1:6])
  count = count+length(plcs)
  return(list(count,mat))
}

for (i in 10454:length(tmp)){
  print(i)
  nums = c(nums,as.numeric(tmp[[i]][2]))
  load(file_names[i])
  alignment_site = matrix(as.numeric(tmp[[i]][2]),dim(muts)[1],1)
  site = matrix(site_to_alignment$site[as.numeric(tmp[[i]][2])],dim(muts)[1],1)
  muts = cbind(muts,alignment_site,site)
  
  # for (i in 1:dim(muts)[1]){
  #   if (children_bases[i,muts$before[i]]>0){
  #     problematic_muts = c(problematic_muts,i)
  #   }
  # }
  problematic_muts = which(is.element(muts$parent,muts$base))
  problematic_muts = c(problematic_muts,which(is.element(muts$base,muts$parent)))
  
  
  insertions_plcs = which(muts$before=="-")
  if (length(insertions_plcs)>0){
    outputs = update_matrix(mat = insertions,plcs = insertions_plcs,count = insertions_count,muts = muts)
    insertions_count = outputs[[1]]; insetions = outputs[[2]]
  }
  
  deletions_plcs = which(muts$after=="-")
  if(length(deletions_plcs)>0){
    outputs = update_matrix(mat = deletions,plcs = deletions_plcs,count = deletions_count,muts = muts)
    deletions_count = outputs[[1]]; deletions = outputs[[2]]
  }
  
  subs_plcs = which(muts$after!="-" & muts$before!="-")
  if(length(subs_plcs)>0){
    outputs = update_matrix(mat = subs,plcs = subs_plcs,count = subs_count,muts = muts)
    subs_count = outputs[[1]]; subs = outputs[[2]]
  }
 }


