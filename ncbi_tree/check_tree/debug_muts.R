####3
# Functions:
get_muts = function(bases,node_depth,site,save_path){
  parents = node_depth$parent
  mut_plcs = (which(bases!=bases[parents]))
  if(length(mut_plcs)>0){
    muts = cbind(bases[parents][mut_plcs],bases[mut_plcs],parents[mut_plcs],names(bases)[mut_plcs])
    muts = data.frame(muts)
    colnames(muts) = c("before","after","parent","base")
    setwd(save_path)
    save(muts,file = paste0("muts.",site,".",dim(muts)[1],collapse = ""))
  }else{
    print("no_muts")
    print(site)
    muts = NULL
  }
  return(muts)
}

get_problematic_muts = function(muts,bases,node_depth,site){
  parents = node_depth$parent
  mut_plcs = (which(bases!=bases[parents]))
  # problematic_muts = which(is.element(muts$parent,muts$base))
  # problematic_muts = c(problematic_muts,which(is.element(muts$base,muts$parent)))
  problematic_muts = NULL
  children = node_depth$children[mut_plcs]
  children_bases = matrix(0,length(mut_plcs),5)
  colnames(children_bases) = c("A","C","G","T","-")
  for (i in 1:length(mut_plcs)){
    curr_children = strsplit(children[i],split = "\\*")[[1]]
    tmp = table(bases[curr_children])
    if(length(tmp)>0){
      for (j in 1:length(names(tmp))){
        children_bases[i,names(tmp)[j]] = tmp[j]
      }
    }
  }
  for (i in 1:length(mut_plcs)){
    if (children_bases[i,muts$before[i]]>0){
      problematic_muts = c(problematic_muts,i)
    }
  }
  problematic_muts = unique(problematic_muts)
  problematic_muts_mat = muts[problematic_muts,]
  return(problematic_muts)
}

update_matrix = function(mat,plcs,count,muts){
  mat[((count+1):(count+length(plcs))),] = as.matrix(muts[plcs,1:7])
  count = count+length(plcs)
  # colnames(mat) = colnames(muts)
  # mat = data.frame(mat)
  return(list(count,mat))
}
##################

begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd("/Users/keren/Desktop/covid_files/new2/ncbi trees/")
file_names = dir()
tmp = strsplit(file_names,"\\.")
nums = NULL
for (i in 1:length(tmp)){
  nums = c(nums,as.numeric(tmp[[i]][2]))
}
tree_list_nums = unique(nums)

setwd(paste0(begin_path,"vars/"))
load("node_depth")
load("site_to_alignment")
load("site_details")
setwd(paste0(begin_path,"vars/muts/"))
file_names = dir()
for (i in 1:length(file_names)){
  tryCatch(load(file_names[i]),error=function(cond){file.remove(file_names[i])})
}
file_names = dir()
tmp = strsplit(file_names,"\\.")
num_muts = NULL
for (i in 1:length(tmp)){
  # print(i)
  num_muts = c(num_muts,as.numeric(tmp[[i]][2]))
}

missing_muts = tree_list_nums[which(is.element(tree_list_nums,unique(num_muts))==FALSE)]
for (i in 1:length(missing_muts)){
  setwd("/Users/keren/Desktop/covid_files/new2/ncbi trees/")
  print(i)
  load(paste0("tree_list.",missing_muts[i]))
  muts = get_muts(bases = bases,node_depth = node_depth,site = missing_muts[i],save_path = 
                    paste0(begin_path,"vars/muts/"))
}

# # Not all tree_list files have corresponding muts files because some result in no muts, example:
# setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
# load("all_site_patterns")
# load("all_site_patterns_plcs")
# table(strsplit(site_patterns[which(site_patterns_plcs==9862)],split ="")[[1]])

#######################
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/muts/"))
file_names = dir()
tmp = strsplit(file_names,"\\.")
for (i in 1:length(tmp)){
  # Check file names make sense
  if(is.na(as.numeric(tmp[[i]][2]))){
    print(i)
  }
}
num_muts = NULL
for (i in 1:length(tmp)){
  # print(i)
  num_muts = c(num_muts,as.numeric(tmp[[i]][3]))
}
nums = NULL
insertions = deletions = subs = data.frame(matrix(0,700000,7))
colnames(insertions) = colnames(deletions) = colnames(subs) = 
  c("before","after","parent","base","alignment_site","ref_site","problem")
insertions_count = deletions_count = subs_count = 0

codon_sites = regions$ref_site[which(regions$gene!=0)]
for (i in 1:length(tmp)){
  start_time = Sys.time()
  print(i)
  nums = c(nums,as.numeric(tmp[[i]][2]))
  setwd(paste0(begin_path,"vars/muts/"))
  load(file_names[i])
  alignment_site = matrix(as.numeric(tmp[[i]][2]),dim(muts)[1],1)
  ref_site = matrix(site_to_alignment$site[as.numeric(tmp[[i]][2])],dim(muts)[1],1)
  setwd("/Users/keren/Desktop/covid_files/new2/ncbi trees/")
  if (file.exists(paste0("tree_list.",tmp[[i]][2]))){
    load(paste0("tree_list.",tmp[[i]][2]))
  }else{
    load(paste0("tree_list.",tmp[[i]][2],".filepart"))
  }
  
  if (is.element(ref_site[1],codon_sites)){
    problematic_muts = get_problematic_muts(muts,bases,node_depth,alignment_site)
    problem = matrix(0,dim(muts)[1],1)
    problem[problematic_muts,1] = 1
    muts = cbind(muts,alignment_site,ref_site,problem)
    insertions_plcs = which(muts$before=="-")
    if (length(insertions_plcs)>0){
      outputs = update_matrix(mat = insertions,plcs = insertions_plcs,count = insertions_count,muts = muts)
      insertions_count = outputs[[1]]; insertions = outputs[[2]]
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
  rm(bases)
  rm(muts)
  print(Sys.time()-start_time)
}

setwd(paste0(begin_path,"vars/debug_muts"))
deletions = deletions[1:(which(deletions$before==0)[1]-1),]
save(deletions,file = "deletions")
insertions = insertions[1:(which(insertions$before==0)[1]-1),]
save(insertions,file = "insertions")
subs = subs[1:(which(subs$before==0)[1]-1),]
save(subs,file = "subs")


