# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biostrings")

cluster = 0
if (cluster ==1){
  begin_path = "/a/home/cc/math/kerenlh/covid_new/"
}else if(cluster == 0) {
  begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
}else if (cluster ==2) {
  begin_path = "/a/home/cc/math/saharon/keren/"
  setwd(paste0(begin_path,"vars/"))
  files<-list.files(full.names = T)
  vect_size <- sapply(files, file.size)
  size_files <- sum(vect_size)
  if ((size_files/1000000)>250){
    stop("no more space on Sahron's disk")
  }
}else if (cluster ==3) {
  begin_path = "C:/Users/Keren/Dropbox/covid/new2/ncbi_tree/"
}

setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
load("all_site_patterns")
load("all_site_patterns_plcs")
unique_site_patterns = site_patterns[which(duplicated(site_patterns)==FALSE)]
site_patterns_numbers = site_patterns_plcs[which(duplicated(site_patterns)==FALSE)]
rm(site_patterns)

setwd(paste0(begin_path,"vars/"))
load("node_depth")
load("tree_list.1_trimmed")

wait = runif(1, 0, 10)
Sys.sleep(wait)
tmp = as.numeric(Sys.time())
set.seed(tmp)

setwd(paste0(begin_path,"vars/trees/",collapse = ""))
file_names = dir()
if(identical(file_names,character(0))){
  file_names = "tmp"
}
nums = NULL
for (i in 1:length(file_names)){
  if (strsplit(file_names[i],split = "[.]")[[1]][1]=="tree_list"){
    nums = c(nums,strsplit(file_names[i],split = "[.]")[[1]][2])
  }
}
nums = as.numeric(nums)

if (file.exists("in_progress_files")){
  load("in_progress_files")
}else{
  in_progress_files = NULL
  save(in_progress_files,file = "in_progress_files")
}

load("cluster_nums")
# if (cluster ==1){
#   cluster_nums = site_patterns_numbers[which(site_patterns_numbers<20000)]
# }else if(cluster == 0) {
#   cluster_nums = site_patterns_numbers[which(site_patterns_numbers>=20000 |
#                                                site_patterns_numbers<=15000)]
# }else if (cluster ==2) {
#   cluster_nums = site_patterns_numbers[which(site_patterns_numbers>10000)]
# }else if(cluster == 3) {
#   cluster_nums = site_patterns_numbers[which(site_patterns_numbers>15000 |
#                                                site_patterns_numbers<=10000)]
#   cluster_nums = NULL
# }

nums = unique(c(nums,in_progress_files,cluster_nums))
missing_nums_index = which(is.element(site_patterns_numbers,nums)==FALSE)
missing_nums = site_patterns_numbers[missing_nums_index]
# print("missing_files:")
# print(missing_nums)
print("in_progress:")
#print(in_progress_files)
print(paste("length(missing_nums):",length(missing_nums)))
curr_nums = missing_nums[sample(1:length(missing_nums),50)]
print(curr_nums)
in_progress_files = c(in_progress_files,curr_nums)
save(in_progress_files,file = "in_progress_files")

get_probs_from_bases = function(bases){
  probs = matrix(0,length(bases),5)
  colnames(probs) = c("A","C","G","T","-")
  for (i in 1:length(colnames(probs))){
    plcs = which(bases==colnames(probs)[i])
    probs[plcs,colnames(probs)[i]] = 1
  }
  plcs = which(bases=="R")
  probs[plcs,c("A","G")] = 0.5
  plcs = which(bases=="Y")
  probs[plcs,c("C","T")] = 0.5
  plcs = which (bases=="S")
  probs[plcs,c("C","G")] = 0.5
  plcs = which(bases=="W")
  probs[plcs,c("A","T")] = 0.5
  plcs = which (bases=="K")
  probs[plcs,c("T","G")] = 0.5
  plcs = which(bases=="M")
  probs[plcs,c("A","C")] = 0.5
  plcs = which(bases=="B")
  probs[plcs,c("C","T","G")] = 1/3
  plcs = which(bases=="D")
  probs[plcs,c("A","T","G")] = 1/3
  plcs = which(bases=="H")
  probs[plcs,c("C","T","A")] = 1/3
  plcs = which(bases=="V")
  probs[plcs,c("C","A","G")] = 1/3
  plcs = which(bases=="N")
  probs[plcs,c("C","A","G","T")] = 1/4
  return(probs)
}

leaves_plcs = which(is.na(node_depth$children))
inner_nodes_plcs = which(is.na(node_depth$children)==FALSE)

curr_nums = missing_nums
for (k in 1:length(curr_nums)){
  start_time = Sys.time()
  print(k)
  print(curr_nums[k])
  probs = data.frame(matrix("m",dim(node_depth)[1],5))
  row.names(probs) = row.names(node_depth)
  colnames(probs) = c("A","C","G","T","-")
  j = which(site_patterns_numbers==curr_nums[k])
  print(j)
  curr_site_pattern = strsplit(unique_site_patterns[j],split = "")[[1]]
  curr_bases = curr_site_pattern[node_depth$seqs_names_index[leaves_plcs]]
  probs[leaves_plcs,] = get_probs_from_bases(curr_bases)
  # check that node_depth is sorted by depth, decreasing = T
  for (i in 1:length(inner_nodes_plcs)){
    # print(i)
    children = tree_list[[node_depth$name[inner_nodes_plcs[i]]]]@children_names
    probs[inner_nodes_plcs[i],] = apply(matrix(as.numeric(unlist(probs[children,])),length(children),5),2,mean)
  }
  
  print(Sys.time()-start_time)
  probs = matrix(as.numeric(unlist(probs)),dim(node_depth)[1],5)
  colnames(probs) = c("A","C","G","T","-")
  row.names(probs) = row.names(node_depth)
  
  for (i in (length(inner_nodes_plcs)-1):1){ # going from the root (shortest depth to the leaves)
    # The root itself is not updated because it doesn't have a parent
    parent = node_depth$parent[inner_nodes_plcs[i]]
    epsilon = 1/(node_depth$num_children[inner_nodes_plcs[i]])
    probs[inner_nodes_plcs[i],] = (probs[inner_nodes_plcs[i],]+epsilon*probs[parent,])/(1+epsilon)  
    # print(i)
  }
  # The update order in the leaves doesn't matter
  parents = node_depth$parent[leaves_plcs]
  probs[leaves_plcs,] = (probs[leaves_plcs,]+probs[parents,])/2
  
  bases = colnames(probs)[apply(probs,1,which.max)]
  names(bases) = row.names(probs)
  tmp = apply(probs,1,max)
  # Check if there are plcs with equal probs
  count = 0
  for (i in 1:length(tmp)){
    if(length(which(probs[i,]==tmp[i]))>1){
      count = count+1
      print(i)
    }
  }
  
  parents = node_depth$parent
  mut_plcs = (which(bases!=bases[parents]))
  if (length(mut_plcs)>0){
    for (i in 1:length(mut_plcs)){
      curr_parent = parents[mut_plcs[i]]
      curr_children = strsplit(node_depth$children[mut_plcs[i]],split = "\\*")[[1]]
      if (is.na(curr_children[1])==FALSE){ # If there are no children don't take majority, there was already an epsilon correction from the father
        tmp = table(c(bases[curr_parent],bases[curr_children]))
        if(length(which(tmp==max(tmp)))==1){
          # if there is more than 1 maximum do nothing
          bases[mut_plcs[i]] = names(which.max(tmp))
        } 
      }
    }
    
    mut_plcs = (which(bases!=bases[parents]))
  }
  print(Sys.time()-start_time)
  # check # of mutations:
  if(length(mut_plcs)>0){
    muts = cbind(bases[parents][mut_plcs],bases[mut_plcs],parents[mut_plcs],names(bases)[mut_plcs])
    muts = data.frame(muts)
    colnames(muts) = c("before","after","parent","base")
    children = node_depth$children[mut_plcs]
    row.names(muts) = muts$base
    ######
    # Check if one of the sons is like the parent (the son's grandparent)
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
    problematic_muts = NULL
    for (i in 1:length(mut_plcs)){
      if (children_bases[i,muts$before[i]]>0){
        problematic_muts = c(problematic_muts,i)
      }
    }
    problematic_muts = c(problematic_muts,which(is.element(muts$parent,muts$base)))
    problematic_muts = c(problematic_muts,which(is.element(muts$base,muts$parent)))
    problematic_muts = unique(problematic_muts)
  }
  setwd(paste0(begin_path,"vars/trees/",collapse = ""))
  save(bases,file = paste0("tree_list.",curr_nums[k],collapse=""))
  if (length(mut_plcs)>0){
    save(muts,file = paste0("muts.",curr_nums[k],".",dim(muts)[1],collapse = ""))
    if(length(problematic_muts)>0){
      problematic_muts_mat = muts[problematic_muts,]
      save(problematic_muts_mat,
           file = paste0("problematic_muts.",curr_nums[k],".",length(problematic_muts),collapse = ""))
    }
  }
  end_time = Sys.time()
  print(end_time-start_time)
  print(paste0("tree_list.",curr_nums[k],collapse=""))
  print(k)
  # load("in_progress_files")
  # in_progress_files = in_progress_files[-which(in_progress_files==curr_nums[k])]
  # save(in_progress_files,file = "in_progress_files")
}



# Debug:
# file_names = dir()
# num_muts = NULL
# for (i in 1:length(file_names)){
#   if (strsplit(file_names[i],split = "\\.")[[1]][1]=="muts"){
#     load(file_names[i])
#     num_muts = c(num_muts,dim(muts)[1])
#     rm(muts)
#   }
# }
# setwd("/Users/keren/Desktop/covid_files/covid_new/trees/")
# load("tree_list.100")
# count = 0
# tree_list_muts = 0
# for (i in 1:length(tree_list)){
#   curr_name = tree_list[[i]]@name
#   parent = tree_list[[i]]@parent_name
#   if (tree_list[[i]]@base!=bases[curr_name]){
#     print(paste(i,curr_name,"old: ",tree_list[[i]]@base,"new: ",bases[curr_name]))
#     count = count+1
#   }
#   if (tree_list[[i]]@base!=tree_list[[parent]]@base){
#     tree_list_muts = tree_list_muts+1
#   }
# }
