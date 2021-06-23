# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biostrings")

cluster = 3
if (cluster ==1){
  begin_path = "/a/home/cc/math/kerenlh/covid_new/"
}else if(cluster == 0) {
  begin_path = "/Users/keren/Dropbox/covid/new2/"
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
  begin_path = "C:/Users/ifog/Dropbox/covid/new/"
}

setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
load("all_site_patterns")
load("all_site_patterns_plcs")
unique_site_patterns = site_patterns[which(duplicated(site_patterns)==FALSE)]
site_patterns_numbers = site_patterns_plcs[which(duplicated(site_patterns)==FALSE)]
rm(site_patterns)

setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/",collapse = ""))
load("all_aligned_seqs")
seqs_names = names(aligned_seqs)
rm(aligned_seqs)
tree = read.tree("ft_SH.tree")
tree = short_names_tree(tree)

wait = runif(1, 0, 10)
Sys.sleep(wait)
tmp = as.numeric(Sys.time())
set.seed(tmp)

setwd(paste0(begin_path,"vars/trees/",collapse = ""))
# setwd("/Users/keren/Desktop/covid_files/covid_new/trees/")
file_names = dir()
if(isEmpty(file_names)){
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
}

if (cluster ==1){
  cluster_nums = site_patterns_numbers[which(site_patterns_numbers<20000)]
}else if(cluster == 0) {
  cluster_nums = site_patterns_numbers[which(site_patterns_numbers>=20000 |
                                               site_patterns_numbers<=15000)]
}else if (cluster ==2) {
  cluster_nums = site_patterns_numbers[which(site_patterns_numbers>10000)]
}else if(cluster == 3) {
  cluster_nums = site_patterns_numbers[which(site_patterns_numbers>15000 |
                                               site_patterns_numbers<=10000)]
  cluster_nums = NULL
}

if (file.exists("dropbox_desktop_files")){
  load("dropbox_desktop_files")
}else{
  dropbox_desktop_files = NULL
}
nums = unique(c(nums,in_progress_files,cluster_nums,dropbox_desktop_files))
missing_nums_index = which(is.element(site_patterns_numbers,nums)==FALSE)
missing_nums = site_patterns_numbers[missing_nums_index]
print("missing_files:")
print(missing_nums)
print("in_progress:")
#print(in_progress_files)
print(paste("length(missing_nums):",length(missing_nums)))
# curr_nums = missing_nums[sample(1:length(missing_nums),50)]
curr_nums = missing_nums
# curr_nums = c(503 , 1322  ,5475,9911, 10118, 10145, 10282, 22184, 
#               22484 ,25811, 26129, 28228)
curr_nums = 10145
print(curr_nums)
in_progress_files = c(in_progress_files,curr_nums)
save(in_progress_files,file = "in_progress_files")

#############
setClass(Class='Tree', 
         representation(name= 'character',node = 'numeric', parent_node = 'numeric',
                        parent_name = 'character',  
                        children_nodes = 'matrix', children_names = 'character', 
                        probs = 'matrix', updated_probs = 'matrix', base = 'character',leaf = 'logical'))


find_details_not_leaf = function(curr_name,curr_site_pattern,
                                 tree,tree_list,root_node){
  # print(curr_name)
  tmp = strsplit(curr_name,split = "")
  if (paste0(tmp[[1]][1:5],collapse = "")=="node_"){
    curr_node = as.numeric(paste0(tmp[[1]][6:length(tmp[[1]])],collapse = ""))
  }else{
    curr_node = which(tree$tip.label==curr_name)
    leaf_bug_check = 1
  }
  curr_children_nodes = as.matrix(tree$edge[which(tree$edge[,1]==curr_node),2])
  if (length(curr_children_nodes)==0){ # This is a leaf
    if (leaf_bug_check==0){
      print("This is a leaf identification bug")
    }
    tree_list = find_details_leaf(curr_name,curr_site_pattern,tree,tree_list)
    return(tree_list)
  }else{
    leaf = FALSE
    curr_parent_node = tree$edge[which(tree$edge[,2]==curr_node),1]
    curr_parent_name = tree$tip.label[curr_parent_node]
    if (curr_node == root_node){ # This is the root
      curr_parent_name = "node_0"
    }
    if(is.na(curr_parent_name)){
      curr_parent_name = paste0("node_",curr_parent_node,collapse = "")
    }
    curr_children_names = NULL
    for (i in 1:length(curr_children_nodes)){
      if (curr_children_nodes[i]<=length(tree$tip.label)){
        curr_children_names = c(curr_children_names,tree$tip.label[curr_children_nodes[i]])
      }else{
        curr_children_names = c(curr_children_names,paste0("node_",curr_children_nodes[i]))
      }
    }
    all_probs = matrix(0,1,5)
    for (i in 1:length(curr_children_names)){
      if (is.null(tree_list[[curr_children_names[i]]])){
        tree_list = find_details_not_leaf(curr_children_names[i],curr_site_pattern,
                                          tree,tree_list,root_node)
      } 
      all_probs = tree_list[[curr_children_names[i]]]@probs + all_probs
    }
    probs = all_probs/length(curr_children_names)
    updated_probs = matrix(0,1,5)
    base = colnames(probs)[which.max(probs)]
    if(length(which(probs==max(probs)))>1){
      print(paste(curr_name, probs))
    }
    tt = new("Tree", name = curr_name, node = curr_node, parent_node = curr_parent_node,
             parent_name = curr_parent_name, children_nodes = curr_children_nodes,
             children_names = curr_children_names, base = base, probs = probs,
             updated_probs = updated_probs, leaf = leaf)
    tree_list[[curr_name]] = tt
    return(tree_list)
  }
}

find_details_leaf = function(curr_name,curr_site_pattern,tree,
                             tree_list){
  #print(curr_name)
  curr_node = which(tree$tip.label==curr_name)
  curr_parent_node = tree$edge[which(tree$edge[,2]==curr_node),1]
  curr_parent_name = tree$tip.label[curr_parent_node]
  if(is.na(curr_parent_name)){
    curr_parent_name = paste0("node_",curr_parent_node,collapse = "")
  }
  curr_children_nodes = tree$edge[which(tree$edge[,1]==curr_node),2]
  if (length(curr_children_nodes)>0 | curr_node>length(tree$tip.label)){
    print("bug, this is not a leaf!")
  }else{
    base = curr_site_pattern[which(seqs_names==curr_name)]
    probs = get_probs_from_base(base)
    updated_probs = matrix(0,1,5)
    leaf = TRUE
  }
  tt = new("Tree", name = curr_name, node = curr_node, parent_node = curr_parent_node,
           parent_name = curr_parent_name, children_nodes = matrix(0,1,1),
           children_names = "_", base = base, probs = probs,
           updated_probs = updated_probs, leaf = leaf)
  tree_list[[curr_name]] = tt
  return(tree_list)
}

get_probs_from_base = function(base){
  probs = matrix(0,1,5)
  colnames(probs) = c("A","C","G","T","-")
  if (is.element(base,colnames(probs))){
    probs[,base] = 1
  }else if (base=="R"){
    probs[,c("A","G")] = 0.5
  }else if (base=="Y"){
    probs[,c("C","T")] = 0.5
  }else if (base=="S"){
    probs[,c("C","G")] = 0.5
  }else if (base=="W"){
    probs[,c("A","T")] = 0.5
  }else if (base=="K"){
    probs[,c("T","G")] = 0.5
  }else if (base=="M"){
    probs[,c("A","C")] = 0.5
  }else if (base=="B"){
    probs[,c("C","T","G")] = 1/3
  }else if (base=="D"){
    probs[,c("A","T","G")] = 1/3
  }else if (base=="H"){
    probs[,c("C","T","A")] = 1/3
  }else if (base=="V"){
    probs[,c("C","A","G")] = 1/3
  }else if (base=="N"){
    probs[,c("C","A","G","T")] = 1/4
  }
  return(probs)
}

update_probs = function(tree_list,curr_name,epsilon){
  if (tree_list[[curr_name]]@leaf==FALSE){
    children = tree_list[[curr_name]]@children_names
    for (i in 1:length(children)){
      tree_list[[children[i]]]@updated_probs = 
        (tree_list[[children[i]]]@probs + epsilon*tree_list[[curr_name]]@probs)/
        (1+epsilon)
      tmp = tree_list[[children[i]]]@probs
      if(length(which(tmp==max(tmp)))>1){
        # print(paste(children[i], tmp, tree_list[[children[i]]]@updated_probs))
        tree_list[[children[i]]]@base = 
          colnames(tmp)[which.max(tree_list[[children[i]]]@updated_probs)]
        #print(tree_list[[children[i]]]@base)
      }
      tree_list = update_probs(tree_list,children[i],epsilon)        
    }
  }
  return(tree_list) 
}

root_node = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
root_name = paste0(c("node_",root_node),collapse = "")
setwd(paste0(begin_path,"vars/old_fitch_trees/",collapse = ""))
epsilon = 0.01
curr_nums = c(3,11,1400,1528,1625,1742,3067,24227,25739,26411,27381)
for (k in 1:length(curr_nums)){
  j = which(site_patterns_numbers==curr_nums[k])
  print(j)
  tree_list = list()
  curr_site_pattern = strsplit(unique_site_patterns[j],split = "")[[1]]
  start_time = Sys.time()
  tree_list = find_details_not_leaf(root_name,curr_site_pattern,
                                    tree,tree_list,root_node)
  tree_list = update_probs(tree_list, root_name,epsilon)
  save(tree_list,file = paste0("tree_list.",curr_nums[k],collapse=""))
  end_time = Sys.time()
  print(end_time-start_time)
  print(paste0("tree_list.",curr_nums[k],collapse=""))
  print(k)
  # load("in_progress_files")
  # in_progress_files = in_progress_files[-which(in_progress_files==curr_nums[k])]
  # save(in_progress_files,file = "in_progress_files")
}

######
# Debug:
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/old_fitch_trees/",collapse = ""))
load("tree_list.27381")
setwd(paste0(begin_path,"vars/"))
load("node_depth")
load("all_al")
tree_list_fitch_non_binary = tree_list
setwd("/Users/keren/Desktop/covid_files/new2/Lanfear trees/")
load("tree_list.27381")
old_alg = matrix(0,length(bases),1)
bases = cbind(bases,old_alg)
colnames(bases) = c("tree_reconstruction","fitch")
bases = data.frame(bases)
bases_rows = row.names(bases)
for (i in 1:dim(bases)[1]){
  print(i)
  bases[i,"fitch"] = tree_list_fitch_non_binary[[bases_rows[i]]]@base
}

