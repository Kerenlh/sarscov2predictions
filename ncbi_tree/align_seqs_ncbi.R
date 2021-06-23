begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
seqs_lanfear = readDNAStringSet(file = "global.fa")
seqs_lanfear = short_names_seqs(seqs_lanfear)
# Ref_seq name at Lanfear's tree: NC_045512.2_|Severe_acute_respiratory_syndrome_coronavirus_2_isolate_Wuhan-Hu-1
setwd("/Users/keren/Dropbox/covid/new2/ncbi_tree/vars/")
tree = read.tree("tree_18_2.nwk") 
tree = short_names_tree(tree)
aligned_seqs = return_only_tree_seqs(seqs_lanfear,tree)
save(aligned_seqs,file = "all_aligned_seqs")

#################
# Trim tree!
setwd("/Users/keren/Dropbox/covid/new2/")
source("COVID_functions.R")
setwd("/Users/keren/Dropbox/covid/new2/ncbi_tree/vars/")
load("all_aligned_seqs")
tree = read.tree("tree_18_2.nwk") 
tree = short_names_tree(tree)
tmp = sort(table(tree$edge[,1])) # There are at least 2 children for each node
missing_leaves = tree$tip.label[which(is.element(tree$tip.label,names(aligned_seqs))==FALSE)]
load("tree_list.1")

delete_node = function(curr_name,tree,tree_list){
  parent_name = tree_list[[curr_name]]@parent_name
  missing_child_plc = which(tree_list[[parent_name]]@children_names==curr_name)
  missing_node_plc = which(tree_list[[parent_name]]@children_nodes==tree_list[[curr_name]]@node)
  tree_list[[parent_name]]@children_names = tree_list[[parent_name]]@children_names[-missing_child_plc]
  tree_list[[parent_name]]@children_nodes = as.matrix(tree_list[[parent_name]]@children_nodes[-missing_node_plc,])
  rel_edge_row = which(tree$edge[,2]==tree_list[[curr_name]]@node)
  tree$edge = tree$edge[-rel_edge_row,]
  if (length(tree_list[[curr_name]]@children_names)==0){
    tree_list_names = names(tree_list)
    tree_list = tree_list[-which(tree_list_names==curr_name)]
    return(list(tree,tree_list,parent_name))
  }
  else if (tree_list[[curr_name]]@children_names=="_"){
    tree_list_names = names(tree_list)
    tree_list = tree_list[-which(tree_list_names==curr_name)]
    return(list(tree,tree_list,parent_name))
  }else{
    print(paste("This seq has children: ", curr_name))
    bug_count = bug_count+1
    break
  }
}

remove_inner_node = function(curr_name,tree,tree_list){
  curr_parent = tree_list[[curr_name]]@parent_name
  tree_list[[curr_parent]]@children_names[which(tree_list[[curr_parent]]@children_names==curr_name)] = 
    curr_children
  tree_list[[curr_parent]]@children_nodes[which(tree_list[[curr_parent]]@children_nodes==
                                                  tree_list[[curr_name]]@node),] =
    tree_list[[curr_children]]@node
  rel_edge_row = which(tree$edge[,2]==tree_list[[curr_name]]@node)
  tree$edge = tree$edge[-rel_edge_row,]
  return(list(tree,tree_list,curr_parent))
}

bug_count = 0
parents = NULL
for (i in 1: length(missing_leaves)){
  print(i)
  curr_name = missing_leaves[i]
  tmp = delete_node(curr_name = missing_leaves[i],tree = tree,tree_list = tree_list)
  tree = tmp[1][[1]]
  tree_list = tmp[2][[1]]
  parents = c(parents,tmp[3][[1]])
}
# save(tree,file = "trimmed_tree")
# save(tree_list,file = "tree_list.1_trimmed")
# 
# load("trimmed_tree")
# load("tree_list.1_trimmed")
parents = unique(parents)
new_parents = NULL
for (i in 1:length(parents)){
  print(i)
  curr_name = parents[i]
  curr_children = tree_list[[curr_name]]@children_names
  while(length(curr_children)<2){
    print(curr_name)
    if (length(curr_children)==0){
      tmp = delete_node(curr_name = curr_name,tree = tree,tree_list = tree_list)
      tree = tmp[1][[1]]
      tree_list = tmp[2][[1]]
      # new_parents = c(new_parents,tmp[3][[1]])
    }else if (length(curr_children)==1){
      tmp = remove_inner_node(curr_name = curr_name,tree = tree,tree_list = tree_list)
      tree = tmp[1][[1]]
      tree_list = tmp[2][[1]]
      # new_parents = c(new_parents,tmp[3][[1]])
    }
    curr_name = tmp[3][[1]] #new_parent
    curr_children = tree_list[[curr_name]]@children_names
  }
}

save(tree,file = "trimmed_tree")
save(tree_list,file = "tree_list.1_trimmed")
length(which(is.element(tree$edge[,2],tree$edge[,1])==FALSE))

# Trim edge.length:
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("trimmed_tree")
trimmed_tree = tree
orig_tree = read.tree("tree_18_2.nwk")
setwd(begin_path)
source("COVID_functions.R")
orig_tree = short_names_tree(orig_tree)
orig_edges = apply(orig_tree$edge,1,paste0,collapse = ".")
trimmed_edges = apply(trimmed_tree$edge,1,paste0,collapse = ".")

for (i in 1:dim(trimmed_tree$edge)[1]){
  if (is.element(trimmed_edges[i],orig_edges)==FALSE){
    print(i)
    tmp = strsplit(trimmed_edges[i],split = "\\.")[[1]]
    trimmed_parent = tmp[1]; trimmed_child = tmp[2]
    curr_child = trimmed_child
    total_edge_length = 0
    while(curr_child!=trimmed_parent){
      curr_edge_num = which(orig_tree$edge[,2] == curr_child)
      print(curr_edge_num)
      total_edge_length = total_edge_length+orig_tree$edge.length[curr_edge_num]
      curr_child = orig_tree$edge[curr_edge_num,1]
    }
    trimmed_tree$edge.length[i] = total_edge_length
  }
}
trimmed_tree$edge.length = trimmed_tree$edge.length[1:dim(trimmed_tree$edge)[1]]
tree = trimmed_tree
setwd(paste0(begin_path,"vars/"))
save(tree,file = "trimmed_tree")

#####################
# Find site patterns:
#####################
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
setwd(paste0(begin_path,"vars/site_patterns/",collapse = ""))
num_seqs = length(aligned_seqs)
seq_length = length(aligned_seqs[[1]])
for (j in 0:38){
  print(j)
  split_seqs = matrix("m",1000,length(aligned_seqs[[1]]))
  count = 0
  for (i in (j*1000+1):min(length(aligned_seqs),(j+1)*1000)){
    print(i)
    count = count+1
    split_seqs[count,] = strsplit(as.character(aligned_seqs[[i]]),split = "")[[1]]
  }
  if (j==38){
    split_seqs = split_seqs[1:(num_seqs-j*1000),]
  }
  save(split_seqs,file = paste0("split_seqs_",(j*1000+1),"_",(j+1)*1000,collapse=""))
}

##################
# Unite split_seqs: 
##################
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
setwd(paste0(begin_path,"vars/site_patterns/",collapse = ""))
file_names = dir()
################3
# Order file_names:
tmp = strsplit(file_names,split = "_")
file_nums = NULL
for (i in 1:length(file_names)){
  file_nums = c(file_nums,as.numeric(tmp[[i]][3]))
}
file_names = file_names[order(file_nums,decreasing = FALSE)]
file_names = file_names[-c(40,41)]  # to remove all,all2
##############
num_seqs = length(aligned_seqs)
seq_length = length(aligned_seqs[[1]])
rm(aligned_seqs)
count = 0
num_sites = 5000
print(file_names)
for (j in 1:6){
  plcs = c((1+(j-1)*num_sites):min(seq_length,(j*num_sites)))
  all_split_seqs = matrix("M",num_seqs,length(plcs))
  count = 0
  for (i in 1:length(file_names)){
    print(i)
    setwd(paste0(begin_path,"/vars/site_patterns/"))
    load(file_names[i])
    print(count)
    print(dim(split_seqs))
    all_split_seqs[(count+1):(dim(split_seqs)[1]+count),] = split_seqs[,plcs]
    count = count + dim(split_seqs)[1]
    rm(split_seqs)
  }
  split_seqs = all_split_seqs
  setwd(paste0(begin_path,"/vars/site_patterns/all"))
  save(split_seqs,file = paste0("all_split_seqs_",plcs[length(plcs)],collapse =""))
}
###########
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
seq_length = length(aligned_seqs[[1]])
rm(aligned_seqs)
setwd(paste0(begin_path,"vars/site_patterns/all/",collapse = ""))
file_names = dir()
################3
# Order file_names:
tmp = strsplit(file_names,split = "_")
file_nums = NULL
for (i in 1:length(file_names)){
  file_nums = c(file_nums,as.numeric(tmp[[i]][4]))
}
file_names = file_names[order(file_nums,decreasing = FALSE)]
##############
count = 0
for (j in 1:6){
  setwd(paste0(begin_path,"vars/site_patterns/all",collapse = ""))
  load(file_names[j])
  site_patterns = site_patterns_plcs = no_subs_plcs = NULL
  for (i in 1:dim(split_seqs)[2]){
    print(i)
    curr_site_pattern = split_seqs[,i]
    if ( length(table(curr_site_pattern))==1) { 
      no_subs_plcs = c(no_subs_plcs,(i+count))
    }else{
      site_patterns = c(site_patterns,paste0(curr_site_pattern, collapse = ""))
      site_patterns_plcs = c(site_patterns_plcs,(i+count)) 
    }
  }
  count = count+dim(split_seqs)[2]
  file_name = paste0("site_patterns_",j,collapse = "")
  setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
  save(site_patterns,file = file_name)
  save(site_patterns_plcs,file = paste0("site_patterns_plcs_",j,collapse = ""))
  save(no_subs_plcs,file = paste0("no_subs_plcs_",j,collapse = ""))
}

######33
# Unite site_patterns:
setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
file_names = dir()
all_site_patterns = NULL
for (i in 7:12){
  print(i)
  load(file_names[i])
  all_site_patterns = c(all_site_patterns,site_patterns)
}
site_patterns = all_site_patterns
save(site_patterns,file = "all_site_patterns")

all_site_patterns_plcs = NULL
for (i in 13:18){
  print(i)
  load(file_names[i])
  all_site_patterns_plcs = c(all_site_patterns_plcs,site_patterns_plcs)
}
site_patterns_plcs = all_site_patterns_plcs
save(site_patterns_plcs,file = "all_site_patterns_plcs")

all_no_subs_plcs = NULL
for (i in 1:6){
  print(i)
  load(file_names[i])
  all_no_subs_plcs = c(all_no_subs_plcs,no_subs_plcs)
}
no_subs_plcs = all_no_subs_plcs
save(no_subs_plcs,file = "all_no_subs_plcs")



