node_plc_name = function(tree_list,curr_name){
  print(curr_name)
  if (tree_list[[curr_name]]@leaf==FALSE){
    children = tree_list[[curr_name]]@children_names
    print(children)
    for (i in 1:length(children)){
      tree_list[[children[i]]]@base =
        paste0(tree_list[[curr_name]]@base,i,collapse = "")
      print(tree_list[[curr_name]]@base)
      tree_list = node_plc_name(tree_list,children[i])
    }
  }
  return(tree_list)
}


begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("tree_list.1_trimmed")
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("trimmed_tree")
root_node = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
root_name = paste0(c("node_",root_node),collapse = "")
tree_list[[root_name]]@base = ""
tree_list = node_plc_name(tree_list, root_name)
save(tree_list,file = "tree_list_node_plcs")
tree = short_names_tree(tree)

node_depth = data.frame(matrix(0,length(tree_list),2))
node_depth[,1] = c(1:length(tree_list))
colnames(node_depth) = c("name","depth")
for (i in 1:length(tree_list)){
  print(i)
  node_depth$depth[i] = length(strsplit(tree_list[[i]]@base,split = "")[[1]])
  node_depth$name[i] = tree_list[[i]]@name
}
rownames(node_depth) = node_depth$name
node_depth = node_depth[order(node_depth$depth,decreasing = TRUE),]

load("all_aligned_seqs")
seqs_names = names(aligned_seqs)
rm(aligned_seqs)
seqs_names_index = matrix(0,dim(node_depth)[1],1)
node_depth = cbind(node_depth,seqs_names_index)
for (i in 1:dim(node_depth)[1]){
  tmp = which(seqs_names==node_depth$name[i])
  if (length(tmp)>0){
    node_depth$seqs_names_index[i] = tmp
  }
}
parent = children = matrix("m",dim(node_depth)[1],1)
node_depth = cbind(node_depth,parent,children)
for (i in 1:dim(node_depth)[1]){
  print(i)
  node_depth$parent[i] = tree_list[[node_depth$name[i]]]@parent_name
  curr_children = tree_list[[node_depth$name[i]]]@children_names
  if(curr_children=="_"){
    node_depth$children[i] = NA
  }else{
    node_depth$children[i] = paste0(curr_children,collapse = "*")
  }
}
num_children = matrix(0,dim(node_depth)[1],1)
for (i in 1:dim(node_depth)[1]){
  if (is.na(node_depth$children[i])==FALSE){
    num_children[i] = length(strsplit(node_depth$children[i],split = "\\*")[[1]])
    print(i)
    print(num_children[i])
  }
}
node_depth = cbind(node_depth,num_children)
setwd(paste0(begin_path,"vars/"))
save(node_depth,file = "node_depth")
