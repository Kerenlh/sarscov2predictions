begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/site_patterns/all2/",collapse = ""))
load("all_site_patterns")
load("all_site_patterns_plcs")

plc = which(site_patterns_plcs==19984)
site_patterns[plc]
deletion_site_pattern = strsplit(site_patterns[plc],split = "")[[1]]
del_rm_plcs = which(deletion_site_pattern!="-" & deletion_site_pattern!="A")
deletion_site_pattern = deletion_site_pattern[-del_rm_plcs]

correlated_plcs = NULL
for (i in 1:length(site_patterns)){
  curr_pattern = strsplit(site_patterns[i],split = "")[[1]][-del_rm_plcs]
  table_2_2 = table(curr_pattern,deletion_site_pattern)
  if (sum(dim(table_2_2)==c(2,2))==2){
    if (table_2_2[1,2]<=20 & table_2_2[2,1]<=20){
      print(i)
      correlated_plcs = c(correlated_plcs,i)
    }
  }
  sparse_vals = rownames(table_2_2)[which(table_2_2<=20,arr.ind = TRUE)[,1]]
  sparse_vals2 = rownames(table_2_2)[which(apply(table_2_2,1,sum)<=100)]
  sparse_vals = c(intersect(sparse_vals,sparse_vals2),"N","K","Y","R","S","W","M","B","D","H","V")
  rm_plcs = c(which(is.element(curr_pattern,sparse_vals)))
  table_2_2 = table(curr_pattern[-rm_plcs],deletion_site_pattern[-rm_plcs])
  if (dim(table_2_2)[1]>1){
    tmp = chisq.test(curr_pattern[-rm_plcs],deletion_site_pattern[-rm_plcs])
    if (tmp$p.value<0.001){
      print(i)
      correlated_plcs = c(correlated_plcs,i)
    }
  }
}
setwd(paste0(begin_path,"check_tree/"))
save(correlated_plcs,file = "correlated_plcs")

setwd(paste0(begin_path,"vars/"))
load("codons_table")
setwd(paste0(begin_path,"vars/tables_data/"))
load("deletions_coding")
load("insertions_coding")
load("events_table")
events_table = events_table[-which(events_table$base=="?"),]

j=15
curr_pattern = strsplit(site_patterns[correlated_plcs[j]],split = "")[[1]][-del_rm_plcs]
table_2_2 = table(curr_pattern,deletion_site_pattern)
sparse_vals = rownames(table_2_2)[which(table_2_2<=20,arr.ind = TRUE)[,1]]
sparse_vals2 = rownames(table_2_2)[which(apply(table_2_2,1,sum)<=100)]
sparse_vals = c(intersect(sparse_vals,sparse_vals2),"N","K","Y","R","S","W","M","B","D","H","V")
rm_plcs = c(which(is.element(curr_pattern,sparse_vals)))
table_2_2 = table(curr_pattern[-rm_plcs],deletion_site_pattern[-rm_plcs])
tmp = chisq.test(curr_pattern[-rm_plcs],deletion_site_pattern[-rm_plcs])
print(correlated_plcs[j])
print(table_2_2)
print(tmp)
correlated_site = site_patterns_plcs[correlated_plcs[j]]
correlated_site = 22582
codons_table[which(codons_table$site==correlated_site),]
deletions[which(deletions$start==correlated_site),c("parent_name","child_name")]
insertions[which(insertions$start==correlated_site),]
events_table[which(events_table$site==correlated_site),c("parent_name","child_name","base","A","C","G","T")]

curr_seqs = c("MT834220.1","MT276330.2","MT920008.1","MT632964.1","MT627611.1")
for (i in 1:length(curr_seqs)){
  print(strsplit(as.character(aligned_seqs[[curr_seqs[i]]]),split = "")[[1]][22582])
}


plcs = which(deletions$start==746)
deletions[plcs,c("parent_name","child_name")]

setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
load("site_details")
tree = read.tree("tree (1).nwk")
setwd(paste0(begin_path,"vars/big_trees/"))
load("all_big_trees")

node_names = names(all_tree_list)
bases_22582 = bases_746 = NULL
for (i in 1:length(node_names)){
  print(i)
  bases_22582 = c(bases_22582,strsplit(all_tree_list[[node_names[i]]]@probs,split="")[[1]][22582])
  bases_746 = c(bases_746,strsplit(all_tree_list[[node_names[i]]]@probs,split="")[[1]][746])
}
tmp = cbind(bases_22582,bases_746)
table(bases_22582,bases_746)

setwd("/Users/keren/Desktop/covid_files/covid_new/trees/")
load("tree_list.22582")

plcs = which(deletions$start==746)
deletions[plcs,c("parent_name","child_name")]
for (i in 1:length(plcs)){
  print(paste0("Child:",
               tree_list[[deletions[plcs[i],"child_name"]]]@base))
  print(tree_list[[deletions[plcs[i],"child_name"]]]@updated_probs)
  print(paste("Parent:",tree_list[[deletions[plcs[i],"parent_name"]]]@base))
  print(tree_list[[deletions[plcs[i],"parent_name"]]]@updated_probs)
  print("################")
}

# table(codons_table$site)
tmp = names(sort(table(unlist(events_table$site[which(events_table$codon.pos!=0)])),decreasing = TRUE)[1:100])

setwd("/Users/keren/Desktop/covid_files/covid_new/trees/")
# j=11161
j=29023
load(paste0("tree_list.",j))

plcs = which(events_table$site==j)
events_table[plcs,c("parent_name","child_name")]
for (i in 1:length(plcs)){
  print(paste0("Child:",
               tree_list[[events_table[plcs[i],"child_name"]]]@base))
  print(tree_list[[events_table[plcs[i],"child_name"]]]@updated_probs)
  print(paste("Parent:",tree_list[[events_table[plcs[i],"parent_name"]]]@base))
  print(tree_list[[events_table[plcs[i],"parent_name"]]]@updated_probs)
  print("################")
}

#############################
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

setClass(Class='Tree', 
         representation(name= 'character',node = 'numeric', parent_node = 'numeric',
                        parent_name = 'character',  
                        children_nodes = 'matrix', children_names = 'character', 
                        probs = 'matrix', updated_probs = 'matrix', base = 'character',leaf = 'logical'))

root_node = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
root_name = paste0(c("node_",root_node),collapse = "")
all_tree_list[[root_name]]@base = ""
all_tree_list = node_plc_name(all_tree_list, root_name)
setwd(paste0(begin_path,"check_tree/"))
save(all_tree_list,file = "all_tree_list_node_plcs")

j=15
correlated_site = site_patterns_plcs[correlated_plcs[j]]
plcs_del = which(deletions$start==correlated_site)
plcs_insertion = which(insertions$start==correlated_site)
plcs_subs = which(events$site==correlated_site)
node_plc_names_mat = matrix("-",length(plcs_del),300)
max_length = 0
for (i in 1:length(plcs_del)){
  tmp = strsplit(all_tree_list[[deletions[plcs_del[i],"child_name"]]]@base,split = "")[[1]]
  node_plc_names_mat[i,1:(1+length(tmp))] = c(length(tmp),tmp)
  max_length = max(max_length,length(tmp))
}

node_plc_names_mat = data.frame(node_plc_names_mat[,1:max_length])
node_plc_names_mat = node_plc_names_mat[order(as.numeric(node_plc_names_mat[,1])),]
node_plc_names_mat = cbind(deletions[plcs_del,c("parent_name","child_name")],node_plc_names_mat)
library(googlesheets4)
sheet_write(node_plc_names_mat)





