####
# remove outlier branch!

# begin_path = "C:/Users/ifog/Dropbox/covid/new/"
begin_path = "/Users/keren/Dropbox/covid/new2/"

get_node_names = function(tree_list){
  nodes_names = NULL
  for (i in 1:length(tree_list)){
    print(i)
    nodes_names = rbind(nodes_names,cbind(tree_list[[i]]@name, tree_list[[i]]@node))
  }
  colnames(nodes_names) = c("name","number")
  nodes_names = data.frame(nodes_names)
  nodes_names$number = as.numeric(nodes_names$number)
  return(nodes_names)
}
# setwd(paste0(begin_path,"vars/trees/"))
# load("tree_list.10012")
# nodes_names = get_node_names(tree_list)
# setwd(paste0(begin_path,"vars/"))
# save(nodes_names,file = "nodes_names")

setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("nodes_names")
load("site_details")
tree = read.tree("ft_SH.tree")
tree = short_names_tree(tree)
# load("all_aligned_seqs")
setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
load("all_seqs")


############
# Functions:
update_states_and_outputs = function(curr_states,states,curr_outputs,
                                     outputs,uncertainty_plcs){
  if (length(uncertainty_plcs)>0){
    curr_states = curr_states[-uncertainty_plcs]
    curr_outputs = curr_outputs[-uncertainty_plcs,]
  }
  row.names(curr_outputs) = curr_states
  tmp = is.element(curr_states,states)
  new_states_plcs = which(tmp==FALSE)
  old_states_plcs = which(tmp==TRUE)
  if(length(old_states_plcs)>0){
    outputs[curr_states[old_states_plcs],] = 
      outputs[curr_states[old_states_plcs],] + curr_outputs[old_states_plcs,]
  }
  outputs = rbind(outputs,curr_outputs[new_states_plcs,])
  states = c(states,curr_states[new_states_plcs])
  return(list(states = states,outputs = outputs))
}

find_indels_plcs = function(indel_plcs,curr_edge,rel_seq){
  if (length(indel_plcs)>0){
    diff_indel_plcs = diff(indel_plcs)
    big_diff_plcs = which(diff_indel_plcs>1)
    start_plcs = c(indel_plcs[1],(indel_plcs[big_diff_plcs+1]))
    end_plcs = c(indel_plcs[big_diff_plcs],indel_plcs[length(indel_plcs)])
    muts = NULL
    for (i in 1:length(start_plcs)){
      muts = rbind(muts,paste0(rel_seq[start_plcs[i]:end_plcs[i]],collapse = ""))
    }
    output = cbind(start_plcs,end_plcs,muts,matrix(curr_edge,length(start_plcs),2,byrow = TRUE))
    colnames(output)[3:5] = c("mutation","parent_node","child_node")
    return(output)
  }else{
    return(NULL)
  }
}
##########

insertions = deletions = NULL
for (i in 1:dim(tree$edge)[1]){
  # if (i==44387){
  #   # outlier edge, seq = "MT825091.1"
  #   next
  # }
  print(i)
  curr_edge = tree$edge[i,]
  parent_name = nodes_names$name[which(nodes_names$number==curr_edge[1])]
  child_name = nodes_names$name[which(nodes_names$number==curr_edge[2])]
  edge_length = tree$edge.length[i]
  parent_seq = strsplit(all_seqs[parent_name][[1]],split = "")[[1]]
  child_seq = strsplit(all_seqs[child_name][[1]],split = "")[[1]]
  subs_plcs = which((parent_seq==child_seq)==FALSE)
  insertion_plcs = subs_plcs[which(parent_seq[subs_plcs]=="-" & 
                                     (child_seq[subs_plcs]=="A" | child_seq[subs_plcs]=="C" |
                                        child_seq[subs_plcs]=="G" | child_seq[subs_plcs]=="T"))]
  deletion_plcs = subs_plcs[which(child_seq[subs_plcs]=="-" & 
                                    (parent_seq[subs_plcs]=="A" | parent_seq[subs_plcs]=="C" |
                                       parent_seq[subs_plcs]=="G" | parent_seq[subs_plcs]=="T" ))]
  new_insertions = find_indels_plcs(insertion_plcs,curr_edge,child_seq)
  new_deletions = find_indels_plcs(deletion_plcs,curr_edge,parent_seq)
  insertions = rbind(insertions,new_insertions)
  deletions = rbind(deletions,new_deletions)
}

setwd(paste0(begin_path,"vars/tables_data"))
insertions = data.frame(insertions)
deletions = data.frame(deletions)
save(insertions,file = paste0("insertions",collapse =""))
save(deletions,file = paste0("deletions",collapse =""))
load("insertions")
load("deletions")

add_regions_data = function(data){
  numeric_cols = c("start_plcs","end_plcs","parent_node","child_node")
  data = change_cols_to_numeric(numeric_cols,data)
  missing_amino_acids = missing_seq = missing_codon =
    codon.pos_vals = amino_acid_nums = matrix("0",dim(data)[1],1)
  for (i in 1:dim(data)[1]){
    missing_amino_acids[i] = 
      paste0(regions$ref_amino_acid[data$start_plcs[i]:data$end_plcs[i]],collapse = "")
    missing_seq[i] = 
      paste0(regions$ref_seq[data$start_plcs[i]:data$end_plcs[i]],collapse = "")
    missing_codon[i] = 
      paste0(regions$ref_codon[data$start_plcs[i]:data$end_plcs[i]],collapse = "")
    codon.pos_vals[i] = 
      paste0(regions$codon.pos[data$start_plcs[i]:data$end_plcs[i]],collapse = "")
    amino_acid_nums[i] = 
      paste0(regions$amino_acid_num[data$start_plcs[i]:data$end_plcs[i]],collapse = ".")
  }
  # node_names = nodes_names[order(nodes_names$number),]
  row.names(nodes_names) = nodes_names$number
  data2 = cbind(data[,1:2],
                regions$ref_site[data$start_plcs],
                regions$ref_site[data$end_plcs],
                (as.numeric(data$end_plcs)-as.numeric(data$start_plcs)+1),
                regions$gene[data$start_plcs],regions$gene[data$end_plcs],
                amino_acid_nums,
                data[,3],
                missing_seq,
                # missing_codon,
                missing_amino_acids,
                codon.pos_vals,
                regions$codon.pos[data$start_plcs],regions$codon.pos[data$end_plcs],
                data[,4:5],nodes_names[as.character(data$parent_node),"name"],nodes_names[as.character(data$child_node),"name"])
  
  colnames(data2)[1:dim(data2)[2]] = 
    c("start","end", "ref_start","ref_end","mut_length","gene_start","gene_end","amino_acid_num",
      "mutation","ref_seq","ref_a.acid","codon.pos_vals","codon.pos_start","codon.pos_end",
      "parent_node","child_node","parent_name","child_name")
  return(data2)
}

insertions2 = add_regions_data(insertions)
deletions2 = add_regions_data(deletions)

del_coding_plcs = which(deletions2$gene_start!="0")
deletions2 = deletions2[del_coding_plcs,]
deletions = deletions2[,-c(13:14)]
save(deletions,file = "deletions_coding")

which(insertions2$gene_start!=insertions2$gene_end)

coding_plcs = which(insertions2$gene_start!=0)
insertions2 = insertions2[coding_plcs,]
insertions = insertions2[,-c(13:14)]
save(insertions,file = "insertions_coding")

tmp = deletions2[which(deletions2$gene_start=="S"),]
tmp2 = tmp[order(tmp$end),]
tmp2 = tmp2[order(tmp2$start),]

write.csv(tmp2,file = "deletions_Spike")

tmp = insertions2[which(insertions2$gene_start=="S"),]
tmp2 = tmp[order(tmp$end),]
tmp2 = tmp2[order(tmp2$start),]

write.csv(tmp2,file = "insertions_Spike")

insertions_muts = NULL
insertion_coding_plcs = which(insertions2$gene_start!="0")
for (i in 1:length(insertion_coding_plcs)){
  # print(i)
  insertions_muts = paste0(insertions_muts,
                           insertions2$mutation[insertion_coding_plcs[i]],collapse = "")
}
table(strsplit(insertions_muts,split = "")[[1]])

deletion_muts = NULL
deletion_coding_plcs = which(deletions2$gene_start!="0")
for (i in 1:length(deletion_coding_plcs)){
  #print(i)
  deletion_muts = paste0(deletion_muts,
                         deletions2$mutation[deletion_coding_plcs[i]],collapse = "")
}
table(strsplit(deletion_muts,split = "")[[1]])


##############3
# Check deletions at site 746:
setwd("/Users/keren/Desktop/covid_files/covid_new/trees/")
load("tree_list.746")
plcs = which(deletions_coding$start==746)
# root_node = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
# root_name = paste0(c("node_",root_node),collapse = "")
root_name = "node_32768"

print_lineage = function(node,lineage,count){
  print(node)
  count = count+1
  print(count)
  lineage = rbind(lineage,cbind(tree_list[[node]]@parent_name,tree_list[[node]]@name,
                                tree_list[[node]]@base))
  if ((tree_list[[node]]@children_names)[1]!="_"){
    for (i in 1:length(tree_list[[node]]@children_names)){
      child_name = tree_list[[node]]@children_names[i]
      lineage = print_lineage(child_name,lineage)
    }
  }else{
    return(lineage)
  }
}

lineage = NULL
lineage = print_lineage(root_name,lineage)


