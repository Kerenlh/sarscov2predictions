begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/reconstructed_seqs"))
# load("Ppart_seqs.1")
# tmp = apply(part_seqs,1,paste0,collapse = "")
# all_seqs = tmp
# rm(tmp)
# rm(part_seqs)
load("all_seqs")
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/",collapse = ""))
load("all_aligned_seqs")
tree = read.tree("tree_18_2.nwk")
tree2 = short_names_tree(tree)
# load("codons_table")
load("site_details")
load("tree_list.1_trimmed")
# setwd("/Users/keren/Dropbox/covid/vars/big_trees/")

# Explore differences between the leaf sequences in all_seqs and 
# the aligned sequences
bad_plcs = NULL
count = 0
orig_seq = my_seq = sites = NULL
for (i in 1:length(tree_list)){
  if (tree_list[[i]]@leaf==TRUE){
    print(i)
    count = count+1
    name = tree_list[[i]]@name
    tmp1 = paste0(aligned_seqs[name][[1]][1:length(aligned_seqs[1][[1]])])
    tmp1 = strsplit(tmp1,split = "")[[1]]
    tmp2 = all_seqs[name][[1]]
    tmp2 = strsplit(tmp2,split = "")[[1]]
    plcs = which((tmp1!=tmp2)==TRUE)
    # tmp1[plcs]
    # tmp2[plcs]
    # plcs
    if(length(plcs)>0){
      bad_plcs = c(bad_plcs,i)
      orig_seq = c(orig_seq,tmp1[plcs])
      my_seq = c(my_seq,tmp2[plcs])
      sites = c(sites,plcs)
    }
  }
}

###########33
# check deletions at site 746
setwd("/Users/keren/Dropbox/covid/vars/")
load("diffs_count")
diffs_count_old = diffs_count
setwd(paste0(begin_path,"vars/tables_data"))
load("deletions_coding")
plcs = which(deletions$start==746)
rel_seqs = deletions$child_name[plcs]
diffs_count_names = colnames(diffs_count)
rel_names = rel_seqs[which(is.element(rel_seqs,diffs_count_names))]

rel_distances = NULL
for (i in 1:(length(rel_names)-1)){
  for (j in (i+1):length(rel_names)){
    rel_distances = rbind(rel_distances,
                          c(rel_names[i],rel_names[j],
                            max(diffs_count[rel_names[i],rel_names[j]],
                                diffs_count[rel_names[j],rel_names[i]])))
  }
}
colnames(rel_distances) = c("seq_1","seq_2","distance")
rel_distances = data.frame(rel_distances)
rel_distances$distance = as.numeric(rel_distances$distance)
summary(rel_distances$distance)

plcs = which(diffs_count!=-1)
summary(diffs_count[plcs])

setwd("/Users/keren/Dropbox/covid/new/vars/site_patterns/")
load("all_site_patterns")
load("all_site_patterns_plcs")
load("all_split_seqs_1_10000")

plcs = which(codons_table$site>=746 & codons_table$site<=754)
codons_table[plcs,]

plcs = which(deletions$start==746)
rel_seqs = deletions$child_name[plcs]
rel_names = rel_seqs[which(is.element(rel_seqs,tree2$tip.label))]
rel_long_names = tree$tip.label[which(is.element(tree2$tip.label,rel_seqs))]
split_leafs = list()
for (i in 1:length(rel_names)){
  print(i)
  split_leafs[[i]] = strsplit(all_tree_list[[rel_names[i]]]@probs,split = "")[[1]]
}

diffs_count = matrix(-1,length(rel_names),length(rel_names))
coding_plcs = which(regions$gene!=0)
for (i in 1:(length(rel_names)-1)){
  print(i)
  start_time = Sys.time()
  for (j in (i+1):length(rel_names)){
    diffs_count[i,j] = length(which(split_leafs[[i]][coding_plcs]!=split_leafs[[j]][coding_plcs]))
  }
  print(Sys.time()-start_time)
}

colnames(diffs_count) = row.names(diffs_count) = rel_names
plcs = which(diffs_count!=-1)
summary(diffs_count[plcs])
min_plcs = which(diffs_count==1,arr.ind = TRUE)
row.names(diffs_count)[min_plcs[1]]
colnames(diffs_count)[min_plcs[2]]

diffs_count_all = matrix(-1,length(rel_names),length(rel_names))
for (i in 1:(length(rel_names)-1)){
  print(i)
  start_time = Sys.time()
  for (j in (i+1):length(rel_names)){
    diffs_count_all[i,j] = length(which(split_leafs[[i]]!=split_leafs[[j]]))
  }
  print(Sys.time()-start_time)
}

colnames(diffs_count_all) = row.names(diffs_count_all) = rel_names
plcs = which(diffs_count_all!=-1)
summary(diffs_count_all[plcs])
min_plcs = which(diffs_count_all==1,arr.ind = TRUE)


############
# Explore leaves that have the same parent:
rel_nodes = NULL
for (i in 1:length(all_tree_list)){
  if (length(all_tree_list[[i]]@children_names)>1){
    count = 0
    for (j in 1:length(all_tree_list[[i]]@children_names)){
      if (strsplit(all_tree_list[[i]]@children_names[j],split="")[[1]][1]!="n"){
        count = count+1
      }
    }
    if (count>=2){
      rel_nodes = c(rel_nodes,i)
    }
  }
}

sites = NULL
subs = list()
for (i in rel_nodes){
  # print(i)
  child_seqs = NULL
  for (j in 1:length(all_tree_list[[i]]@children_names)){
    if (strsplit(all_tree_list[[i]]@children_names[j],split="")[[1]][1]!="n"){
      child = all_tree_list[[i]]@children_names[j]
      child_seqs = rbind(child_seqs,strsplit(all_tree_list[[child]]@probs,split = "")[[1]])
    }
  }
  tmp =  apply(child_seqs,2,unique)
  for (k in 1:length(tmp)){
    if (length(tmp[k][[1]])>1){
      sites = c(sites,k)
      if (k==29744){
        print(i)
        print(tmp[k])
      }
      count = length(sites)
      subs[[count]] = tmp[k][[1]]
    }
  }
}


plcs = NULL
for (i in 1:length(subs)){
  if(length(subs[[i]])==2 & is.element("-",subs[[i]])==FALSE){
    plcs = c(plcs,i)
  }
  if(length(subs[[i]])>2){
    plcs = c(plcs,i)
  }
}

tmp = sort(table(sites),decreasing = TRUE)

setwd("/Users/keren/Dropbox/covid/vars/")
load("codons_table1")
load("base_table1")

codon_sites = unique(codons_table$site)
sites2 = sites[which(is.element(sites,codon_sites))]
tmp= sort(table(sites2),decreasing = TRUE)

tmp2 = codons_table[order((codons_table$y+codons_table$line),decreasing = TRUE),]


####################
# Check ancient subs:

root_edge = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
#root_edge = 12225
root_name = paste0(c("no_name_",root_edge),collapse = "")

compare_seqs = function(name1,name2){
  tmp2 = all_tree_list[[name2]]@probs
  tmp2 = strsplit(tmp2,split = "")[[1]]
  # children_names = all_tree_list[[root_name]]@children_names
  tmp1 = all_tree_list[[name1]]@probs
  tmp1 = strsplit(tmp1,split = "")[[1]]
  plcs = which((tmp1!=tmp2)==TRUE)
  print(c("name1",tmp1[plcs]))
  print(c("name2",tmp2[plcs]))
  print(plcs)
}

name1 = "no_name_20106"
name2 = root_name


compare_seqs("no_name_20106",all_tree_list[["no_name_20106"]]@children_names[2])
which(codons_table$site==21371)


#######3
setwd("/Users/keren/Dropbox/covid/vars/")
load("all_site_patterns")
unique_site_patterns = unique(all_site_patterns)
entropy = NULL
for(i in 1:length(unique_site_patterns)){
  print(i)
  tmp = strsplit(unique_site_patterns[[i]],split = "")[[1]]
  tmp2 = table(tmp)/length(tmp)
  entropy = c(entropy, -sum(tmp2*log(tmp2)))
}

tmp3 = order(entropy,decreasing = TRUE)[1:100]

i = 9781
site = regions$site[which(regions$site_pattern==i)]
site = regions$site[which(regions$site_pattern==tmp3)]
base_table[which(base_table$site==site),]
tmp = strsplit(unique_site_patterns[[i]],split = "")[[1]]
tmp2 = table(tmp)/length(tmp)

#######_23_11_20
# Check Ghana seq MT890240.1
name = "MT890240.1"
tmp1 = paste0(aligned_seqs[name][[1]][1:30139])
Ghana_seq = strsplit(tmp1,split = "")[[1]]
name3 = "MT825091.1"
tmp3 = paste0(aligned_seqs[name3][[1]][1:30139])
long_branch_seq = strsplit(tmp3,split = "")[[1]]

get_splitted_seq = function(name){
  tmp = paste0(aligned_seqs[name][[1]][1:30139])
  splitted_seq = strsplit(tmp,split = "")[[1]]
  return(splitted_seq)
}

Ghana_seq = get_splitted_seq("MT890240.1")
long_branch_seq = get_splitted_seq("MT825091.1")
diffs_count = NULL
for (i in 1:length(all_tree_list)){
  if (all_tree_list[[i]]@leaf==TRUE){
    print(i)
    tmp = all_tree_list[[i]]@probs
    tmp = strsplit(tmp,split = "")[[1]]
    diffs_count = rbind(diffs_count, c(all_tree_list[[i]]@name, length(which(tmp!=Ghana_seq)),
                                       length(which(tmp!=long_branch_seq))))
  }
}

tmp = strsplit(all_tree_list[["MT890240.1"]]@probs,split = "")[[1]]
plcs = which(tmp!=Ghana_seq)


parent_name = "no_name_18149"

name1 = "MT907519.1"
name2 = "MT801000.1"
tmp1 = paste0(aligned_seqs[name1][[1]][1:30139])
tmp1 = strsplit(tmp1,split = "")[[1]]
tmp2 = paste0(aligned_seqs[name2][[1]][1:30139])
tmp2 = strsplit(tmp2,split = "")[[1]]
plcs = which(tmp1!=tmp2)

# Long branch tip = "MT825091.1"
name3 = "MT825091.1"
tmp3 = paste0(aligned_seqs[name3][[1]][1:30139])
tmp3 = strsplit(tmp3,split = "")[[1]]
plcs = which(tmp3!=tmp2)

###########3
# Check diffreneces between seqs:
leaf_names = NULL
split_leafs = list()
for (i in 1:length(all_tree_list)){
  print(i)
  if (all_tree_list[[i]]@leaf==TRUE){
    leaf_names = rbind(leaf_names,c(i,all_tree_list[[i]]@name))
  }
}

leaf_plcs = as.numeric(leaf_names[,1])
split_leafs = list()
for (i in 1:length(leaf_plcs)){
  print(i)
  split_leafs[[i]] = strsplit(all_tree_list[[leaf_plcs[i]]]@probs,split = "")[[1]]
}

setwd("/Users/keren/Dropbox/covid/vars/")
# save(split_leafs,file = "split_leafs")
load("split_leafs")
diffs_count = matrix(-1,length(leaf_plcs),length(leaf_plcs))
coding_plcs = which(regions$gene!=0)
for (i in 1:length(leaf_plcs)){
  print(i)
  start_time = Sys.time()
  for (j in (i+1):length(leaf_plcs)){
    diffs_count[i,j] = length(which(split_leafs[[i]]!=split_leafs[[j]]))
  }
  print(Sys.time()-start_time)
}
colnames(diffs_count) = leaf_names[,2]
rownames(diffs_count) = leaf_names[,2]

# save(diffs_count,file = "diffs_count")
setwd("/Users/keren/Dropbox/covid/vars/")
load("diffs_count")

all_tree_list = gene_tree
siblings_distances = NULL
for (i in 1:length(all_tree_list)){
  print(i)
  if (all_tree_list[[i]]@leaf==TRUE){
    name = all_tree_list[[i]]@name
    parent_name = all_tree_list[[i]]@parent_name
    siblings_names = all_tree_list[[parent_name]]@children_names
    if (length(siblings_names)>=2){
      plcs = which(siblings_names!=name)
      k = which(colnames(diffs_count)==name)
      if(k<dim(diffs_count)[1]){
        closest_seq_plc_col = which.min(diffs_count[name,(k+1):dim(diffs_count)[2]])
      }else{
        closest_seq_plc_col = NULL
      }
      if(k>1){
        closest_seq_plc_row = which.min(diffs_count[1:(k-1),name])
      }else{
        closest_seq_plc_row = NULL
      }
      if(is.null(closest_seq_plc_col)){
        closest_seq_name = colnames(diffs_count)[closest_seq_plc_row]
      }else if (is.null(closest_seq_plc_row)){
        closest_seq_name = colnames(diffs_count)[closest_seq_plc_col+k]
      }else{
        tmp = which.min(c(diffs_count[k, (k+closest_seq_plc_col)],diffs_count[closest_seq_plc_row,k]))
        if (tmp==1){
          closest_seq_name = colnames(diffs_count)[closest_seq_plc_col+k]
        }else{
          closest_seq_name = colnames(diffs_count)[closest_seq_plc_row]
        }
      }
      while(length(plcs)>=1){
        if (paste0(strsplit(siblings_names[plcs[1]],split = "")[[1]][1:8],collapse = "")==
            "no_name_"){
          plcs = plcs[-1]
        }else{
          break
        }
      }
      if(length(plcs)>=1){
        siblings_distances = 
          rbind(siblings_distances,
                c(name,siblings_names[plcs[1]],
                  max(diffs_count[name,siblings_names[plcs[1]]],
                      diffs_count[siblings_names[plcs[1]],name]),
                  closest_seq_name,
                  max(diffs_count[name,closest_seq_name],
                      diffs_count[closest_seq_name,name])))    
        
      }
    }
  }
}

colnames(siblings_distances) = c("seq","sibling_seq","diffs_from_sibling",
                                 "closest_seq","diffs_from_closest_seq")
siblings_distances = data.frame(siblings_distances)
siblings_distances$diffs_from_sibling = 
  as.numeric(siblings_distances$diffs_from_sibling)
siblings_distances$diffs_from_closest_seq= 
  as.numeric(siblings_distances$diffs_from_closest_seq)
tmp = siblings_distances$diffs_from_sibling- siblings_distances$diffs_from_closest_seq
distances_comparison = NULL
for (i in 1:dim(siblings_distances2)[1]){
  plc = which(siblings_distances$seq==siblings_distances2$seq[i])
  if(length(plc)>0){
    distances_comparison = rbind(distances_comparison,
                                 c(siblings_distances[plc,],siblings_distances2[i,2:3]))
  }
}
colnames(distances_comparison)[6:7] = c("gene_sibling_seq","diffs_from_gene_sibling")
distances_comparison = apply(distances_comparison,2,unlist)
distances_comparison = data.frame(distances_comparison)
distances_comparison$diffs_from_closest_seq = 
  as.numeric(distances_comparison$diffs_from_closest_seq)
distances_comparison$diffs_from_sibling = 
  as.numeric(distances_comparison$diffs_from_sibling)
distances_comparison$diffs_from_gene_sibling = 
  as.numeric(distances_comparison$diffs_from_gene_sibling)
tmp = (distances_comparison$diffs_from_sibling-distances_comparison$diffs_from_gene_sibling)
distances_comparison[which.min(tmp),]
distances_comparison[which.max(tmp),]
summary(siblings_distances$diffs_from_closest_seq)
summary(siblings_distances$diffs_from_sibling)

siblings_distances_gene_tree = siblings_distances2
setwd("/Users/keren/Dropbox/covid/check_tree/")
save(distances_comparison,file = "distances_comparison")
save(siblings_distances,file = "sibling_distances")
save(siblings_distances_gene_tree,file = "siblings_distances_gene_tree")
load("distances_comparison")
load("sibling_distances")
load("siblings_distances_gene_tree")

coding_plcs = which(regions$gene!=0)
which.min(tmp)
distances_comparison[634,]
tmp1 = strsplit(paste0(all_tree_list[["MT079843.1"]]@probs),split = "")[[1]]
tmp2 = strsplit(paste0(all_tree_list[["MT079853.1"]]@probs),split = "")[[1]]
plcs = which(tmp1!=tmp2)

distances_comparison[3088,]
tmp1 = strsplit(paste0(all_tree_list[["MT786843.1"]]@probs),split = "")[[1]]
tmp2 = strsplit(paste0(all_tree_list[["MT375449.1"]]@probs),split = "")[[1]]
plcs = which(tmp1!=tmp2)

distances_comparison[7756,]
tmp1 = strsplit(paste0(all_tree_list[["MT786843.1"]]@probs),split = "")[[1]]
tmp2 = strsplit(paste0(all_tree_list[["MT375449.1"]]@probs),split = "")[[1]]
plcs = which(tmp1!=tmp2)


get_ancestors = function(all_tree_list,name,root_name){
  ancestors = NULL
  curr_name = name
  while(all_tree_list[[curr_name]]@name!=root_name){
    ancestors = c(ancestors,all_tree_list[[curr_name]]@parent_name)
    curr_name = all_tree_list[[curr_name]]@parent_name
    # print(curr_name)
  }
  return(ancestors)
}

root_name = "no_name_12225"
tmp = get_ancestors(all_tree_list,name = "MT831451.1",root_name = root_name)
tmp2 = get_ancestors(all_tree_list,name = "MT612185.1",root_name = root_name)
both = intersect(tmp,tmp2)
which(tmp==both[1])
which(tmp2==both[1])
num_seperating_nodes = which(tmp==both[1])+which(tmp2==both[1])
