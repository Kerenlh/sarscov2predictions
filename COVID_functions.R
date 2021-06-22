# General functions:
 library("msa")
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("Biostrings")
 BiocManager::install("DECIPHER")

 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("msa")

library("seqinr")
library("ape")
library("DECIPHER")
library("TreeSearch")
library("phangorn")



short_names_tree = function(tree){
  # shortens tree$tip.label
  new_names_tree = NULL
  for (i in 1:length(tree$tip.label)){
    # tmp = strsplit(tree$tip.label[i],split = "\\|")[[1]][1]
    tmp = strsplit(tree$tip.label[i],split = "\\.")[[1]][1]
    tmp2 = strsplit(tmp,"")
    if(tmp2[[1]][length(tmp2[[1]])]=="_"){
      tmp = paste0(tmp2[[1]][1:(length(tmp2[[1]])-1)],collapse = "")
    }
    new_names_tree = c(new_names_tree, tmp)
  }
  tree$tip.label = new_names_tree
  return(tree)
}

short_names_seqs = function(seqs){
  # shortens names_seqs
  new_names_seqs = NULL
  names_seqs = names(seqs)
  for (i in 1:length(names_seqs)){
    # print(i)
    tmp = strsplit(names_seqs[i],split = "\\.")[[1]][1]
    tmp2 = strsplit(tmp,"")
    if(tmp2[[1]][length(tmp2[[1]])]==" "){
      tmp = paste0(tmp2[[1]][1:(length(tmp2[[1]])-1)],collapse = "")
    }
    # print(tmp)
    new_names_seqs = c(new_names_seqs, tmp)
  }
  names(seqs) = new_names_seqs
  return(seqs)
}

short_names_phyDat = function(phyDat_seqs){
  curr_names_seqs = names(phyDat_seqs)
  new_names_seqs = NULL
  for (i in 1:length(phyDat_seqs)){
    new_names_seqs = c(new_names_seqs, strsplit(curr_names_seqs[i],split = "\\ ")[[1]][1])
  }
  names(phyDat_seqs) = new_names_seqs
  return(phyDat_seqs)
}

return_only_tree_seqs = function(seqs,tree){
  relevant_seqs_index = matrix(0,length(tree$tip.label),1)
  missing_seqs = NULL
  for (i in 1:length(tree$tip.label)){
    print(i)
    tmp = which(names(seqs)==tree$tip.label[i])
    if(length(tmp)>0){
      relevant_seqs_index[i,1] = tmp 
    }else{
      missing_seqs = c(missing_seqs,tree$tip.label[i])
    }
  }
  print(missing_seqs)
  new_seqs = seqs[relevant_seqs_index]
  return(new_seqs)
}

#find phyDat_seqs site patterns for msa output
find_site_patterns = function(aligned_seqs){
  site_patterns = NULL
  for (i in 1:length(aligned_seqs@unmasked[[1]])){
    print(i)
    curr_site_pattern = NULL
    for (j in 1:length(aligned_seqs@unmasked)){
      curr_site_pattern = c(curr_site_pattern, as.character(aligned_seqs@unmasked[j][[1]][i]))
    }
    site_patterns = c(site_patterns,paste0(curr_site_pattern, collapse = ""))
  }
  return(site_patterns)
}

change_to_letters = function(sub_col){
  sub_col[which(sub_col==15)] = "-"
  sub_col[which(sub_col==1)] = "A"
  sub_col[which(sub_col==2)] = "C"
  sub_col[which(sub_col==4)] = "G"
  sub_col[which(sub_col==8)] = "T"
  sub_col[which(sub_col==5)] = "R"
  sub_col[which(sub_col==6)] = "GC"
  sub_col[which(sub_col==9)] = "AT"
  sub_col[which(sub_col==10)] = "Y"
  sub_col[which(sub_col==12)] = "TG"
  return(sub_col)
}

change_cols_to_numeric = function(numeric_cols,data){
  for (i in 1:length(numeric_cols)){
    data[,numeric_cols[i]] = as.numeric(data[,numeric_cols[i]])
  }
  return(data)
}

