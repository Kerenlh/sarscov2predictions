begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/old_fitch_trees"))

file_names = dir()
for (i in 1:length(file_names)){
  load(file_names[i])
  tree_list_old_fitch = tree_list
  setwd("/Users/keren/Desktop/covid_files/new2/Lanfear trees/")
  load(file_names[i])
  bases = cbind(bases,matrix("0",length(bases),1))
  colnames(bases) = c("recontruct_tree","old_fitch")
  bases = data.frame(bases)
  bases_rows = row.names(bases)
  for (j in 1:dim(bases)[1]){
    print(j)
    bases$old_fitch[j] = tree_list[[bases_rows[[j]]]]@base
  }
  setwd(paste0(begin_path,"vars/old_fitch_trees"))
  save(bases,file = paste0("bases_comaprison",file_names[i]))
  
}

get_muts = function(curr_bases,node_depth){
  parents = node_depth$parent
  mut_plcs = (which(curr_bases!=curr_bases[parents]))
  if (length(mut_plcs)>0){
    for (i in 1:length(mut_plcs)){
      curr_parent = parents[mut_plcs[i]]
      curr_children = strsplit(node_depth$children[mut_plcs[i]],split = "\\*")[[1]]
      if (is.na(curr_children[1])==FALSE){ # If there are no children don't take majority, there was already an epsilon correction from the father
        tmp = table(c(curr_bases[curr_parent],curr_bases[curr_children]))
        if(length(which(tmp==max(tmp)))==1){
          # if there is more than 1 maximum do nothing
          curr_bases[mut_plcs[i]] = names(which.max(tmp))
        } 
      }
    }
    
    mut_plcs = (which(curr_bases!=curr_bases[parents]))
  }
  # check # of mutations:
  if(length(mut_plcs)>0){
    muts = cbind(curr_bases[parents][mut_plcs],curr_bases[mut_plcs],parents[mut_plcs],names(curr_bases)[mut_plcs])
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
      tmp = table(curr_bases[curr_children])
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
    return(list(muts,problematic_muts))
  }
}

setwd(paste0(begin_path,"vars/"))
load("node_depth")
setwd(paste0(begin_path,"vars/old_fitch_trees"))
file_names = dir()
for (i in 1:11){
  print(file_names[i])
  load(file_names[i])
  which(bases[,1]!=bases[,2])
  curr_bases = bases$recontruct_tree
  names(curr_bases) = row.names(bases)
  tmp = get_muts(bases = curr_bases,node_depth = node_depth)
  muts.reconstruct_tree = tmp[[1]]
  problematic_muts.reconstruct_tree = tmp[[2]]
  print(dim(muts.reconstruct_tree)[1])
  curr_bases = bases$old_fitch
  names(curr_bases) = row.names(bases)
  tmp = get_muts(bases = curr_bases,node_depth = node_depth)
  muts.old_fitch = tmp[[1]]
  problematic_muts.old_fitch = tmp[[2]]
  print(dim(muts.old_fitch)[1])
  if(dim(muts.old_fitch)[1]<dim(muts.reconstruct_tree)[1]){
    print(paste("more muts in new alg, file_name = ",file_names[i],"i = ",i))
  }
  if(length(problematic_muts.old_fitch)[1]<length(problematic_muts.reconstruct_tree)[1]){
    print(paste("more problematic muts in new alg, file_name = ",file_names[i],"i = ",i))
  }
  
   # site_num = strsplit(file_names[i],split = "\\.")[[1]][2]
  # setwd(paste0(begin_path,"vars/muts"))
  # load("muts.11.4069")
  # load(paste0("muts.",site_num))
  
  
}