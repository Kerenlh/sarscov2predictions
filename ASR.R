# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biostrings")

ASR = function(begin_path,site_patterns_numbers,node_depth,unique_site_patterns,tree_list){
  
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
  
  curr_nums = site_patterns_numbers
  for (k in 1:length(curr_nums)){
    start_time = Sys.time()
    # print(k)
    # print(curr_nums[k])
    probs = data.frame(matrix("m",dim(node_depth)[1],5))
    row.names(probs) = row.names(node_depth)
    colnames(probs) = c("A","C","G","T","-")
    j = which(site_patterns_numbers==curr_nums[k])
    # print(j)
    curr_site_pattern = strsplit(unique_site_patterns[j],split = "")[[1]]
    curr_bases = curr_site_pattern[node_depth$seqs_names_index[leaves_plcs]]
    probs[leaves_plcs,] = get_probs_from_bases(curr_bases)
    # check that node_depth is sorted by depth, decreasing = T
    print(paste("Site pattern number:", curr_nums[k]))
    cat("\n")
    print("Going over inner nodes:")
    pb <- txtProgressBar(min = 0, max = length(inner_nodes_plcs), style = 3, width = 50, char = "=")
    for (i in 1:length(inner_nodes_plcs)){
      # print(i)
      children = tree_list[[node_depth$name[inner_nodes_plcs[i]]]]@children_names
      probs[inner_nodes_plcs[i],] = apply(matrix(as.numeric(unlist(probs[children,])),length(children),5),2,mean)
      setTxtProgressBar(pb, i)
    }
    
    # print(Sys.time()-start_time)
    probs = matrix(as.numeric(unlist(probs)),dim(node_depth)[1],5)
    colnames(probs) = c("A","C","G","T","-")
    row.names(probs) = row.names(node_depth)
    
    cat("\n")
    print("Going from the root down:")
    pb <- txtProgressBar(min = 0, max = length(inner_nodes_plcs), style = 3, width = 50, char = "=")
    for (i in (length(inner_nodes_plcs)-1):1){ # going from the root (shortest depth to the leaves)
      # The root itself is not updated because it doesn't have a parent
      parent = node_depth$parent[inner_nodes_plcs[i]]
      epsilon = 1/(node_depth$num_children[inner_nodes_plcs[i]])
      probs[inner_nodes_plcs[i],] = (probs[inner_nodes_plcs[i],]+epsilon*probs[parent,])/(1+epsilon)  
      # print(i)
      setTxtProgressBar(pb, (length(inner_nodes_plcs)-i))
    }
    # The update order in the leaves doesn't matter
    parents = node_depth$parent[leaves_plcs]
    probs[leaves_plcs,] = (probs[leaves_plcs,]+probs[parents,])/2
    
    bases = colnames(probs)[apply(probs,1,which.max)]
    names(bases) = row.names(probs)
    
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
    # print(Sys.time()-start_time)
    # check # of mutations:
    if(length(mut_plcs)>0){
      muts = cbind(bases[parents][mut_plcs],bases[mut_plcs],parents[mut_plcs],names(bases)[mut_plcs])
      muts = data.frame(muts)
      colnames(muts) = c("before","after","parent","base")
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
      problematic_muts = NULL
      for (i in 1:length(mut_plcs)){
        if (children_bases[i,muts$before[i]]>0){
          problematic_muts = c(problematic_muts,i)
        }
      }
      problematic_muts = c(problematic_muts,which(is.element(muts$parent,muts$base)))
      problematic_muts = c(problematic_muts,which(is.element(muts$base,muts$parent)))
    }
    setwd(paste0(begin_path,"vars/trees/",collapse = ""))
    save(bases,file = paste0("tree_list.",curr_nums[k],collapse=""))
    if (length(mut_plcs)>0){
      setwd(paste0(begin_path,"vars/muts/",collapse = ""))
      save(muts,file = paste0("muts.",curr_nums[k],".",dim(muts)[1],collapse = ""))
      if(length(problematic_muts)>0){
        problematic_muts_mat = muts[problematic_muts,]
        save(problematic_muts_mat,file = paste0("problematic_muts.",curr_nums[k],".",length(problematic_muts),collapse = ""))
      }
    }
    end_time = Sys.time()
    # print(end_time-start_time)
    cat("\n")
    print(paste0("tree_list.",curr_nums[k],collapse=""))
    # print(k)
  }
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
