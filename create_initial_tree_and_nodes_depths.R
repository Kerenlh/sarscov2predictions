setClass(Class='Tree', 
         representation(name= 'character',node = 'numeric', parent_node = 'numeric',
                        parent_name = 'character',  
                        children_nodes = 'matrix', children_names = 'character', 
                        probs = 'matrix', updated_probs = 'matrix', base = 'character',leaf = 'logical'))

######### 
# Create initial tree list:
create_initial_tree_list = function(tree, begin_path){
  find_details_not_leaf = function(curr_name,tree,tree_list,root_node){
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
        # print("This is a leaf identification bug")
      }
      tree_list = find_details_leaf(curr_name,tree,tree_list)
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
          tree_list = find_details_not_leaf(curr_children_names[i],
                                            tree,tree_list,root_node)
        } 
        # all_probs = tree_list[[curr_children_names[i]]]@probs + all_probs
      }
      probs = matrix(0,1,5)
      colnames(probs) = c("A","C","G","T","-")
      updated_probs = matrix(0,1,5)
      base = colnames(probs)[which.max(probs)]
      if(length(which(probs==max(probs)))>1){
        # print(paste(curr_name, probs))
      }
      tt = new("Tree", name = curr_name, node = curr_node, parent_node = curr_parent_node,
               parent_name = curr_parent_name, children_nodes = curr_children_nodes,
               children_names = curr_children_names, base = base, probs = probs,
               updated_probs = updated_probs, leaf = leaf)
      tree_list[[curr_name]] = tt
      return(tree_list)
    }
  }
  
  find_details_leaf = function(curr_name,tree,tree_list){
    #print(curr_name)
    curr_node = which(tree$tip.label==curr_name)
    curr_parent_node = tree$edge[which(tree$edge[,2]==curr_node),1]
    curr_parent_name = tree$tip.label[curr_parent_node]
    if(is.na(curr_parent_name)){
      curr_parent_name = paste0("node_",curr_parent_node,collapse = "")
    }
    curr_children_nodes = tree$edge[which(tree$edge[,1]==curr_node),2]
    if (length(curr_children_nodes)>0 | curr_node>length(tree$tip.label)){
      # print("bug, this is not a leaf!")
    }else{
      base = "A"
      probs = matrix(0,1,5)
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
  
  tree_list = list()
  start_time = Sys.time()
  root_node = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
  root_name = paste0(c("node_",root_node),collapse = "")
  print("Creating initial tree_list (takes about 15 minutes)")
  tree_list = find_details_not_leaf(root_name,tree,tree_list,root_node)
  setwd(paste0(begin_path,"vars/",collapse = ""))
  save(tree_list,file = paste0("tree_list.initial"))
  print(Sys.time()-start_time)
}

######################
# Create node_depth: #
######################

create_node_depth = function(begin_path, tree, tree_list, seqs_names){
  start_time = Sys.time()
  node_depth = data.frame(matrix(0,length(tree_list),2))
  node_depth[,1] = c(1:length(tree_list))
  colnames(node_depth) = c("name","depth")
  cat("\n")
  print("Part 1 out of 4:")
  pb <- txtProgressBar(min = 0, max = length(tree_list), style = 3, width = 50, char = "=")
  for (i in 1:length(tree_list)){
    # print(i)
    node_depth$depth[i] = length(strsplit(tree_list[[i]]@base,split = "")[[1]])
    node_depth$name[i] = tree_list[[i]]@name
    setTxtProgressBar(pb, i)
  }
  rownames(node_depth) = node_depth$name
  node_depth = node_depth[order(node_depth$depth,decreasing = TRUE),]
  
  seqs_names_index = matrix(0,dim(node_depth)[1],1)
  node_depth = cbind(node_depth,seqs_names_index)
  cat("\n")
  print("Part 2 out of 4:")
  pb <- txtProgressBar(min = 0, max = dim(node_depth)[1], style = 3, width = 50, char = "=")
  for (i in 1:dim(node_depth)[1]){
    tmp = which(seqs_names==node_depth$name[i])
    if (length(tmp)>0){
      node_depth$seqs_names_index[i] = tmp
    }
    setTxtProgressBar(pb, i)
  }
  parent = children = matrix("m",dim(node_depth)[1],1)
  node_depth = cbind(node_depth,parent,children)
  cat("\n")
  print("Part 3 out of 4:")
  pb <- txtProgressBar(min = 0, max = dim(node_depth)[1], style = 3, width = 50, char = "=")
  for (i in 1:dim(node_depth)[1]){
    # print(i)
    node_depth$parent[i] = tree_list[[node_depth$name[i]]]@parent_name
    curr_children = tree_list[[node_depth$name[i]]]@children_names
    if(curr_children=="_"){
      node_depth$children[i] = NA
    }else{
      node_depth$children[i] = paste0(curr_children,collapse = "*")
    }
    setTxtProgressBar(pb, i)
  }
  num_children = matrix(0,dim(node_depth)[1],1)
  cat("\n")
  print("Part 4 out of 4:")
  pb <- txtProgressBar(min = 0, max = dim(node_depth)[1], style = 3, width = 50, char = "=")
  for (i in 1:dim(node_depth)[1]){
    if (is.na(node_depth$children[i])==FALSE){
      num_children[i] = length(strsplit(node_depth$children[i],split = "\\*")[[1]])
      # print(i)
      # print(num_children[i])
    }
    setTxtProgressBar(pb, i)
  }
  node_depth = cbind(node_depth,num_children)
  save(node_depth,file = "node_depth")
  print(Sys.time()-start_time)
}