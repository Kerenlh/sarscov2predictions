get_node_names = function(tree_list){
  nodes_names = NULL
  for (i in 1:length(tree_list)){
    # print(i)
    nodes_names = rbind(nodes_names,cbind(tree_list[[i]]@name, tree_list[[i]]@node))
  }
  colnames(nodes_names) = c("name","number")
  nodes_names = data.frame(nodes_names)
  nodes_names$number = as.numeric(nodes_names$number)
  return(nodes_names)
}

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

create_codon_states_and_outputs = function(tree_list,regions,tree,all_seqs,node_depth,leaves_flag){
  leaves_names = node_depth$name[which(node_depth$num_children==0)]
  nodes_names = get_node_names(tree_list)
  codon_plcs = which(regions$codon.pos!=0)
  codon_states = "tmp"
  init_outputs = matrix(0,dim(regions)[1],10)
  colnames(init_outputs) = c("A","C","G","T","line","q","s","M","exposure","zero_branch")
  init_outputs = data.frame(init_outputs)
  codon_outputs = init_outputs[1,]
  edges_with_no_subs = NULL
  start_time = Sys.time()
  if (leaves_flag){
    num_iter = length(leaves_names)
  }else{
    num_iter = dim(tree$edge)[1]
  }
  cat("\n")
  pb <- txtProgressBar(min = 0, max = num_iter, style = 3, width = 50, char = "=")
  cat("\n")
  for (i in 1:num_iter){
    curr_outputs = init_outputs
    if(leaves_flag){
      curr_name = leaves_names[i]	
      edge_length = 1	
      curr_seq = strsplit(all_seqs[curr_name][[1]],split = "")[[1]]
      parent_seq = curr_seq
    }else{
      curr_edge = tree$edge[i,]
      parent_name = nodes_names$name[which(nodes_names$number==curr_edge[1])]
      child_name = nodes_names$name[which(nodes_names$number==curr_edge[2])]
      edge_length = tree$edge.length[i]
      parent_seq = strsplit(all_seqs[parent_name][[1]],split = "")[[1]]
      child_seq = strsplit(all_seqs[child_name][[1]],split = "")[[1]]
      subs_plcs = which((parent_seq==child_seq)==FALSE)
      if(length(subs_plcs)>0){
        for (j in 1:length(subs_plcs)){
          if (child_seq[subs_plcs[j]]=="-"){
            curr_outputs[subs_plcs[j],"line"]=1
          }else if (child_seq[subs_plcs[j]]=="?"){
            curr_outputs[subs_plcs[j],"q"]=1
          }else if (child_seq[subs_plcs[j]]=="*"){
            curr_outputs[subs_plcs[j],"s"]=1
          }else{
            curr_outputs[subs_plcs[j],child_seq[subs_plcs[j]]] = 1
          }
        }
      }else{
        edges_with_no_subs = c(edges_with_no_subs,i)
      }
    }
    curr_outputs$exposure = edge_length
    if (edge_length==0){
      curr_outputs$zero_branch = 1
    }
    
    gaps_plcs = which(parent_seq=="-")
    short_parent_sites = c(1:length(parent_seq))
    if (length(gaps_plcs)>0){
      short_parent = parent_seq[-gaps_plcs]
      short_parent_sites = short_parent_sites[-gaps_plcs]
      tmp_check = which(is.element(codon_plcs,gaps_plcs))
      if (length(tmp_check)>0){
        curr_codon_plcs = codon_plcs[-which(is.element(codon_plcs,gaps_plcs))]
      }else{
        curr_codon_plcs = codon_plcs
      }
      curr_sites = regions$site[-gaps_plcs]
    }else{
      short_parent = parent_seq
      curr_codon_plcs = codon_plcs
      curr_sites = regions$site
    }
    short_parent = cbind(short_parent_sites,short_parent)
    colnames(short_parent) = c("sites","seq")
    short_parent = data.frame(short_parent)
    short_parent$sites = as.numeric(short_parent$sites)
    codon_short_parent = parent_seq[curr_codon_plcs]
    curr_codon_sites = regions$site[curr_codon_plcs]
    curr_codon.pos = regions$codon.pos[curr_codon_plcs]
    codons = NULL
    k=1
    while((k+2) <=length(curr_codon_sites)){
      if (sum(curr_codon.pos[k:(k+2)]==c(1,2,3))==3){
        codons = c(codons,rep(paste0(codon_short_parent[k:(k+2)],collapse = ""),3))
        k = k+3
      }else{
        codons = c(codons,"missing")
        # print(paste("k=",k,"i=",i))
        k= k+1
      }
    }
    L.N. = c("*",short_parent$seq[1:(length(short_parent$seq)-1)])
    R.N. = c(short_parent$seq[2:length(short_parent$seq)],"*")
    L.N.codon = short_parent$seq[which(is.element(short_parent$sites,curr_codon_sites)==TRUE)-1]
    R.N.codon = short_parent$seq[which(is.element(short_parent$sites,curr_codon_sites)==TRUE)+1]
    
    curr_codon_states = paste0(L.N.codon,R.N.codon,
                               parent_seq[curr_codon_plcs],codons,
                               regions$site[curr_codon_plcs])
    codon_uncertainty_plcs = which(parent_seq[curr_codon_plcs]=="?" | 
                                     curr_outputs$q[curr_codon_plcs]=="1")
    tmp = update_states_and_outputs(curr_states = curr_codon_states,
                                    states = codon_states,
                                    curr_outputs = curr_outputs[curr_codon_plcs,],
                                    outputs = codon_outputs,
                                    uncertainty_plcs = codon_uncertainty_plcs)
    codon_states = tmp$states
    codon_outputs = tmp$outputs
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  print(Sys.time()-start_time)
  codon_states = codon_states[-1]
  codon_outputs = codon_outputs[-1,]
  
  setwd(paste0(begin_path,"vars/"))
  if(leaves_flag){
    save(codon_states,file = "codon_states_leaves")
    save(codon_outputs,file = "codon_outputs_leaves")
  }else{
    save(codon_states,file = "codon_states")
    save(codon_outputs,file = "codon_outputs")
  }
}