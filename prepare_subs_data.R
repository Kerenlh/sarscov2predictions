####
# remove outlier branch!

for (p in 8){
parallel_nums = p
  # begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
  begin_path = "/Users/keren/Dropbox/covid/new2/"
  # begin_path = "C:/Users/Keren/Dropbox/covid/new2/"
  cluster_delta = 7500
  
  # get_node_names = function(tree_list){
  #   nodes_names = NULL
  #   for (i in 1:length(tree_list)){
  #     print(i)
  #     nodes_names = rbind(nodes_names,cbind(tree_list[[i]]@name, tree_list[[i]]@node))
  #   }
  #   colnames(nodes_names) = c("name","number")
  #   nodes_names = data.frame(nodes_names)
  #   nodes_names$number = as.numeric(nodes_names$number)
  #   return(nodes_names)
  # }
  # setwd(paste0(begin_path,"vars/"))
  # load("tree_list.1")
  # nodes_names = get_node_names(tree_list)
  # setwd(paste0(begin_path,"vars/"))
  # save(nodes_names,file = "nodes_names")
  
  setwd(begin_path)
  source("COVID_functions.R")
  setwd(paste0(begin_path,"vars/"))
  load("nodes_names")
  load("site_details")
  tree = read.tree("ft_SH.tree")
  setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  load("all_seqs")
  
  codon_plcs = which(regions$codon.pos!=0)
  
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
  ##########
  
  base_states = "tmp"
  codon_states = "tmp"
  init_outputs = matrix(0,dim(regions)[1],10)
  colnames(init_outputs) = c("A","C","G","T","line","q","s","M","exposure","zero_branch")
  init_outputs = data.frame(init_outputs)
  base_outputs = init_outputs[1,]
  codon_outputs = init_outputs[1,]
  edges_with_no_subs = NULL
  start_time = Sys.time()
  for (i in c((1+(parallel_nums-1)*cluster_delta):
              min((parallel_nums*cluster_delta),dim(tree$edge)[1]))){
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
    curr_outputs = init_outputs
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
        print(paste("k=",k,"i=",i))
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
    
    curr_base_states = paste0(L.N.,R.N.,parent_seq,regions$site)
    base_uncertainty_plcs = which(parent_seq=="?" | 
                                    curr_outputs$q=="1")
    tmp = update_states_and_outputs(curr_states = curr_base_states,
                                    states = base_states,
                                    curr_outputs = curr_outputs[regions$site,],
                                    outputs = base_outputs,
                                    uncertainty_plcs = base_uncertainty_plcs)
    base_states = tmp$states
    base_outputs = tmp$outputs
    # end_time = Sys.time()
    # print(end_time-start_time)
    
  }
  print(Sys.time()-start_time)
  base_outputs = base_outputs[-1,]
  base_states = base_states[-1]
  codon_states = codon_states[-1]
  codon_outputs = codon_outputs[-1,]
  
  setwd(paste0(begin_path,"vars/tables_data"))
  save(base_states,file = paste0("base_states_",parallel_nums,collapse =""))
  save(codon_states,file = paste0("codon_states_",parallel_nums,collapse =""))
  save(base_outputs,file =  paste0("base_outputs_",parallel_nums,collapse =""))
  save(codon_outputs,file = paste0("codon_outputs_",parallel_nums,collapse =""))
}
# ###########3
# # Unite states and outputs:
begin_path = "/Users/keren/Dropbox/covid/new2/"
# setwd(paste0(begin_path,"vars/"))
# load("site_details")
setwd(paste0(begin_path,"vars/tables_data/"))
file_names = dir()

load("base_states_1")
load("base_outputs_1")
all_base_states = base_states
all_base_outputs = base_outputs
for (i in 2:11){
  load(paste0("base_states_",i))
  load(paste0("base_outputs_",i))
  tmp = update_states_and_outputs(curr_states = base_states,
                                  states = all_base_states,
                                  curr_outputs = base_outputs,
                                  outputs = all_base_outputs,
                                  uncertainty_plcs = NULL)
  all_base_states = tmp$states
  all_base_outputs = tmp$outputs
}

base_outputs = all_base_outputs
base_states = all_base_states
save(base_states,file = "base_states")
save(base_outputs,file =  "base_outputs")

load("codon_outputs_1")
load("codon_states_1")
all_codon_states = codon_states
all_codon_outputs = codon_outputs
for (i in 2:11){
  load(paste0("codon_states_",i))
  load(paste0("codon_outputs_",i))
  tmp = update_states_and_outputs(curr_states = codon_states,
                                  states = all_codon_states,
                                  curr_outputs = codon_outputs,
                                  outputs = all_codon_outputs,
                                  uncertainty_plcs = NULL)
  all_codon_states = tmp$states
  all_codon_outputs = tmp$outputs
}
codon_outputs = all_codon_outputs
codon_states = all_codon_states
save(codon_states,file = "codon_states")
save(codon_outputs,file = "codon_outputs")