# begin_path = "/Users/keren/Dropbox/covid/new2/"
# begin_path = "C:/Users/Keren/Dropbox/covid/new2/"
begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("node_depth")
load("nodes_names")
load("site_details")
setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
load("all_seqs")
leaves_names = node_depth$name[which(node_depth$num_children==0)]

codon_plcs = which(regions$codon.pos!=0)

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

for (p in 61:89){#89){
  parallel_nums = p
  cluster_delta = 500
  codon_states = "tmp"
  init_outputs = matrix(0,dim(regions)[1],10)
  colnames(init_outputs) = c("A","C","G","T","line","q","s","M","exposure","zero_branch")
  init_outputs = data.frame(init_outputs)
  codon_outputs = init_outputs[1,]
  edges_with_no_subs = NULL
  start_time = Sys.time()
  for (i in c((1+(parallel_nums-1)*cluster_delta):
              min((parallel_nums*cluster_delta),length(leaves_names)))){
    print(i)
    curr_name = leaves_names[i]
    edge_length = 1
    curr_seq = strsplit(all_seqs[curr_name][[1]],split = "")[[1]]
    curr_outputs = init_outputs
    curr_outputs$exposure = edge_length

    gaps_plcs = which(curr_seq=="-")
    short_parent_sites = c(1:length(curr_seq))
    if (length(gaps_plcs)>0){
      short_parent = curr_seq[-gaps_plcs]
      short_parent_sites = short_parent_sites[-gaps_plcs]
      tmp_check = which(is.element(codon_plcs,gaps_plcs))
      if (length(tmp_check)>0){
        curr_codon_plcs = codon_plcs[-which(is.element(codon_plcs,gaps_plcs))]
      }else{
        curr_codon_plcs = codon_plcs
      }
      curr_sites = regions$site[-gaps_plcs]
    }else{
      short_parent = curr_seq
      curr_codon_plcs = codon_plcs
      curr_sites = regions$site
    }
    short_parent = cbind(short_parent_sites,short_parent)
    colnames(short_parent) = c("sites","seq")
    short_parent = data.frame(short_parent)
    short_parent$sites = as.numeric(short_parent$sites)
    codon_short_parent = curr_seq[curr_codon_plcs]
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
                               curr_seq[curr_codon_plcs],codons,
                               regions$site[curr_codon_plcs])
    codon_uncertainty_plcs = which(curr_seq[curr_codon_plcs]=="?" | 
                                     curr_outputs$q[curr_codon_plcs]=="1")
    tmp = update_states_and_outputs(curr_states = curr_codon_states,
                                    states = codon_states,
                                    curr_outputs = curr_outputs[curr_codon_plcs,],
                                    outputs = codon_outputs,
                                    uncertainty_plcs = codon_uncertainty_plcs)
    codon_states = tmp$states
    codon_outputs = tmp$outputs
  }
  print(Sys.time()-start_time)
  codon_states = codon_states[-1]
  codon_outputs = codon_outputs[-1,]
  
  setwd(paste0(begin_path,"vars/tables_data_prediction"))
  save(codon_states,file = paste0("leaves_codon_states_",parallel_nums,collapse =""))
  save(codon_outputs,file = paste0("leaves_codon_outputs_",parallel_nums,collapse =""))
}

begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/tables_data_prediction/"))
file_names = dir()
load("leaves_codon_outputs_1")
load("leaves_codon_states_1")
all_codon_states = codon_states
all_codon_outputs = codon_outputs
for (i in c(2:89)){ #77
  load(paste0("leaves_codon_states_",i))
  load(paste0("leaves_codon_outputs_",i))
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
