get_problematic_muts = function(muts,bases,node_depth,site){
  parents = node_depth$parent
  mut_plcs = (which(bases!=bases[parents]))
  # problematic_muts = which(is.element(muts$parent,muts$base))
  # problematic_muts = c(problematic_muts,which(is.element(muts$base,muts$parent)))
  problematic_muts = NULL
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
  for (i in 1:length(mut_plcs)){
    if (children_bases[i,muts$before[i]]>0){
      problematic_muts = c(problematic_muts,i)
    }
  }
  problematic_muts = unique(problematic_muts)
  problematic_muts_mat = muts[problematic_muts,]
  return(problematic_muts)
}
get_muts = function(bases,node_depth,site,save_path){
  parents = node_depth$parent
  mut_plcs = (which(bases!=bases[parents]))
  if(length(mut_plcs)>0){
    muts = cbind(bases[parents][mut_plcs],bases[mut_plcs],parents[mut_plcs],names(bases)[mut_plcs])
    muts = data.frame(muts)
    colnames(muts) = c("before","after","parent","base")
    setwd(save_path)
    save(muts,file = paste0("muts.",site,".",dim(muts)[1],collapse = ""))
  }else{
    print("no_muts")
    print(site)
    muts = NULL
  }
  return(muts)
}
get_siblings_vals = function(muts,node_depth,bases){
  siblings = siblings_vals = siblings_num_children = matrix(0,dim(muts)[1],1)
  for (i in 1:dim(muts)[1]){
    siblings[i] = node_depth$children[which(node_depth$name==muts$parent[i])]
    siblings_vals[i] = paste0(bases[strsplit(siblings[i],split = "\\*")[[1]]],collapse = "*")
    siblings_num_children[i] = sum(node_depth$num_children[which(is.element(node_depth$name,strsplit(siblings[i],split = "\\*")[[1]]))])
  }
  muts = cbind(muts,siblings,siblings_vals,siblings_num_children)
  return(muts)
}



#debug muts at site 21169
# q_site = 21169
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("node_depth")
node_depth_Lanfear = node_depth
load("site_to_alignment")
load("site_details")
q_site = site_to_alignment$alignment_num[which(site_to_alignment$site==24034)]
setwd(paste0(begin_path,"vars/muts/"))
load("muts.21169.177")
alignment_site = matrix(q_site,dim(muts)[1],1)
ref_site = matrix(site_to_alignment$site[q_site],dim(muts)[1],1)
setwd("/Users/keren/Desktop/covid_files/new2/Lanfear trees/")
load(paste0("tree_list.",q_site))
problematic_muts = get_problematic_muts(muts,bases,node_depth,alignment_site)
problem = matrix(0,dim(muts)[1],1)
problem[problematic_muts,1] = 1
muts_Lanfear = cbind(muts,alignment_site,ref_site,problem)


muts_Lanfear = get_siblings_vals(muts = muts_Lanfear,node_depth = node_depth,bases = bases)
muts_Lanfear[which(muts_Lanfear$siblings_num_children==0),]
length(which(muts_Lanfear$siblings_num_children==0))
table(muts_Lanfear$problem)

begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("node_depth")
alignment_site = matrix(q_site,dim(muts)[1],1)
setwd(paste0(begin_path,"vars/muts/"))
# load("muts.21169.133")
load("muts.23891.20")
ref_site = matrix(site_to_alignment$site[q_site],dim(muts)[1],1)
setwd("/Users/keren/Desktop/covid_files/new2/ncbi trees/")
load(paste0("tree_list.",q_site))
problematic_muts = get_problematic_muts(muts,bases,node_depth,alignment_site)
problem = matrix(0,dim(muts)[1],1)
problem[problematic_muts,1] = 1
# muts_ncbi = get_muts(bases = bases,node_depth = node_depth,site = alignment_site[1],
#                     save_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/vars/muts/")
muts_ncbi = cbind(muts,alignment_site,ref_site,problem)

muts_ncbi = get_siblings_vals(muts = muts_ncbi,node_depth = node_depth,bases = bases)
muts_ncbi[which(muts_ncbi$siblings_num_children==0),]
length(which(muts_ncbi$siblings_num_children==0))
table(muts_ncbi$problem)

ancient_nodes.ncbi = c("node_47128","node_47129","node_47130","node_47131","node_75473","node_77902","node_77905","node_77906","node_77931")
rel_node = ancient_nodes.ncbi[which(is.element(ancient_nodes.ncbi,muts_ncbi$parent))]
muts_ncbi[which(muts_ncbi$parent==rel_node),]

setwd("/Users/keren/Dropbox/covid/new2/check_tree/compare_to_gisaid_full_tree/spikeprot0306.tar/spikeprot0306/")
spike_alignment = read.fasta("spikeprot0306.fasta")
spike_sites = c(21563:25384)
names_spike_aln = names(spike_alignment)
tmp = strsplit(names_spike_aln,split = "\\|")
q_a.acid = 677 #= (23593-21562)/3
a.acid_677 = a.acid_260 = short_names = matrix("0",length(spike_alignment),1)
for (i in 1:length(spike_alignment)){
  print(i)
  a.acid_677[i] = spike_alignment[[i]][677]
  a.acid_260[i] = spike_alignment[[i]][260]
  short_names[i] = tmp[[i]][4]
}
names(a.acid_677) = short_names
names(a.acid_260) = short_names
save(a.acid_260,file = "a.acid_260_spike")
save(a.acid_677,file = "a.acid_677_spike")

