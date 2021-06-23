begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
load("all_seqs")
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
tree = read.tree("ft_SH.tree")
load("all_aligned_seqs")
load("codons_table")
# load("base_table")
load("site_details")



setwd("/Users/keren/Dropbox/covid/models/")
load("all_models")
D614G = codons_table[which(codons_table$site>=23469 & codons_table$site<=23471),]
501*3
regions[plcs[1501:1503],]
N501Y = codons_table[which(codons_table$site>=23130 & codons_table$site<=23132),]
del_69_70 = codons_table[which(codons_table$site>=21825 & codons_table$site<=21830),]
del_144_145 = codons_table[which(codons_table$site>=21992 & codons_table$site<=21997),]

get_amino_muts = function(mut_amino_plc){
  S_plcs = which(regions$gene=="S")
  mut_plcs = regions$site[plcs[((mut_amino_plc*3)-2):(mut_amino_plc*3)]]
  return(codons_table[which(is.element(codons_table$site,mut_plcs)),])
}
get_amino_muts(501)
amino_muts_plcs = c(69,70,144,145,501,570,614,681,716,982,1118)
data = NULL
for (i in 1:length(amino_muts_plcs)){
  curr_data = get_amino_muts(amino_muts_plcs[i])
  data = rbind(data,cbind(curr_data,
                          matrix(amino_muts_plcs[i],dim(curr_data)[1],1)))
}
colnames(data)[dim(data)[2]] = "amino_acid_plc"
write.csv(data,file = "VUI_202012_01.csv")
find_num_subs_per_site = function(data){
  data$site = as.numeric(data$site)
  # data[,c("site","A","C","G","T")] = as.numeric(data[,c("site","A","C","G","T")])
  all_sites = sort(unique(data$site))
  num_subs = NULL
  for (i in 1:length(all_sites)){
    print(i)
    num_subs = c(num_subs,sum(data[data$site==all_sites[i], c("A","C","G","T")]))
  }
  return(cbind(all_sites,num_subs))
}


base_table2 = base_table[-which(base_table$base=="-" | base_table$line>0),]
codons_table2 = codons_table[-which(codons_table$base=="-" | codons_table$line>0),]
num_subs_base = find_num_subs_per_site(base_table2)
order(num_subs_base,decreasing = TRUE)
num_subs_base[]
sort(num_subs_base,decreasing = TRUE)

num_subs_codons = find_num_subs_per_site(codons_table)
num_subs_codons = data.frame(num_subs_codons)
plot(num_subs_codons$all_sites,num_subs_codons$num_subs)


data = base_table
data = base_table[which(base_table$gene==0),]
data = codons_table
subs_table = matrix(0,5,5)
colnames(subs_table) = c("A","C","G","T","line")
row.names(subs_table) = c("A","C","G","T","-")
for (i in 1:5){
  for (j in 1:5){
    subs_table[i,j] = sum(as.numeric(data[which(data$base==row.names(subs_table)[i]),
                                          colnames(subs_table)[j]]))
  }
}
Q = subs_table[1:4,1:4]
for (i in 1:4){
  Q[i,i] = -sum(Q[i,])
}
for(i in 1:4){
  Q[i,] = Q[i,]/avg_ACGT[i]
}

tmp = eigen(t(Q))
tmp$vectors[,4]/sum(tmp$vectors[,4])

Q = subs_table[1:5,1:5]
for (i in 1:5){
  Q[i,i] = -sum(Q[i,])
}
tmp = eigen(t(Q))
tmp$vectors[,5]/sum(tmp$vectors[,5])

######
# Transition/transversion:
base_table$transitions = as.numeric(base_table$transitions)
base_table$transversions = as.numeric(base_table$transversions)
sum(base_table$transitions)/sum(base_table$transversions)
non_coding_plcs = which(base_table$gene==0)
sum(base_table$transitions[non_coding_plcs])/
  sum(base_table$transversions[non_coding_plcs])
######

seq = strsplit(all_tree_list[[root_name]]@probs,split = "")[[1]]
table(seq)/length(seq)

load("subs_regions")
subs_regions = data.frame(subs_regions)
data = subs_regions
subs_table = matrix(0,4,4)
colnames(subs_table) = c("A","C","G","T")
row.names(subs_table) = c("A","C","G","T")
for (i in 1:4){
  for (j in 1:4){
    subs_table[i,j] = length(which(data$before==row.names(subs_table)[i] &
                                     data$after==colnames(subs_table)[j]))
  }
}


find_node_times = function(node_times, curr_node_name,parent_time){
  curr_edge = which(tree$edge[,1]==all_tree_list[[ curr_node_name]]@parent_edge &
                      tree$edge[,2]==all_tree_list[[ curr_node_name]]@edge)
  print(curr_node_name)
  if ( curr_node_name=="no_name_12225"){
    edge_time = 0
  }else{
    edge_time = tree$edge.length[curr_edge]
  }
  curr_time = parent_time + edge_time
  node_times[curr_edge,] = c(curr_time,edge_time)
  if (all_tree_list[[ curr_node_name]]@leaf==FALSE){
    for (i in 1:length(all_tree_list[[ curr_node_name]]@children_names)){
      node_times = find_node_times(node_times = node_times,
                                   curr_node_name = all_tree_list[[ curr_node_name]]@children_names[i],
                                   parent_time = curr_time)
    }
  }
  return(node_times)
}

find_ACGT_table = function(ACGT_table, curr_node_name){
  print(curr_node_name)
  tmp = table(strsplit(all_tree_list[[curr_node_name]]@probs,split = "")[[1]])
  ACGT_table[all_tree_list[[curr_node_name]]@edge,names(tmp)] = tmp
  if (all_tree_list[[ curr_node_name]]@leaf==FALSE){
    for (i in 1:length(all_tree_list[[ curr_node_name]]@children_names)){
      ACGT_table = find_ACGT_table(ACGT_table,
                                   curr_node_name = all_tree_list[[ curr_node_name]]@children_names[i])
    }
  }
  return(ACGT_table)
}

# root_edge = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
root_edge = 12225
root_name = paste0(c("no_name_",root_edge),collapse = "")
node_times = matrix(0,dim(tree$edge)[1]+1,2)
node_times = find_node_times(node_times = node_times, 
                             curr_node_name = root_name, parent_time = 0)
colnames(node_times) = c("time","edge_length")

ACGT_table = matrix(0,length(node_times),6)
colnames(ACGT_table) = c("-","?","A","C","G","T")
ACGT_table = find_ACGT_table(ACGT_table = ACGT_table, 
                             curr_node_name = root_name)
tmp = (t(node_times[,"edge_length"])%*%ACGT_table)
avg_ACGT = tmp[3:6]/sum(tmp[3:6])

mean_ACGT = matrix(0,round(max(node_times)),dim(ACGT_table)[2])
colnames(mean_ACGT) = colnames(ACGT_table)
weights = matrix(0,round(max(node_times)),1)
for (i in 1:round(max(node_times))){
  plcs = which(node_times>(i-1) & node_times<i)
  if (length(plcs)>1){
    mean_ACGT[i,] = apply(ACGT_table[plcs,],2,mean)
  }
  else if(length(plcs)==1){
    mean_ACGT[i,] = ACGT_table[plcs,]
  }
  weights[i] = length(plcs)
}
mean_ACGT = data.frame(mean_ACGT)
plcs = which(mean_ACGT$T!=0)
mean_ACGT = mean_ACGT[plcs,]
times = c(1:round(max(node_times)))[plcs]
weights = weights[plcs]

plot(times,mean_ACGT$A,cex = weights/2000,ylim = c(5400,9650))
par(new=TRUE)
plot(times,mean_ACGT$C,cex = weights/2000,col="blue",ylim = c(5400,9650))
par(new=TRUE)
plot(times,mean_ACGT$G,cex = weights/2000,col="green",ylim = c(5400,9650))
par(new=TRUE)
plot(times,mean_ACGT$T,cex = weights/2000,col="red",ylim = c(5400,9650))


ACGT_table2 = ACGT_table[-which(node_times>100),]
node_times2 = node_times[-which(node_times>100)]
plot(node_times,ACGT_table[,"T"])
abline(lm(ACGT_table[,"T"] ~ node_times))
plot(times,mean_ACGT$T,cex = weights/1000)
plot(times,mean_ACGT$T)



scatter(times,mean_ACGT$T,circles = weights)

tmp2 = cbind(1:round(max(tmp$time)),mean_mean_ACGT)
colnames(tmp2) =c("time","\\-","\\?","A","C","G","T")
tmp2 = tmp2[order(tmp2[,"time"]),]
tmp2 = data.frame(tmp2)
tmp2 = tmp2[which(tmp2$T!=0),]



unique_node_times = unique(node_times)
mean_ACGT = matrix(0,length(unique_node_times),dim(ACGT_table)[2])
colnames(mean_ACGT) = colnames(ACGT_table)
for (i in 1:length(unique_node_times)){
  plcs = which(node_times==unique_node_times[i])
  if (length(plcs)>1){
    mean_ACGT[i,] = apply(ACGT_table[plcs,],2,mean)
  }
  else{
    mean_ACGT[i,] = ACGT_table[plcs,]
  }
}

tmp = cbind(unique_node_times,mean_ACGT)
colnames(tmp)[1] = "time"
tmp = tmp[order(tmp[,"time"]),]
tmp = data.frame(tmp)

mean_mean_ACGT = matrix(0,round(max(tmp$time)),dim(mean_ACGT)[2])
for (i in 1:round(max(tmp$time))){
  plcs = which(tmp$time>(i-1) & tmp$time<i)
  if (length(plcs)>1){
    mean_mean_ACGT[i,] = apply(mean_ACGT[plcs,],2,mean)
  }
  else if(length(plcs)==1){
    mean_mean_ACGT[i,] = mean_ACGT[plcs,]
  }
}
tmp2 = cbind(1:round(max(tmp$time)),mean_mean_ACGT)
colnames(tmp2) =c("time","\\-","\\?","A","C","G","T")
tmp2 = tmp2[order(tmp2[,"time"]),]
tmp2 = data.frame(tmp2)
tmp2 = tmp2[which(tmp2$T!=0),]
##########3
# new stat dist. :
find_ACGT_table = function(ACGT_table, curr_node_name){
  print(curr_node_name)
  tmp = table(strsplit(all_tree_list[[curr_node_name]]@probs,split = "")[[1]])
  ACGT_table[all_tree_list[[curr_node_name]]@edge,names(tmp)] = tmp
  if (all_tree_list[[ curr_node_name]]@leaf==FALSE){
    for (i in 1:length(all_tree_list[[ curr_node_name]]@children_names)){
      ACGT_table = find_ACGT_table(ACGT_table,
                                   curr_node_name = all_tree_list[[ curr_node_name]]@children_names[i])
    }
  }
  return(ACGT_table)
}




#########33
# bat and other seqs:
tmp = read.GenBank("MG772934.1")
tmp = read.GenBank("MN611519.1")
tmp = read.GenBank("MN611520.1")
tmp = read.GenBank("MK204393.1")
tmp = read.GenBank("LR757998.1")



# Substitutions in stop codons:
plcs = which(codons_table$amino_acid=="Stop" & codons_table$non_syn>0)
sum(codons_table$non_syn[plcs])
new_stop_codons_plcs = which(codons_table$T>0 & codons_table$T.amino=="Stop")
new_stop_codons_plcs = which(codons_table$C>0 & codons_table$C.amino=="Stop")
new_stop_codons_plcs = which(codons_table$G>0 & codons_table$G.amino=="Stop")
new_stop_codons_plcs = which(codons_table$A>0 & codons_table$A.amino=="Stop")
