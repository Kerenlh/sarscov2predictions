find_num_subs_per_site = function(data){
  data$ref_site = as.numeric(data$ref_site)
  ref_site = sort(unique(data$ref_site))
  num_subs = NULL
  for (i in 1:length(ref_site)){
    # print(i)
    num_subs = c(num_subs,sum(data[data$ref_site==ref_site[i], c("A","C","G","T")]))
  }
  num_subs = data.frame(num_subs)
  return(cbind(ref_site,num_subs))
}

print_top_muts = function(mat,mat_type){
  print(paste("Top ",mat_type))
  mat_per_site = sort(table(mat$ref_site),decreasing = T)
  print("ref_site: ")
  print(mat_per_site[1:30])
  print("alignment_site: ")
  print(sort(table(mat$alignment_site),decreasing = T)[1:30])
  return(mat_per_site)
}

load_files = function(begin_path){
  setwd(paste0(begin_path,"vars/"))
  load("codons_table")
  print(paste("dim(codons_table): ",dim(codons_table)))
  setwd(paste0(begin_path,"vars/debug_muts"))
  load("subs")
  load("insertions")
  load("deletions")
  # setwd(paste0(begin_path,"vars/tables_data"))
  # load("insertions_coding")
  # load("deletions_coding")
  print(paste("dim(insertions): ",dim(insertions)))
  print("problem table:")
  print(table(insertions$problem))
  print(paste("dim(deletions): ",dim(deletions)))
  print("problem table:")
  print(table(deletions$problem))
  print(paste("num subs: ",dim(subs)[1]))
  print("problem table:")
  print(table(subs$problem))
  print(paste("num subs: ",sum(codons_table$y)))
  print(paste("num transitions: ",sum(codons_table$transitions)))
  print(paste("num transversions: ",sum(codons_table$transversions)))
  print(paste("num synonymous subs: ",sum(codons_table$syn)))
  print(paste("num non-synonymous subs: ",sum(codons_table$non_syn)))
  print(paste("ti/tv: ",sum(codons_table$transitions)/sum(codons_table$transversions)))
  print(paste("syn/non_syn: ",sum(codons_table$syn)/sum(codons_table$non_syn)))
  
  # Stationary dist:
  data = codons_table
  subs_table = subs_table2 = matrix(0,5,5)
  colnames(subs_table) = colnames(subs_table2) = c("A","C","G","T","line")
  row.names(subs_table) = row.names(subs_table2) = c("A","C","G","T","-")
  for (i in 1:5){
    for (j in 1:5){
      subs_table2[i,j] = sum(as.numeric(data[which(data$base==row.names(subs_table2)[i]),
                                            colnames(subs_table2)[j]]))
      subs_table[i,j] = length(which(subs$before==row.names(subs_table)[i] & 
                                       subs$after==colnames(subs_table)[j]))
    }
  }
  Q = subs_table[1:4,1:4]
  Q2 = subs_table2[1:4,1:4]
  print("Subs matrix: ")
  print(Q)
  print("Subs matrix, codons_table: ")
  print(Q2)
  for (i in 1:4){
    Q[i,i] = -sum(Q[i,])
    Q2[i,i] = -sum(Q2[i,])
  }
  tmp = eigen(t(Q))
  tmp2 = eigen(t(Q2))
  print("Stationary dist.:")
  print(tmp$vectors[,4]/sum(tmp$vectors[,4]))
  print("Stationary dist., codons_table:")
  print(tmp2$vectors[,4]/sum(tmp2$vectors[,4]))
  
  # Top deletions:
  # top_deletions = sort(table(deletions$ref_start),decreasing = T)[1:30]
  # print(top_deletions)
  # top_deletions = sort(table(deletions$start),decreasing = T)[1:30]
  # deletion_sites = as.numeric(names(top_deletions))
  # print(deletion_sites)
  # table(deletions$ref_start)
  
  # Top subs:
  print("Top subs, codons_table:")
  subs_per_site2 = find_num_subs_per_site(codons_table)
  top_subs = subs_per_site2[order(subs_per_site2$num_subs,decreasing = T)[1:30],]
  print(top_subs)
  print("Top subs:")
  subs_per_site = print_top_muts(subs,"subs")
  deletions_per_site = print_top_muts(deletions,"deletions")
  insertions_per_site = print_top_muts(insertions,"insertions")
  
  return(list(subs,insertions,deletions,subs_per_site2))
}

tmp = load_files( "/Users/keren/Dropbox/covid/new2/")
gisaid.codons_table = tmp[1][[1]]; gisaid.insertions = tmp[2][[1]]; gisaid.deletions = tmp[3][[1]]
gisaid.subs_per_site = tmp[4][[1]]
tmp = load_files( "/Users/keren/Dropbox/covid/new2/ncbi_tree/")
ncbi.codons_table = tmp[1][[1]]; ncbi.insertions = tmp[2][[1]]; ncbi.deletions = tmp[3][[1]]
ncbi.subs_per_site = tmp[4][[1]]
tmp = load_files( "/Users/keren/Dropbox/covid/new/")
old.codons_table = tmp[1][[1]]; old.insertions = tmp[2][[1]]; old.deletions = tmp[3][[1]]
old.subs_per_site = tmp[4][[1]]

# Ref seq dist:
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("site_details")
setwd(paste0(begin_path,"vars/reconstructed_seqs"))
load("all_seqs")
ref_seq = strsplit(all_seqs["NC_045512"][[1]],split = "")[[1]]
print(table(ref_seq)/length(ref_seq))

# Check deletions:
setwd("/Users/keren/Dropbox/covid/new2/ncbi_tree/vars/muts/")
load("muts.7502.432")
length(which(muts$after == "-"))
load("problematic_muts.7502.47")
setwd("/Users/keren/Dropbox/covid/new2/vars/muts/")
load("muts.7502.267")
load("muts.222.345")
old.deletions[which(old.deletions$ref_start<=7594 & old.deletions$ref_end>=7594),]

gisaid.top_subs = gisaid.subs_per_site[order(gisaid.subs_per_site$num_subs,decreasing = T)[1:100],]
ncbi.top_subs = ncbi.subs_per_site[order(ncbi.subs_per_site$num_subs,decreasing = T)[1:100],]
length(intersect(gisaid.top_subs$ref_site,ncbi.top_subs$ref_site))

# Find sites in Spike
setwd(paste0(begin_path,"vars/"))
load("site_details")
summary(regions$site[which(regions$gene=="S")])
gisaid.subs_per_site[which(as.numeric(names(gisaid.subs_per_site))>=21426 & 
        as.numeric(names(gisaid.subs_per_site))<=25234)][1:10]
spike_sites = c(21563:25384)
gisaid.codons_table[which(gisaid.codons_table$ref_site==23593),]
q_site = 23593
q_a.acid = (23593-21562)/3
ncbi.codons_table[which(ncbi.codons_table$ref_site==23593),]

gisaid.deletions[which(gisaid.deletions$ref_site==22342),]
q_site = 22342
q_a.acid = (22342-21562)/3
ncbi.deletions[which(ncbi.deletions$ref_site==22342),]
