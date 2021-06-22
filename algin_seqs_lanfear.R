begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
ncbi_seqs = readDNAStringSet(file = "sequences.fasta")
ref_seq = ncbi_seqs[1]
rm(ncbi_seqs)
seqs = readDNAStringSet(file = "global.fa")
seqs = short_names_seqs(seqs)
tree = read.tree("ft_SH.tree")
tree = short_names_tree(tree)
tree_seqs = return_only_tree_seqs(seqs,tree)
aligned_seqs = tree_seqs
save(aligned_seqs,file = "all_aligned_seqs")

#################
# relevant sites: 
setwd(paste0("/Users/keren/Dropbox/covid/new/vars/"))
load("site_details")
ref_gene_sites = regions$site[which(regions$ref_site!=0 &regions$gene!=0)]
setwd(paste0(begin_path,"vars/"))
library("vcfR")
bad_sites = read.vcfR("problematic_sites_sarsCov2.mask.vcf" , verbose = FALSE )
bad_sites = as.numeric(bad_sites@fix[252:502])
bad_sites = unique(c(bad_sites,1:68))
tmp = 1:29903
tmp = tmp[-bad_sites]
which(is.element(tmp,ref_gene_sites)==FALSE)
ref_seq_lanfear = strsplit(as.character(aligned_seqs[[1]]),split = "")[[1]]
ref_seq2 = strsplit(as.character(ref_seq[[1]]),split = "")[[1]]
ref_seq2 = ref_seq2[tmp]
for (i in 1:length(tmp)){
  if (ref_seq_lanfear[i]!=ref_seq2[i]){
    print(i)
    break
  }
}
site_to_alignment = cbind(1:length(ref_seq2),tmp)
colnames(site_to_alignment) = c("alignment_num","site")
site_to_alignment = data.frame(site_to_alignment)
setwd(paste0(begin_path,"vars/"))
save(site_to_alignment,file = "site_to_alignment")
#####################
# Find site patterns:
#####################
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
setwd(paste0(begin_path,"vars/site_patterns/",collapse = ""))
num_seqs = length(aligned_seqs)
seq_length = length(aligned_seqs[[1]])
for (j in 44){
  print(j)
  split_seqs = matrix("m",1000,length(aligned_seqs[[1]]))
  count = 0
  for (i in (j*1000+1):min(length(aligned_seqs),(j+1)*1000)){
    print(i)
    count = count+1
    split_seqs[count,] = strsplit(as.character(aligned_seqs[[i]]),split = "")[[1]]
  }
  if (j==44){
    split_seqs = split_seqs[1:(num_seqs-j*1000),]
  }
  save(split_seqs,file = paste0("split_seqs_",(j*1000+1),"_",(j+1)*1000,collapse=""))
}

##################
# Unite split_seqs: 
##################
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
setwd(paste0(begin_path,"vars/site_patterns/",collapse = ""))
file_names = dir()
################3
# Order file_names:
tmp = strsplit(file_names,split = "_")
file_nums = NULL
for (i in 1:length(file_names)){
  file_nums = c(file_nums,as.numeric(tmp[[i]][3]))
}
file_names = file_names[order(file_nums,decreasing = FALSE)]
##############
num_seqs = length(aligned_seqs)
seq_length = length(aligned_seqs[[1]])
rm(aligned_seqs)
count = 0
num_sites = 5000
print(file_names)
for (j in 1:6){
  plcs = c((1+(j-1)*num_sites):min(seq_length,(j*num_sites)))
  all_split_seqs = matrix("M",num_seqs,length(plcs))
  count = 0
  for (i in 1:45){#length(file_names)){
    print(i)
    setwd(paste0(begin_path,"/vars/site_patterns/"))
    load(file_names[i])
    print(count)
    print(dim(split_seqs))
    all_split_seqs[(count+1):(dim(split_seqs)[1]+count),] = split_seqs[,plcs]
    count = count + dim(split_seqs)[1]
    rm(split_seqs)
  }
  split_seqs = all_split_seqs
  setwd(paste0(begin_path,"/vars/site_patterns/all"))
  save(split_seqs,file = paste0("all_split_seqs_",plcs[length(plcs)],collapse =""))
}
###########
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
seq_length = length(aligned_seqs[[1]])
rm(aligned_seqs)
setwd(paste0(begin_path,"vars/site_patterns/all/",collapse = ""))
file_names = dir()
################3
# Order file_names:
tmp = strsplit(file_names,split = "_")
file_nums = NULL
for (i in 1:length(file_names)){
  file_nums = c(file_nums,as.numeric(tmp[[i]][4]))
}
file_names = file_names[order(file_nums,decreasing = FALSE)]
##############
count = 0
for (j in 1:6){
  setwd(paste0(begin_path,"vars/site_patterns/all",collapse = ""))
  load(file_names[j])
  site_patterns = site_patterns_plcs = no_subs_plcs = NULL
  for (i in 1:dim(split_seqs)[2]){
    print(i)
    curr_site_pattern = split_seqs[,i]
    if ( length(table(curr_site_pattern))==1) { 
      no_subs_plcs = c(no_subs_plcs,(i+count))
    }else{
      site_patterns = c(site_patterns,paste0(curr_site_pattern, collapse = ""))
      site_patterns_plcs = c(site_patterns_plcs,(i+count)) 
    }
  }
  count = count+dim(split_seqs)[2]
  file_name = paste0("site_patterns_",j,collapse = "")
  setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
  save(site_patterns,file = file_name)
  save(site_patterns_plcs,file = paste0("site_patterns_plcs_",j,collapse = ""))
  save(no_subs_plcs,file = paste0("no_subs_plcs_",j,collapse = ""))
}

######33
# Unite site_patterns:
setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
file_names = dir()
all_site_patterns = NULL
for (i in 7:12){
    print(i)
    load(file_names[i])
    all_site_patterns = c(all_site_patterns,site_patterns)
}
site_patterns = all_site_patterns
save(site_patterns,file = "all_site_patterns")

all_site_patterns_plcs = NULL
for (i in 13:18){
  print(i)
  load(file_names[i])
  all_site_patterns_plcs = c(all_site_patterns_plcs,site_patterns_plcs)
}
site_patterns_plcs = all_site_patterns_plcs
save(site_patterns_plcs,file = "all_site_patterns_plcs")

all_no_subs_plcs = NULL
for (i in 1:6){
  print(i)
  load(file_names[i])
  all_no_subs_plcs = c(all_no_subs_plcs,no_subs_plcs)
}
no_subs_plcs = all_no_subs_plcs
save(no_subs_plcs,file = "all_no_subs_plcs")

