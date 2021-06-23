# begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
ref_seq = aligned_seqs[which(names(aligned_seqs)=="NC_045512")]
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
seqs = readDNAStringSet(file = "global_1_4_21.fa")
seqs = short_names_seqs(seqs)
old_seqs_names = names(aligned_seqs)
rm(aligned_seqs)
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
new_seqs_plcs = which(is.element(names(seqs),old_seqs_names)==FALSE)
new_seqs = seqs[new_seqs_plcs]
save(new_seqs,file = "new_seqs_for_prediction")
rm(seqs)
#################
# relevant sites: 
setwd(paste0("/Users/keren/Dropbox/covid/new/vars/"))
load("site_details")
ref_gene_sites = regions$site[which(regions$ref_site!=0 &regions$gene!=0)]
setwd(paste0(begin_path,"vars/"))
library("vcfR")
bad_sites = read.vcfR("prediction_problematic_sites_sarsCov2.mask.vcf" , verbose = FALSE )
bad_sites = as.numeric(bad_sites@fix[252:502])
bad_sites = unique(c(bad_sites,1:68))
tmp = 1:29903
tmp = tmp[-bad_sites]
which(is.element(tmp,ref_gene_sites)==FALSE)
ref_seq_lanfear = strsplit(as.character(seqs[which(names(seqs)=="NC_045512")]),split = "")[[1]]
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
save(site_to_alignment,file = "site_to_alignment_prediction")
#####################
# Find site patterns:
#####################
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("new_seqs_for_prediction")
setwd(paste0(begin_path,"vars/site_patterns_prediction/",collapse = ""))
aligned_seqs = new_seqs
num_seqs = length(aligned_seqs)
seq_length = length(aligned_seqs[[1]])
for (j in 0:39){
  print(j)
  split_seqs = matrix("m",1000,length(aligned_seqs[[1]]))
  count = 0
  for (i in (j*1000+1):min(length(aligned_seqs),(j+1)*1000)){
    print(i)
    count = count+1
    split_seqs[count,] = strsplit(as.character(aligned_seqs[[i]]),split = "")[[1]]
  }
  if (j==39){
    split_seqs = split_seqs[1:(num_seqs-j*1000),]
  }
  save(split_seqs,file = paste0("split_seqs_",(j*1000+1),"_",(j+1)*1000,collapse=""))
}

##################
# Unite split_seqs: 
##################
setwd(paste0(begin_path,"vars/"))
load("new_seqs_for_prediction")
setwd(paste0(begin_path,"vars/site_patterns_prediction/",collapse = ""))
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
aligned_seqs = new_seqs
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
  for (i in 1:length(file_names)){
    print(i)
    setwd(paste0(begin_path,"/vars/site_patterns_prediction/"))
    load(file_names[i])
    print(count)
    print(dim(split_seqs))
    all_split_seqs[(count+1):(dim(split_seqs)[1]+count),] = split_seqs[,plcs]
    count = count + dim(split_seqs)[1]
    rm(split_seqs)
  }
  split_seqs = all_split_seqs
  setwd(paste0(begin_path,"/vars/site_patterns_prediction/all"))
  save(split_seqs,file = paste0("all_split_seqs_",plcs[length(plcs)],collapse =""))
}
###############

setwd(paste0(begin_path,"vars/"))
load("new_seqs_for_prediction")
aligned_seqs = new_seqs
seq_length = length(aligned_seqs[[1]])
rm(aligned_seqs)
setwd(paste0(begin_path,"vars/site_patterns_prediction/all/",collapse = ""))
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
  setwd(paste0(begin_path,"vars/site_patterns_prediction/all",collapse = ""))
  load(file_names[j])
  site_patterns = site_patterns_plcs = no_subs_plcs = NULL
  for (i in 1:dim(split_seqs)[2]){
    print(i)
    curr_site_pattern = split_seqs[,i]
    if ( length(table(curr_site_pattern))==1) { 
      no_subs_plcs = c(no_subs_plcs,(i+count))
    }else if (length(table(curr_site_pattern))==2 & is.element("-",names(table(curr_site_pattern))) ){
      no_subs_plcs = c(no_subs_plcs,(i+count))
    }
    else{
      site_patterns = c(site_patterns,paste0(curr_site_pattern, collapse = ""))
      site_patterns_plcs = c(site_patterns_plcs,(i+count)) 
    }
  }
  count = count+dim(split_seqs)[2]
  file_name = paste0("site_patterns_",j,collapse = "")
  setwd(paste0(begin_path,"vars/site_patterns_prediction/all2",collapse = ""))
  save(site_patterns,file = file_name)
  save(site_patterns_plcs,file = paste0("site_patterns_plcs_",j,collapse = ""))
  save(no_subs_plcs,file = paste0("no_subs_plcs_",j,collapse = ""))
}

######33
# Unite site_patterns:
setwd(paste0(begin_path,"vars/site_patterns_prediction/all2",collapse = ""))
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

###########
# Find sites with no subs in the tree dataset

begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("new_seqs_for_prediction")
load("codons_table")
subs_ref_sites = unique(codons_table$ref_site[which(codons_table$y>0)])
load("site_details")
codon_ref_sites = unique(regions$ref_site[which(regions$codon.pos>0)])
no_subs_ref_sites = codon_ref_sites[which(is.element(codon_ref_sites,subs_ref_sites)==FALSE)]

save(no_subs_ref_sites,file = "no_subs_ref_sites")

load("site_to_alignment_prediction")
setwd(paste0(begin_path,"vars/site_patterns_prediction/all2",collapse = ""))
load("all_no_subs_plcs")
load("all_site_patterns_plcs")
load("all_site_patterns")
no_subs_plcs_prediction_ref_sites = site_to_alignment$site[no_subs_plcs]

rel_sites_no_subs = no_subs_ref_sites[which(is.element(no_subs_ref_sites,no_subs_plcs_prediction_ref_sites))]
rel_sites_subs = no_subs_ref_sites[which(is.element(no_subs_ref_sites,no_subs_plcs_prediction_ref_sites)==FALSE)]

subs = problematic_site_patterns =  NULL
for (i in 1:length(rel_sites_subs)){
  # print(i)
  curr_alignment_site = site_to_alignment$alignment_num[which(site_to_alignment$site==rel_sites_subs[i])]
  split_site_pattern = strsplit(site_patterns[which(site_patterns_plcs==curr_alignment_site)],split = "")[[1]]
  table_site_pattern = table(split_site_pattern)
  table_site_pattern = table_site_pattern[which(is.element(names(table_site_pattern),c("A","C","G","T")))]
  if (length(table_site_pattern)==2 & min(table_site_pattern)==1){
    before = names(table_site_pattern)[which.max(table_site_pattern)]
    after = names(table_site_pattern)[which.min(table_site_pattern)]
    sub_seq = strsplit(as.character(new_seqs[which(split_site_pattern==after)]),split = "")[[1]]
    if (sub_seq[curr_alignment_site]!=after){
      print("bug!")
    }
    subs = rbind(subs,c(before,sub_seq[(curr_alignment_site-2):(curr_alignment_site+2)],
                        rel_sites_subs[i],curr_alignment_site))
  }else{
    problematic_site_patterns = c(problematic_site_patterns,rel_sites_subs[i])
    print("problematic site pattern:")
    print(rel_sites_subs[i])
    print(table_site_pattern)
  }
}

colnames(subs) = c("before","L2","L1","after","R1","R2","ref_site","aln_site")
subs = data.frame(subs)

setwd(paste0(begin_path,"vars/prediction/"))
new_subs_prediction = subs
save(new_subs_prediction,file = "new_subs_prediction")

rel_ref_sites_no_subs = rel_sites_no_subs
save(rel_ref_sites_no_subs,file = "rel_ref_sites_no_subs")







