# begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("all_aligned_seqs")
tree = read.tree("ft_SH.tree")
tree = short_names_tree(tree)
load("site_to_alignment")

ref_seq = "NC_045512"
# plcs = which(strsplit(as.character(aligned_seqs[[ref_seq]]),split = "")[[1]]!="-")
ref_regions = matrix(0,29903,8)
colnames(ref_regions) = c("ref_site","codon.pos","UTR","gene",
                          "mat_peptide","stem_loop","notes","notes2")
ref_regions = data.frame(ref_regions)
ref_regions$ref_site = 1:29903

ref_regions$UTR[1:265] = "5UTR"
ref_regions$gene[266:21555] = "ORF1ab"
ref_regions$codon.pos[13468:21555] = matrix(1:3, length(13468:21555),1)
ref_regions$codon.pos[266:13468] = matrix(1:3, length(266:13468),1)
ref_regions$mat_peptide[266:805] = "ORF1ab"
ref_regions$notes[266:805] = "nsp1"
ref_regions$mat_peptide[806:2719] = "ORF1ab"
ref_regions$notes[806:2719] = "nsp2"
ref_regions$notes[2720:8554] = "nsp3"
ref_regions$mat_peptide[2720:8554] = "ORF1ab"
ref_regions$notes[8555:10054] = "nsp4"
ref_regions$mat_peptide[8555:10054] = "ORF1ab"
ref_regions$notes[10055:10972] = "nsp5"
ref_regions$mat_peptide[10055:10972] = "ORF1ab"
ref_regions$notes[10973:11842] = "nsp6"
ref_regions$mat_peptide[10973:11842] = "ORF1ab"
ref_regions$notes[11843:12091] = "nsp7"
ref_regions$mat_peptide[11843:12091] = "ORF1ab"
ref_regions$notes[12092:12685] = "nsp8"
ref_regions$mat_peptide[12092:12685] = "ORF1ab"
ref_regions$notes[12686:13024] = "nsp9"
ref_regions$mat_peptide[12686:13024] = "ORF1ab"
ref_regions$notes[13025:13441] = "nsp10"
ref_regions$mat_peptide[13025:13441] = "ORF1ab"
ref_regions$notes[13442:16236] = "nsp12"
ref_regions$mat_peptide[13442:16236] = "ORF1ab"
ref_regions$notes[16237:18039] = "nsp13"
ref_regions$mat_peptide[16237:18039] = "ORF1ab"
ref_regions$notes[18040:19620] = "nsp14"
ref_regions$mat_peptide[18040:19620] = "ORF1ab"
ref_regions$notes[19621:20658] = "nsp15"
ref_regions$mat_peptide[19621:20658] = "ORF1ab"
ref_regions$notes[20659:21552] = "nsp16"
ref_regions$mat_peptide[20659:21552] = "ORF1ab"
ref_regions$codon.pos[266:13483] = matrix(1:3, length(266:13483),1)
ref_regions$notes2[13442:13480] = "nsp11"
ref_regions$mat_peptide[13442:13480] = "ORF1ab"
ref_regions$stem_loop[13476:13503] = "ORF1ab"
ref_regions$notes[13476:13503] = "stem_loop_1"
ref_regions$notes[13488:13542] = "stem_loop_2"
ref_regions$stem_loop[13488:13542] = "ORF1ab"
ref_regions$gene[21563:25384] = "S"
ref_regions$gene[25393:26220] = "ORF3a"
ref_regions$codon.pos[21563:25384] = matrix(1:3, length(21563:25384),1)
ref_regions$codon.pos[25393:26220] = matrix(1:3, length(25393:26220),1)
ref_regions$gene[26245:26472] = "E"
ref_regions$codon.pos[26245:26472] = matrix(1:3, length(26245:26472),1)
ref_regions$gene[26523:27191] = "M"
ref_regions$codon.pos[26523:27191] = matrix(1:3, length(26523:27191),1)
ref_regions$gene[27202:27387] = "ORF6"
ref_regions$codon.pos[27202:27387] = matrix(1:3, length(27202:27387),1)
ref_regions$gene[27394:27759] = "ORF7a"
ref_regions$codon.pos[27756:27887] = matrix(1:3, length(27756:27887),1)
ref_regions$codon.pos[27394:27759] = matrix(1:3, length(27394:27759),1)
ref_regions$gene[27756:27758] = "ORF7a_ORF7b"
ref_regions$gene[27756:27887] = "ORF7b"
ref_regions$gene[27894:28259] = "ORF8"
ref_regions$codon.pos[27894:28259] = matrix(1:3, length(27894:28259),1)
ref_regions$gene[28274:29533] = "N"
ref_regions$codon.pos[28274:29533] = matrix(1:3, length(28274:29533),1)
ref_regions$gene[29558:29674] = "ORF10"
ref_regions$codon.pos[29558:29674] = matrix(1:3, length(29558:29674),1)
ref_regions$stem_loop[29609:29644] = "ORF10"
ref_regions$notes[29609:29644] = "stem_loop_1"
ref_regions$stem_loop[29629:29657] = "ORF10"
ref_regions$notes2[29629:29657] = "stem_loop_2"
ref_regions$UTR[29675:29903] = "3UTR"
ref_regions$stem_loop[29728:29768] = "stem_loop_2_like"

regions = matrix(0,length(aligned_seqs[1][[1]]),9)
colnames(regions) = c("site",colnames(ref_regions))
regions = data.frame(regions)
regions[,2:9]=ref_regions[site_to_alignment$site,]
regions$site = 1:dim(regions)[1]

ref_seq = read.GenBank("NC_045512.2")
tmp = NULL
for (i in 1:length(ref_seq$NC_045512.2)){
  if (ref_seq$NC_045512.2[i]=="88"){
    letter = "A"
  }else if (ref_seq$NC_045512.2[i]=="18"){
    letter = "T"
  }else if(ref_seq$NC_045512.2[i]=="48"){
    letter = "G"
  }else if(ref_seq$NC_045512.2[i]=="28"){
    letter = "C"
  }
  tmp = c(tmp,letter)
}

ref_seq = matrix(0,dim(regions)[1],1)
ref_seq[which(regions$ref_site>0)] = tmp[site_to_alignment$site]
regions = cbind(regions,ref_seq)

codons = NULL
k=1
curr_codon_sites = which(regions$codon.pos>0)
curr_codon.pos = regions$codon.pos[curr_codon_sites]
seq = regions$ref_seq[unique(curr_codon_sites)]
while(k <=length(curr_codon_sites)){
  if (sum(curr_codon.pos[k:(k+2)]==c(1,2,3))==3){
    codons = c(codons,rep(paste0(seq[k:(k+2)],collapse = ""),3))
    k = k+3
  }else{
    codons = c(codons,"missing")
    print(paste("k=",k))
    k= k+1
  }
}

ref_codon = matrix(0,dim(regions)[1],1)
ref_codon[curr_codon_sites] = codons
regions = cbind(regions,ref_codon)

load("codons_to_amino_acids")
amino_acid = matrix(0,dim(regions)[1],1)
for (i in 1:dim(codons_to_amino_acids)[1]){
  print(i)
  curr_codon = rownames(codons_to_amino_acids)[i]
  amino_acid[which(regions$ref_codon==curr_codon)] = codons_to_amino_acids[i,2]
}

ref_amino_acid = amino_acid
regions = cbind(regions,ref_amino_acid)

setwd(paste0(begin_path,"vars/"))
save (regions,file = "site_details")
load("site_details")

########
# Add amino acid number (numbers start at 1 for each gene)
load("site_details")
gene_names = unique(regions$gene)
amino_acid_num = matrix("0",dim(regions)[1],1)
for (i in 1:length(gene_names)){
  if (gene_names[i]!="0"){
    print(i)
    gene_plcs = which(regions$ref_seq!="0" & regions$gene==gene_names[i])
    full_gene_plcs = c(regions$ref_site[gene_plcs[1]]:regions$ref_site[gene_plcs[length(gene_plcs)]])
    missing_plcs = which(is.element(full_gene_plcs,
                                    regions$ref_site[gene_plcs])==FALSE)
  }
}




load("site_details")
gene_names = unique(regions$gene)
amino_acid_num = matrix("0",dim(regions)[1],1)
for (i in 3:length(gene_names)){
  if (gene_names[i]!="0"){
    print(i)
    gene_plcs = which(regions$ref_seq!="0" & regions$gene==gene_names[i])
    count = 0
    j=1
    while (j<=(length(gene_plcs)-2)){
      if (count > 1439){
        bug = 1
        break
      }
      count = count+1
      # print(count)
      if (sum(regions$codon.pos[gene_plcs[j:(j+2)]]==c(1,2,3))==3){
        amino_acid_num[gene_plcs[j:(j+2)],1] = count    
        j = j+3
        print(j)
      }else{
        while(regions$codon.pos[gene_plcs[j]]!=1){
          amino_acid_num[gene_plcs[j]] = count
          j = j+1
          print(j)
        }
      }
    }
    while(j<=length(gene_plcs)){
      count = count+1
      amino_acid_num[gene_plcs[j]] = count
      j = j+1
    }
  }
}
regions = cbind(regions,amino_acid_num)
save(regions,file = "site_details")

plcs = which(regions$ref_seq=="0" & regions$codon.pos!="0")

regions[23620:23640,
        c("amino_acid_num","gene","codon.pos","ref_site","site")]
###########3 # Old version:
# regions$codon.pos[11118:11120] = c(1,2,3)
# regions$codon.pos[11123:11125] = c(3,1,2)
# regions$codon.pos[11132:11134] = c(3,1,2)
# regions$codon.pos[11137:11139] = c(2,3,1)
# 
# regions$gene[11118:11120] = "ORF1ab"
# regions$gene[11123:11125] = "ORF1ab"
# regions$gene[11132:11134] = "ORF1ab"
# regions$gene[11137:11139] = "ORF1ab"
# 
# regions$mat_peptide[11118:11120] = "ORF1ab"
# regions$mat_peptide[11123:11125] = "ORF1ab"
# regions$mat_peptide[11132:11134] = "ORF1ab"
# regions$mat_peptide[11137:11139] = "ORF1ab"
# regions$notes[11118:11120] = "nsp6"
# regions$notes[11123:11125] = "nsp6"
# regions$notes[11132:11134] = "nsp6"
# regions$notes[11137:11139] = "nsp6"
# 
# 
# regions$codon.pos[22411:22419] = c(2,3,1,2,3,1,2,3,1)
# regions$gene[22411:22419] = "S"
# 
# regions$codon.pos[25759:25761] = c(3,1,2)
# regions$gene[25759:25761] = "ORF3a"
# 
# regions$codon.pos[21441:21443] = c(2,3,1)
# regions$gene[21441:21443] = "ORF1ab"
# regions$mat_peptide[21441:21443] = "ORF1ab"
# regions$notes[21441:21443] = "nsp16"
# 
# regions$codon.pos[26156:26158] = c(1,2,3)
# regions$gene[26156:26158] = "ORF3a"
# 
# regions$codon.pos[28324:28326] = c(1,2,3)
# regions$gene[28324:28326] = "ORF8"
# 
# setwd("/Users/keren/Dropbox/covid/vars/")
# save (regions,file = "site_details")
# setwd("/Users/keren/Dropbox/EC2_files/")
# save (regions,file = "site_details")
# 
# 
# 
# #gaps_in_coding = c(11118:11120,11123:11125,11132:11134,11137:11139,
# # 21441:21443,22411:22419,25759:25761,26156:26158,
# # 28324:28326)
# # 29646:29647: This is also a gap of 2, an insertion appears only in one seq, 
# # relevant node is i=5350,  parent_name = "no_name_14339", child_name "MT520337.1"
# 
# 
# codon_pos_1_plcs = which(regions$codon.pos==1)
# codon_pos_2_plcs = which(regions$codon.pos==2)
# codon_pos_3_plcs = which(regions$codon.pos==3)
# codon_pos_plcs = cbind(codon_pos_1_plcs[4411:length(codon_pos_1_plcs)],
#                        codon_pos_2_plcs[4412:length(codon_pos_2_plcs)],
#                        codon_pos_3_plcs[4412:length(codon_pos_3_plcs)])
# colnames(codon_pos_plcs) = c("pos1","pos2","pos3")
# codon_pos_plcs = data.frame(codon_pos_plcs)
# codon_pos_plcs = cbind(codon_pos_plcs$pos1[4729:length(codon_pos_plcs$pos1)],
#                        codon_pos_plcs$pos2[4730:length(codon_pos_plcs$pos2)],
#                        codon_pos_plcs$pos3[4730:length(codon_pos_plcs$pos3)])
# colnames(codon_pos_plcs) = c("pos1","pos2","pos3")
# codon_pos_plcs = data.frame(codon_pos_plcs)
# 
# tmp = which(codon_pos_plcs$pos2-codon_pos_plcs$pos1!=1)
# tmp = which(codon_pos_3_plcs-codon_pos_1_plcs!=2)
# 
# codon_plcs = which(regions$codon.pos!=0)
# 
# # Sites 13468:13484 have double codon meaning - it should be 13468:21555 
# # but it overlaps with 266:13483 so 13468 is kept as codon.pos=3 but 
# # according to 13468:21555 it should be 1.
# # Sites 27756:27759 have the same problem (27394:27759,27756:27887)
# 
# ###########
# Check : #
load("ORF1ab1")
plcs = which(regions$gene=="ORF1ab")
ncbi_amino_acid = matrix(0,dim(regions),1)
j=1
count = 0
for (i in 1:length(plcs)){
  count = count+1
  ncbi_amino_acid[plcs[i]] = ORF1ab1[j]
  if (count==3){
    count = 0
    j = j+1
  }
}

check = cbind(ncbi_amino_acid,regions$ref_amino_acid)
colnames(check) = c("ncbi","regions")
plcs = which(check[,1]!=check[,2] & check[,1]!="0")

load("ORF1ab2")
plcs = which(regions$gene=="ORF1ab")
j=1
count = 0
for (i in 12244:length(plcs)){
  count = count+1
  ncbi_amino_acid[plcs[i]] = ORF1ab2[j]
  if (count==3){
    count = 0
    j = j+1
  }
}

ncbi_amino_acid = ncbi_amino_acid[-13524]
check = cbind(ncbi_amino_acid,regions$ref_amino_acid,regions$ref_amino_acid2)
colnames(check) = c("ncbi","regions","codon.pos2")
plcs = which(check[,1]!=check[,2] & check[,1]!="0")
plcs = which(check[,1]!=check[,3] & check[,1]!="0")

ref_codon.pos2 = ref_regions$codon.pos
ref_codon.pos2[13471:21555] = matrix(1:3,length(13471:21555))
codon.pos2 = matrix(0,dim(regions)[1],1)
codon.pos2[which(regions$ref_site!="0")] = ref_codon.pos2

load("site_details")
gene_names = unique(regions$gene)
amino_num = matrix(0,dim(regions)[1],1)
for (i in c(3:7,10:12)){
  print(gene_names[i])
  gene_start = min(regions$ref_site[which(regions$gene==gene_names[i])])
  gene_end = max(regions$ref_site[which(regions$gene==gene_names[i])])
  print((gene_end-gene_start+1)/3)
  gene_plcs = which(regions$gene==gene_names[i])
  amino_num[gene_plcs] = (((regions$ref_site[gene_plcs]-regions$codon.pos[gene_plcs])-(gene_start-1))/3)+1
  print(summary(amino_num[gene_plcs]))
  print(gene_start)
  print(gene_end)
}

# ORF1ab
gene_start = 13468
gene_end = 21555
print((gene_end-gene_start+1)/3)
gene_plcs = which(regions$ref_site>=gene_start & regions$ref_site<=gene_end)
amino_num[gene_plcs] = (((regions$ref_site[gene_plcs]-regions$codon.pos[gene_plcs])-(gene_start-1))/3)+1+4401
print(summary(amino_num[gene_plcs]))

gene_start = 266
gene_end = 13483
print((gene_end-gene_start+1)/3)
gene_plcs = which(regions$ref_site>=gene_start & regions$ref_site<=gene_end)
amino_num[gene_plcs] = (((regions$ref_site[gene_plcs]-regions$codon.pos[gene_plcs])-(gene_start-1))/3)+1
print(summary(amino_num[gene_plcs]))


# ORF7ab
gene_start = 27756
gene_end = 27887
print((gene_end-gene_start+1)/3)
gene_plcs = which(regions$ref_site>=gene_start & regions$ref_site<=gene_end)
amino_num[gene_plcs] = (((regions$ref_site[gene_plcs]-regions$codon.pos[gene_plcs])-(gene_start-1))/3)+1
print(summary(amino_num[gene_plcs]))

gene_start = 27394
gene_end = 27759
print((gene_end-gene_start+1)/3)
gene_plcs = which(regions$ref_site>=gene_start & regions$ref_site<=gene_end)
amino_num[gene_plcs] = (((regions$ref_site[gene_plcs]-regions$codon.pos[gene_plcs])-(gene_start-1))/3)+1
print(summary(amino_num[gene_plcs]))

regions$gene[which(regions$ref_site>=27756 & regions$ref_site<=27759)]="ORF7a"

regions = cbind(regions,amino_num)
save(regions,file = "site_details")
