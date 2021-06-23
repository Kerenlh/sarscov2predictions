############################
# Synonymous substitutions #
############################
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("codons_table")
load("syn_ti_const2")
load("syn_tv_const2")
load("non_syn_ti_const2")
load("non_syn_tv_const2")
# codons_table$transitions = as.numeric(codons_table$transitions)
# codons_table$transversions = as.numeric(codons_table$transversions)
# transversion_to_transitions =
#   0.5*sum(codons_table$transversions)/sum(codons_table$transitions)
# save(transversion_to_transitions,file = "transversion_to_transitions")
# ratio = transversion_to_transitions

edges = matrix("black",nrow = 64, ncol = 64)
colnames(edges) <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG",
                     "ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC",
                     "CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA",
                     "GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT",
                     "TTA","TTC","TTG","TTT")
rownames(edges) = colnames(edges)
codon_names = colnames(edges)
groups = edges
g1 = c("GCT", "GCC", "GCA", "GCG")
g2 = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")
g3 = c("AAT", "AAC")
g4 = c("GAT", "GAC")
g5 = c("TGT", "TGC")
g6 = c("CAA", "CAG")
g7 = c("GAA", "GAG")
g8 = c("GGT", "GGC", "GGA", "GGG")
g9 = c("CAT", "CAC")
g10 = c("ATT", "ATC", "ATA")
g11 = c("ATG")
g12 = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")
g13 = c("AAA", "AAG")
g14 = c("TTT", "TTC")
g15 = c("CCT", "CCC", "CCA", "CCG")
g16 = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")
g17 = c("ACT", "ACC", "ACA", "ACG")
g18 = c("TGG")
g19 = c("TAT", "TAC")
g20 = c("GTT", "GTC", "GTA", "GTG")
g21 = c("TAA", "TGA", "TAG")

# .1/2/3 number of possible mutations in a syn group in the i'th codon position
get_group_sizes = function(ratio,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21){
  group_sizes.3 = matrix(0,64,1)
  row.names(group_sizes.3) = codon_names
  group_sizes.2 = group_sizes.1 = group_sizes.3
  group_sizes.3[c(g1,g8,g15,g16,g17,g20),1] = 1+2*ratio
  group_sizes.3[c(g2,g12),1] = 1+2*ratio
  group_sizes.3[c("AGA", "AGG","TAA","TAG"),1] = 1
  group_sizes.3[c(g3,g4,g5,g6,g7,g9,g13,g14,g19),1] = 1
  group_sizes.3[c(g10),1] = 1+ratio
  group_sizes.3["ATA",1] = ratio
  group_sizes.3[c("TTA", "TTG","AGT", "AGC"),1] = 1
  
  group_sizes.2[c("TAA", "TGA"),1] = 1
  group_sizes.1[c("CGA", "CGG", "AGA", "AGG"),1] = ratio
  group_sizes.1[c("TTA", "TTG", "CTA", "CTG"),1] = 1
  return(list(group_sizes.1,group_sizes.2,group_sizes.3))
}

tmp = get_group_sizes(ratio=syn_tv_const,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21)
syn_group_sizes.1 = tmp[[1]];syn_group_sizes.2 = tmp[[2]];syn_group_sizes.3 = tmp[[3]]
tmp = get_group_sizes(ratio=non_syn_tv_const,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21)
non_syn_group_sizes.1 = tmp[[1]];non_syn_group_sizes.2 = tmp[[2]];non_syn_group_sizes.3 = tmp[[3]]

group_names = c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Start/Met",
                "Leu","Lys","Phe","Pro","Ser","Thr","Trp","Tyr","Val","Stop")
group_sizes = c(4,6,2,2,2,2,2,4,2,3,1,6,2,2,4,6,4,1,2,4,3)

group_names2 = c("A","R","N","D","C","Q","E","G","H","I","M",
                 "L","K","F","P","S","T","W","Y","V","Stop")
codons_to_amino_acids = matrix(0,64,2)
row.names(codons_to_amino_acids) = codon_names
for (j in 1:length(group_names)){
  gr = paste("g",j,sep="")
  curr.g = eval(parse(text = gr))
  codons_to_amino_acids[curr.g,1] = group_names[j]
  codons_to_amino_acids[curr.g,2] = group_names2[j]
}

save(groups,file = "groups")
save(group_sizes,file = "group_sizes")
save(codons_to_amino_acids,file = "codons_to_amino_acids")
save(syn_group_sizes.1,file = "syn_group_sizes.1")
save(syn_group_sizes.2,file = "syn_group_sizes.2")
save(syn_group_sizes.3,file = "syn_group_sizes.3")
save(non_syn_group_sizes.1,file = "non_syn_group_sizes.1")
save(non_syn_group_sizes.2,file = "non_syn_group_sizes.2")
save(non_syn_group_sizes.3,file = "non_syn_group_sizes.3")


##############
# Functions: #
##############
add_CG = function(data){
  CG = matrix(0,dim(data)[1],1)
  CG[which(data$base=="C" & data$R.neighbor=="G")] = 1
  CG[which(data$base=="G" & data$L.neighbor=="C")] = 2
  data = cbind(data,CG)
  return(data)
}

add_transitions_transversions = function(data){
  transitions = transversions = matrix(0,dim(data),1)
  transitions[which(data$base=="A")] = data$G[which(data$base=="A")]
  transitions[which(data$base=="G")] = data$A[which(data$base=="G")]
  transitions[which(data$base=="C")] = data$T[which(data$base=="C")]
  transitions[which(data$base=="T")] = data$C[which(data$base=="T")]
  
  transversions[which(data$base=="A")] = data$C[which(data$base=="A")]+
    data$T[which(data$base=="A")]
  transversions[which(data$base=="C")] = data$A[which(data$base=="C")]+
    data$G[which(data$base=="C")]
  transversions[which(data$base=="G")] = data$C[which(data$base=="G")]+
    data$T[which(data$base=="G")]
  transversions[which(data$base=="T")] = data$A[which(data$base=="T")]+
    data$G[which(data$base=="T")]
  
  return(cbind(data,transitions,transversions))
}


add_amino_acid_exposure = function(data,codons_to_amino_acids,ratio,group_sizes.1,
                                   group_sizes.2,group_sizes.3){
  amino_acid = syn_exposure =  non_syn_exposure = matrix(0,dim(data)[1],1)
  # plcs = which(is.element(data$codon,rownames(codons_to_amino_acids)))
  # no_valid_codon_plcs = which(is.element(data$codon,rownames(codons_to_amino_acids)==FALSE))
  for (i in 1:dim(codons_to_amino_acids)[1]){
    print(i)
    curr_codon = rownames(codons_to_amino_acids)[i]
    amino_acid[which(data$codon==curr_codon)] = codons_to_amino_acids[i,2]
    syn_exposure[which(data$codon==curr_codon & data$codon.pos==1)] = 
      syn_group_sizes.1[curr_codon,1]
    non_syn_exposure[which(data$codon==curr_codon & data$codon.pos==1)] = 
      non_syn_group_sizes.1[curr_codon,1]
    
    syn_exposure[which(data$codon==curr_codon & data$codon.pos==2)] = 
      syn_group_sizes.2[curr_codon,1]
    non_syn_exposure[which(data$codon==curr_codon & data$codon.pos==2)] = 
      non_syn_group_sizes.2[curr_codon,1]
    
    syn_exposure[which(data$codon==curr_codon & data$codon.pos==3)] = 
      syn_group_sizes.3[curr_codon,1]
    non_syn_exposure[which(data$codon==curr_codon & data$codon.pos==3)] = 
      non_syn_group_sizes.3[curr_codon,1]
  }
  non_syn_exposure = (1+2*non_syn_tv_const)-non_syn_exposure
  return(cbind(data,amino_acid,syn_exposure,non_syn_exposure))
}

add_syn_non_syn = function(data){
  output_codons = matrix(0,dim(data)[1],4)
  input_codons = unlist(data$codon)
  colnames(output_codons) = c("A","C","G","T")
  output_codons = data.frame(output_codons)
  data$codon.pos = as.numeric(data$codon.pos)
  for (i in 1:dim(output_codons)[1]){
    #print(i)
    for(j in 1:4){
      tmp = strsplit(input_codons[i],split = "")[[1]]
      tmp[data$codon.pos[i][[1]]] = colnames(output_codons)[j]
      output_codons[i,j] = paste0(tmp,collapse = "")
    }
  }
  output_amino_acid = matrix(0,dim(data)[1],4)
  colnames(output_amino_acid) = c("A.amino","C.amino","G.amino","T.amino")
  syn_flag = output_amino_acid
  plcs = which(is.element(input_codons,row.names(codons_to_amino_acids)))
  for (i in 1:4){
    output_amino_acid[plcs,i] = codons_to_amino_acids[output_codons[plcs,i],2]
    syn_flag[which(data$amino_acid==output_amino_acid[,i] & data$amino_acid!=0),i] = 1
  }
  syn = apply(syn_flag*data[,c("A","C","G","T")],1,sum)
  non_syn = apply( (syn_flag-1)*(-1) * data[,c("A","C","G","T")],1,sum)
  syn[-plcs] = NA
  non_syn[-plcs] = NA
  return(cbind(data,syn,non_syn,output_amino_acid))
}

unlist_table = function(old_table){
  new_table = NULL
  for (i in c(1:dim(old_table)[2])){
    print(i)
    tmp = unlist(old_table[,i])
    print(length(tmp))
    new_table = cbind(new_table,tmp)
  }
  colnames(new_table) = colnames(old_table)
  new_table = data.frame(new_table)
  return(new_table)
}

change_cols_to_numeric = function(numeric_cols,data){
  for (i in 1:length(numeric_cols)){
    data[,numeric_cols[i]] = as.numeric(data[,numeric_cols[i]])
  }
  return(data)
}

add_ACGT_exposure = function(data){
  A_exposure = C_exposure = G_exposure = T_exposure = data$exposure
  A_exposure[which(data$base=="A")] = 0
  C_exposure[which(data$base=="C")] = 0
  G_exposure[which(data$base=="G")] = 0
  T_exposure[which(data$base=="T")] = 0
  data = cbind(data,A_exposure,C_exposure,G_exposure,T_exposure)
  return(data)
}

add_amino_num = function(data,regions){
  tmp = regions[,c("ref_site","amino_num")]
  row.names(tmp) = regions$ref_site
  amino_num = tmp[as.character(data$ref_site),"amino_num"]
  data = cbind(data,amino_num)
  return(data)
}

##############
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
# setwd(paste0(begin_path,"/vars/tables_data_prediction/"))
setwd(paste0(begin_path,"/vars/tables_data/"))
# load("base_outputs")
# load("base_states")
load("codon_states")
load("codon_outputs")

setwd(begin_path)
source("COVID_functions.R")
setwd(paste0(begin_path,"vars/"))
load("site_details")

tmp = strsplit(codon_states,split = "")
bad_rows = NULL
codons_table = matrix(0,length(codon_states),17)
for (i in 1:length(tmp)){
  print(i)
  if (tmp[i][[1]][7]=="s"){ # codon is "missing"
    site = as.numeric(paste0(tmp[i][[1]][11:length(tmp[i][[1]])],collapse = ""))
  }else{
    site = as.numeric(paste0(tmp[i][[1]][7:length(tmp[i][[1]])],collapse = ""))
  }
  if(is.na(site)){
    print(paste("row ",i, " is omitted!"))
    bad_rows = c(bad_rows,i)
  }else{
    codons_table[i,] = unlist(c(tmp[i][[1]][1],tmp[i][[1]][2],tmp[i][[1]][3],
                                paste0(tmp[i][[1]][4:6],collapse = ""),
                                site,regions[regions$site==site,]))
  }
}
colnames(codons_table) = c("L.neighbor","R.neighbor","base","codon","site",
                           colnames(regions))
codons_table = data.frame(codons_table)
missing_codons_plcs = which(codons_table$codon=="mis")
codons_table = cbind(codons_table,codon_outputs)
codons_table = add_CG(codons_table)
codons_table = add_transitions_transversions(codons_table)
codons_table = unlist_table(codons_table)

setwd(paste0(begin_path,"/vars/"))
load("groups")
load("group_sizes")
load("codons_to_amino_acids")
load("syn_group_sizes.1")
load("syn_group_sizes.2")
load("syn_group_sizes.3")
load("non_syn_group_sizes.1")
load("non_syn_group_sizes.2")
load("non_syn_group_sizes.3")
load("syn_ti_const2")
load("syn_tv_const2")
load("non_syn_ti_const2")
load("non_syn_tv_const2")

# codons_table = codons_table[which(codons_table$L.neighbor!="?" & 
#                                     codons_table$R.neighbor!="?"),]
codons_table = add_amino_acid_exposure(codons_table,
                                       codons_to_amino_acids,ratio,group_sizes.1,group_sizes.2,group_sizes.3)

numeric_cols = c("site","ref_site","codon.pos",
                 "A","C","G","T","line","exposure","zero_branch",
                 "CG","transitions","transversions","syn_exposure",
                 "non_syn_exposure")
codons_table = change_cols_to_numeric(numeric_cols,data = codons_table)

rm_cols = which(is.element(colnames(codons_table),c("M","q","s","site.1","UTR")))

codons_table = codons_table[,-rm_cols]
y = codons_table$transitions+codons_table$transversions
codons_table = cbind(codons_table,y)
codons_table = add_syn_non_syn(codons_table)
codons_table = add_amino_num(data = codons_table,regions = regions)
codons_table[which(codons_table$ref_site>=27756 & codons_table$ref_site<=27759),"gene"] = "ORF7a"

# codons_table = add_ACGT_exposure(codons_table)

plcs = which(is.na(codons_table$non_syn))
missing_codons_table = codons_table[plcs,]
codons_table = codons_table[-plcs,]
plcs = which(codons_table$exposure==0)
missing_codons_table = rbind(missing_codons_table,codons_table[plcs,])
codons_table = codons_table[-plcs,]
plcs = which(codons_table$gene=="0")
#codons_table$gene[plcs] = "ORF1ab"
plcs = which(codons_table$R.neighbor=="-" | codons_table$R.neighbor=="?" |
               codons_table$L.neighbor=="-" | codons_table$L.neighbor=="?")
# codons_table = codons_table[-plcs,]
plcs = which(codons_table$base=="0")
plcs = which(is.na(codons_table$syn))
table(codons_table$codon[plcs])
# codons_table = codons_table[-plcs,]
sum((codons_table$transitions+codons_table$transversions)!=codons_table$y)
sum((codons_table$syn+codons_table$non_syn)!=codons_table$y)
bad_plc = which(codons_table$codon=="ATG" & codons_table$codon.pos==1 
                & codons_table$base!="A")
codons_table = codons_table[-bad_plc,]
setwd(paste0(begin_path,"vars/"))
# save(codons_table,file = "codons_table")
codons_table_prediction = codons_table
save(codons_table_prediction,file = "codons_table_prediction_all")
save(missing_codons_table,file = "missing_codons_table_preditcion_all")

###########3
# Get only sites with no subs in the tree dataset:
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("codons_table_prediction_all")
setwd(paste0(begin_path,"vars/prediction/"))
load("new_subs_prediction")
load("rel_ref_sites_no_subs")

codons_table_prediction = 
  codons_table_prediction[which(is.element(codons_table_prediction$ref_site,
                                           c(rel_ref_sites_no_subs,new_subs_prediction$ref_site))),]
problematic_plcs = NULL
diff_seqs_num = matrix(0,dim(codons_table_prediction)[1],1)
for (i in 1:dim(new_subs_prediction)[1]){
  print(i)
  plcs = which(codons_table_prediction$ref_site==new_subs_prediction$ref_site[i])
  codon.pos = codons_table_prediction$codon.pos[plcs[1]]
  diff_seqs_num[plcs] = new_subs_prediction$diff_seqs_num[i]
  if (length(plcs)>0){
    if (codon.pos==1){
      codon = paste0(new_subs_prediction[i,c("before","R1","R2")],collapse = "")
    }else if (codon.pos==2){
      codon = paste0(new_subs_prediction[i,c("L1","before","R1")],collapse = "")
    }else{
      codon = paste0(new_subs_prediction[i,c("L2","L1","before")],collapse = "")
    }
    rel_plc = which(codons_table_prediction$codon==codon & 
                      codons_table_prediction$ref_site==new_subs_prediction$ref_site[i] &
                      codons_table_prediction$L.neighbor==new_subs_prediction$L1[i] &
                      codons_table_prediction$R.neighbor==new_subs_prediction$R1[i])
    if (length(rel_plc)==1){
      codons_table_prediction[rel_plc,new_subs_prediction$after[i]] = 1
    }else{
      print("bug!")
      print(i)
      problematic_plcs = c(problematic_plcs,i)
    }
  }
}

codons_table = codons_table_prediction[,1:23]
codons_table = add_transitions_transversions(codons_table)
codons_table = add_amino_acid_exposure(codons_table,
                                       codons_to_amino_acids,ratio,group_sizes.1,group_sizes.2,group_sizes.3)
numeric_cols = c("site","ref_site","codon.pos",
                 "A","C","G","T","line","exposure","zero_branch",
                 "CG","transitions","transversions","syn_exposure",
                 "non_syn_exposure")
codons_table = change_cols_to_numeric(numeric_cols,data = codons_table)
y = codons_table$transitions+codons_table$transversions
codons_table = cbind(codons_table,y)
codons_table = add_syn_non_syn(codons_table)
codons_table = cbind(codons_table,diff_seqs_num)
setwd(paste0(begin_path,"vars/"))
codons_table_prediction = codons_table
save(codons_table_prediction,file = "codons_table_prediction")

all_sites = (unique(codons_table_prediction$ref_site))
no_subs_sites = (unique(codons_table_prediction$ref_site
                        [which(codons_table_prediction$y==0 & 
                                 codons_table_prediction$diff_seqs_num==0)]))
seq_1_sites = (unique(codons_table_prediction$ref_site
                      [which(codons_table_prediction$y>0 & 
                               codons_table_prediction$diff_seqs_num==1)]))
more_than_1_sites = (unique(codons_table_prediction$ref_site
                            [which(codons_table_prediction$y>0 & 
                                     codons_table_prediction$diff_seqs_num>1)]))
# Sites in which a sub occurred in but the initial state is not in the codons_table
problematic_sites = unique(all_sites[which(is.element(all_sites,no_subs_sites)==FALSE & 
                                             is.element(all_sites,more_than_1_sites)==FALSE & 
                                             is.element(all_sites,seq_1_sites)==FALSE)])
codons_table_prediction = codons_table_prediction[-which(is.element(codons_table_prediction$ref_site,problematic_sites)),]
bad_plc = which(codons_table_prediction$codon=="ATG" & codons_table_prediction$codon.pos==1 
                & codons_table_prediction$base!="A")
codons_table_prediction = codons_table_prediction[-bad_plc,]
save(codons_table_prediction,file = "codons_table_prediction")

####################3
###########3
# Get only sites with no non_syn subs in the tree dataset:
get_non_syn_sites = function(codons_table,regions){
  non_syn_sites = unique(codons_table$ref_site[which(codons_table$non_syn>0)])
  codon.pos_non_syn_sites = regions$codon.pos[non_syn_sites]
  non_syn_sites = c(non_syn_sites,(non_syn_sites[which(codon.pos_non_syn_sites==1)]+1),
                    (non_syn_sites[which(codon.pos_non_syn_sites==1)]+2),
                    (non_syn_sites[which(codon.pos_non_syn_sites==2)]-1),
                    (non_syn_sites[which(codon.pos_non_syn_sites==2)]+1),
                    (non_syn_sites[which(codon.pos_non_syn_sites==3)]-1),
                    (non_syn_sites[which(codon.pos_non_syn_sites==3)]-2))
  return(unique(non_syn_sites))
}

begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("site_details")
load("codons_table")
Lanfear_non_syn_sites = get_non_syn_sites(codons_table = codons_table,regions = regions)
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/"))
load("codons_table")
ncbi_non_syn_sites = get_non_syn_sites(codons_table = codons_table,regions = regions)
non_syn_sites = unique(c(Lanfear_non_syn_sites,ncbi_non_syn_sites))
begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
save(non_syn_sites,file = "non_syn_sites")
load("codons_table_prediction")

plcs = which(is.element(codons_table_prediction$ref_site,non_syn_sites))
codons_table_prediction = codons_table_prediction[-plcs,]
bad_plc = which(codons_table_prediction$codon=="ATG" & codons_table_prediction$codon.pos==1 
                & codons_table_prediction$base!="A")
codons_table_prediction = codons_table_prediction[-bad_plc,]

save(codons_table_prediction,file = "codons_table_prediction_only_syn")
#######################
#####3
# remove sites for which the initial state is not in the codons_table
codons_table = codons_table_prediction
length(unique(codons_table$ref_site[which(codons_table$y>0)]))
length(unique(codons_table$ref_site))
length(unique(codons_table$ref_site[which(codons_table$diff_seqs_num==1 & codons_table$y>0)]))
length(unique(codons_table$ref_site[which(codons_table$diff_seqs_num>1 & codons_table$y>0)]))
length(which(is.element(unique(codons_table$ref_site[which(codons_table$y==0)]),
                        unique(codons_table$ref_site[which(codons_table$y>0)]))==FALSE))
tmp = unique(codons_table$ref_site[which(codons_table$diff_seqs_num>0)])
problematic_sites = tmp[which(is.element(tmp,unique(codons_table$ref_site[which(codons_table$y>0)]))==FALSE)]
codons_table = codons_table[-which(is.element(codons_table$ref_site,problematic_sites)),]
codons_table_prediction = codons_table
save(codons_table_prediction,file = "codons_table_prediction")

################
load("codons_table")

amino_acid_num = regions$amino_acid_num[codons_table$site]
amino_acid_before = codons_to_amino_acids[codons_table$codon,2]
amino_acid_after = codons_table$site*0
codons_table = cbind(codons_table,amino_acid_num,amino_acid_before)


###########3
# #check:
# plcs = NULL
# tmp = strsplit(codons_table$codon,split = "")
# for (i in 1: dim(codons_table)[1]){
#   if(tmp[[i]][codons_table$codon.pos[i]]!=codons_table$base[i]){
#     plcs = c(plcs,i)
#   }
# }

