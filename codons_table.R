codon_names = c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG",
                "ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC",
                "CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA",
                "GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT",
                "TTA","TTC","TTG","TTT")
group_names = c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Start/Met",
                "Leu","Lys","Phe","Pro","Ser","Thr","Trp","Tyr","Val","Stop")
group_names2 = c("A","R","N","D","C","Q","E","G","H","I","M",
                 "L","K","F","P","S","T","W","Y","V","Stop")

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

codons_to_amino_acids = matrix(0,64,2)
row.names(codons_to_amino_acids) = codon_names
for (j in 1:length(group_names)){
  gr = paste("g",j,sep="")
  curr.g = eval(parse(text = gr))
  codons_to_amino_acids[curr.g,1] = group_names[j]
  codons_to_amino_acids[curr.g,2] = group_names2[j]
}
setwd(paste0(begin_path,"vars/"))
save(codons_to_amino_acids,file = "codons_to_amino_acids")


# Functions:
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

add_syn_non_syn = function(data){
  output_codons = matrix(0,dim(data)[1],4)
  input_codons = unlist(data$codon)
  colnames(output_codons) = c("A","C","G","T")
  output_codons = data.frame(output_codons)
  data$codon.pos = as.numeric(data$codon.pos)
  cat("\n")
  print("Adding syn, non_syn, output_amino columns")
  cat("\n")
  pb <- txtProgressBar(min = 0, max = dim(output_codons)[1], style = 3, width = 50, char = "=")
  for (i in 1:dim(output_codons)[1]){
    #print(i)
    for(j in 1:4){
      tmp = strsplit(input_codons[i],split = "")[[1]]
      tmp[data$codon.pos[i][[1]]] = colnames(output_codons)[j]
      output_codons[i,j] = paste0(tmp,collapse = "")
    }
    setTxtProgressBar(pb, i)
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

find_syn_non_syn_plcs = function(codons_table,before,after,syn_flag){
  if (syn_flag){
    plcs = which(codons_table$amino_acid==codons_table[,paste0(after,".amino",collapse = "")] &
                   codons_table$base==before)
  }else{
    plcs = which(codons_table$amino_acid!=codons_table[,paste0(after,".amino",collapse = "")] &
                   codons_table$base==before)
  }
  output = codons_table[plcs,after]
  return(cbind(codons_table[plcs,],output))
}

get_syn_non_syn_exposure_data = function(codons_table,syn_flag){
  tmp1 = find_syn_non_syn_plcs(codons_table = codons_table,before = "A",after = "G",syn_flag = syn_flag)
  tmp2 = find_syn_non_syn_plcs(codons_table = codons_table,before = "G",after = "A",syn_flag = syn_flag)
  tmp3 = find_syn_non_syn_plcs(codons_table = codons_table,before = "C",after = "T",syn_flag = syn_flag)
  tmp4 = find_syn_non_syn_plcs(codons_table = codons_table,before = "T",after = "C",syn_flag = syn_flag)
  ti = rbind(tmp1,tmp2,tmp3,tmp4)
  
  tmp5 = find_syn_non_syn_plcs(codons_table = codons_table,before = "A",after = "C",syn_flag = syn_flag)
  tmp6 = find_syn_non_syn_plcs(codons_table = codons_table,before = "A",after = "T",syn_flag = syn_flag)
  tmp7 = find_syn_non_syn_plcs(codons_table = codons_table,before = "C",after = "A",syn_flag = syn_flag)
  tmp8 = find_syn_non_syn_plcs(codons_table = codons_table,before = "C",after = "G",syn_flag = syn_flag)
  tmp9 = find_syn_non_syn_plcs(codons_table = codons_table,before = "G",after = "C",syn_flag = syn_flag)
  tmp10 = find_syn_non_syn_plcs(codons_table = codons_table,before = "G",after = "T",syn_flag = syn_flag)
  tmp11 = find_syn_non_syn_plcs(codons_table = codons_table,before = "T",after = "A",syn_flag = syn_flag)
  tmp12 = find_syn_non_syn_plcs(codons_table = codons_table,before = "T",after = "G",syn_flag = syn_flag)
  tv = rbind(tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12)
  
  data = rbind(ti,tv)
  ti_flag = matrix(0,dim(data)[1],1)
  ti_flag[1:dim(ti)[1],1] = 1
  data = cbind(data,ti_flag)
  return(data)
}

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

add_amino_acid_exposure = function(data,codons_to_amino_acids,syn_group_sizes.1,
                                   syn_group_sizes.2,syn_group_sizes.3,
                                   non_syn_group_sizes.1,
                                   non_syn_group_sizes.2,non_syn_group_sizes.3,
                                   non_syn_tv_const){
  amino_acid = syn_exposure =  non_syn_exposure = matrix(0,dim(data)[1],1)
  # plcs = which(is.element(data$codon,rownames(codons_to_amino_acids)))
  # no_valid_codon_plcs = which(is.element(data$codon,rownames(codons_to_amino_acids)==FALSE))
  for (i in 1:dim(codons_to_amino_acids)[1]){
    # print(i)
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

unlist_table = function(old_table){
  new_table = NULL
  for (i in c(1:dim(old_table)[2])){
    # print(i)
    tmp = unlist(old_table[,i])
    # print(length(tmp))
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

get_tv_ti_ratio = function(data){
  syn_data = get_syn_non_syn_exposure_data(codons_table = data,syn_flag = 1)
  non_syn_data = get_syn_non_syn_exposure_data(codons_table = data,syn_flag = 0)
  curr_model_names = output~ti_flag+offset(log(exposure))
  
  syn.model.p = glm(curr_model_names,data = syn_data,family = "poisson")
  table(exp(predict(syn.model.p,syn_data))/syn_data$exposure)
  non_syn.model.p = glm(curr_model_names,data = non_syn_data,family = "poisson")
  table(exp(predict(non_syn.model.p,non_syn_data))/non_syn_data$exposure)
  
  tmp = unique(exp(predict(syn.model.p,syn_data))/syn_data$exposure)
  syn_ti_const = tmp[1]
  syn_tv_const = tmp[length(tmp)]
  syn_tv_const = syn_tv_const/syn_ti_const
  syn_ti_const = 1
  
  tmp = unique(exp(predict(non_syn.model.p,non_syn_data))/non_syn_data$exposure)
  non_syn_ti_const = tmp[1]
  non_syn_tv_const = tmp[length(tmp)]
  non_syn_tv_const = non_syn_tv_const/non_syn_ti_const
  non_syn_ti_const = 1
  
  setwd(paste0(begin_path,"vars/"))
  save(syn_ti_const,file = "syn_ti_const")
  save(syn_tv_const,file = "syn_tv_const")
  save(non_syn_ti_const,file = "non_syn_ti_const")
  save(non_syn_tv_const,file = "non_syn_tv_const")
  return(c(syn_ti_const,syn_tv_const,non_syn_ti_const,non_syn_tv_const))
}

add_amino_num_to_regions = function(regions,codons_to_amino_acids){
  gene_names = unique(regions$gene)
  amino_acid_num = matrix("0",dim(regions)[1],1)
  for (i in 1:length(gene_names)){
    if (gene_names[i]!="0"){
      # print(i)
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
          # print(j)
        }else{
          while(regions$codon.pos[gene_plcs[j]]!=1){
            amino_acid_num[gene_plcs[j]] = count
            j = j+1
            # print(j)
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
  return(regions)
}

add_amino_num = function(data,regions){
  tmp = regions[,c("ref_site","amino_acid_num")]
  tmp = tmp[-which(tmp$ref_site==0),]
  row.names(tmp) = tmp$ref_site
  amino_num = tmp[as.character(data$ref_site),"amino_acid_num"]
  data = cbind(data,amino_num)
  return(data)
}

##############

create_codons_table = function(codon_states,codon_outputs,regions,testing_flag){
  tmp = strsplit(codon_states,split = "")
  bad_rows = NULL
  codons_table = matrix(0,length(codon_states),(5+dim(regions)[2]))
  cat("\n")
  print("Creating initial codon_table")
  cat("\n")
  pb <- txtProgressBar(min = 0, max = length(tmp), style = 3, width = 50, char = "=")
  for (i in 1:length(tmp)){
    # print(i)
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
    setTxtProgressBar(pb, i)
  }
  colnames(codons_table) = c("L.neighbor","R.neighbor","base","codon","site",
                             colnames(regions))
  codons_table = data.frame(codons_table)
  codons_table = cbind(codons_table,codon_outputs)
  plcs_to_remove = which(codons_table$codon=="mis" | codons_table$exposure==0)
  codons_table = codons_table[-plcs_to_remove,]
  rm_cols = which(is.element(colnames(codons_table),c("M","q","s","site.1","UTR")))
  codons_table = codons_table[,-rm_cols]
  amino_acid = codons_to_amino_acids[codons_table$codon,2]
  ref_amino_acid = codons_to_amino_acids[codons_table$ref_codon,2]
  codons_table = cbind(codons_table,amino_acid,ref_amino_acid)
  codons_table = add_syn_non_syn(codons_table)
  codons_table = add_CG(codons_table)
  codons_table = add_transitions_transversions(codons_table)
  if(testing_flag){
    setwd(paste0(begin_path,"vars/"))
    load("syn_tv_const"); load("non_syn_tv_const")
  }else{
    tmp = get_tv_ti_ratio(codons_table)
    syn_ti_const = tmp[1]; syn_tv_const = tmp[2]; non_syn_ti_const = tmp[3]; non_syn_tv_const = tmp[4]
  }
  codons_table = unlist_table(codons_table)
  
  tmp = get_group_sizes(ratio=syn_tv_const,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21)
  syn_group_sizes.1 = tmp[[1]];syn_group_sizes.2 = tmp[[2]];syn_group_sizes.3 = tmp[[3]]
  tmp = get_group_sizes(ratio=non_syn_tv_const,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21)
  non_syn_group_sizes.1 = tmp[[1]];non_syn_group_sizes.2 = tmp[[2]];non_syn_group_sizes.3 = tmp[[3]]
  
  codons_table = add_amino_acid_exposure(codons_table,codons_to_amino_acids,syn_group_sizes.1,
                                         syn_group_sizes.2,syn_group_sizes.3,
                                         non_syn_group_sizes.1,non_syn_group_sizes.2,
                                         non_syn_group_sizes.3,non_syn_tv_const)
  
  numeric_cols = c("site","ref_site","codon.pos",
                   "A","C","G","T","line","exposure","zero_branch",
                   "CG","transitions","transversions","syn_exposure",
                   "non_syn_exposure","syn","non_syn")
  codons_table = change_cols_to_numeric(numeric_cols,data = codons_table)
  
  y = codons_table$transitions+codons_table$transversions
  codons_table = cbind(codons_table,y)
  regions = add_amino_num_to_regions(regions,codons_to_amino_acids)
  codons_table = add_amino_num(data = codons_table,regions = regions)
  
  # Check:
  # plcs = which(is.na(codons_table$non_syn))
  # missing_codons_table = codons_table[plcs,]
  # codons_table = codons_table[-plcs,]
  # plcs = which(codons_table$exposure==0)
  # missing_codons_table = rbind(missing_codons_table,codons_table[plcs,])
  # codons_table = codons_table[-plcs,]
  # plcs = which(codons_table$gene=="0")
  #codons_table$gene[plcs] = "ORF1ab"
  # plcs = which(codons_table$R.neighbor=="-" | codons_table$R.neighbor=="?" |
  #                codons_table$L.neighbor=="-" | codons_table$L.neighbor=="?")
  # codons_table = codons_table[-plcs,]
  # plcs = which(codons_table$base=="0")
  # plcs = which(is.na(codons_table$syn))
  # table(codons_table$codon[plcs])
  # codons_table = codons_table[-plcs,]
  # bad_plc = which(codons_table$codon=="ATG" & codons_table$codon.pos==1 
  #                 & codons_table$base!="A")
  # codons_table = codons_table[-bad_plc,]
  sum((codons_table$transitions+codons_table$transversions)!=codons_table$y)
  sum((codons_table$syn+codons_table$non_syn)!=codons_table$y)
  
  codons_table = codons_table[,c("L.neighbor","R.neighbor","base","codon",
                                 "site","ref_site","codon.pos","gene","mat_peptide","stem_loop",       
                                 "notes","notes2","ref_seq","ref_codon","ref_amino_acid",
                                 "A","C","G","T","line","exposure","zero_branch","CG",
                                 "transitions","transversions","amino_acid","syn_exposure",
                                 "non_syn_exposure","y","syn","non_syn","A.amino","C.amino",
                                 "G.amino","T.amino","amino_num")]
  return(codons_table)
}

create_codons_table_testing = function(codons_table,new_subs_testing,rel_ref_sites_no_subs,
                                       syn_tv_const,non_syn_tv_const){
  codons_table = codons_table[which(is.element(codons_table$ref_site,
                                               c(rel_ref_sites_no_subs,new_subs_testing$ref_site))),]
  problematic_plcs = NULL
  diff_seqs_num = matrix(0,dim(codons_table)[1],1)
  cat("\n")
  pb <- txtProgressBar(min = 0, max = dim(new_subs_testing)[1], style = 3, width = 50, char = "=")
  for (i in 1:dim(new_subs_testing)[1]){
    # print(i)
    plcs = which(codons_table$ref_site==new_subs_testing$ref_site[i])
    codon.pos = codons_table$codon.pos[plcs[1]]
    diff_seqs_num[plcs] = new_subs_testing$diff_seqs_num[i]
    if (length(plcs)>0){
      if (codon.pos==1){
        codon = paste0(new_subs_testing[i,c("before","R1","R2")],collapse = "")
      }else if (codon.pos==2){
        codon = paste0(new_subs_testing[i,c("L1","before","R1")],collapse = "")
      }else{
        codon = paste0(new_subs_testing[i,c("L2","L1","before")],collapse = "")
      }
      rel_plc = which(codons_table$codon==codon & 
                        codons_table$ref_site==new_subs_testing$ref_site[i] &
                        codons_table$L.neighbor==new_subs_testing$L1[i] &
                        codons_table$R.neighbor==new_subs_testing$R1[i])
      if (length(rel_plc)==1){
        codons_table[rel_plc,new_subs_testing$after[i]] = 1
      }else{
        # print("bug!")
        # print(i)
        problematic_plcs = c(problematic_plcs,i)
      }
    }
    setTxtProgressBar(pb, i)
  }
  
  amino_num = codons_table[,36]
  codons_table = codons_table[,c(1:23)]
  codons_table = add_transitions_transversions(codons_table)
  tmp = get_group_sizes(ratio=syn_tv_const,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21)
  syn_group_sizes.1 = tmp[[1]];syn_group_sizes.2 = tmp[[2]];syn_group_sizes.3 = tmp[[3]]
  tmp = get_group_sizes(ratio=non_syn_tv_const,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21)
  non_syn_group_sizes.1 = tmp[[1]];non_syn_group_sizes.2 = tmp[[2]];non_syn_group_sizes.3 = tmp[[3]]
  
  codons_table = add_amino_acid_exposure(codons_table,codons_to_amino_acids,syn_group_sizes.1,
                                         syn_group_sizes.2,syn_group_sizes.3,
                                         non_syn_group_sizes.1,non_syn_group_sizes.2,
                                         non_syn_group_sizes.3,non_syn_tv_const)
  
  numeric_cols = c("site","ref_site","codon.pos",
                   "A","C","G","T","line","exposure","zero_branch",
                   "CG","transitions","transversions","syn_exposure",
                   "non_syn_exposure")
  codons_table = change_cols_to_numeric(numeric_cols,data = codons_table)
  y = codons_table$transitions+codons_table$transversions
  codons_table = cbind(codons_table,y)
  codons_table = add_syn_non_syn(codons_table)
  codons_table = cbind(codons_table,diff_seqs_num,amino_num)
  codons_table$mat_peptide[which(codons_table$mat_peptide!="0")] = 1
  codons_table$mat_peptide = as.numeric(codons_table$mat_peptide)
  
  all_sites = (unique(codons_table$ref_site))
  no_subs_sites = (unique(codons_table$ref_site
                          [which(codons_table$y==0 & 
                                   codons_table$diff_seqs_num==0)]))
  seq_1_sites = (unique(codons_table$ref_site
                        [which(codons_table$y>0 & 
                                 codons_table$diff_seqs_num==1)]))
  more_than_1_sites = (unique(codons_table$ref_site
                              [which(codons_table$y>0 & 
                                       codons_table$diff_seqs_num>1)]))
  # Sites in which a sub occurred in but the initial state is not in the codons_table
  problematic_sites = unique(all_sites[which(is.element(all_sites,no_subs_sites)==FALSE & 
                                               is.element(all_sites,more_than_1_sites)==FALSE & 
                                               is.element(all_sites,seq_1_sites)==FALSE)])
  codons_table = codons_table[-which(is.element(codons_table$ref_site,problematic_sites)),]
  return(codons_table)
}