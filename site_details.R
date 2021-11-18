create_site_details = function(num_sites,site_to_alignment){
   # Reference sequence regions:
  ref_seq = "NC_045512"
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
  ref_regions$gene[27756:27759] = "ORF7a"
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
  
  # Sites 13468:13484 have double codon meaning - it should be 13468:21555
  # but it overlaps with 266:13483 so 13468 is kept as codon.pos=3 but
  # according to 13468:21555 it should be 1.
  # Sites 27756:27759 have the same problem (27394:27759,27756:27887)
  
  ###########
  regions = matrix(0,num_sites,9)
  colnames(regions) = c("site",colnames(ref_regions))
  regions = data.frame(regions)
  regions[which(site_to_alignment$site!=0),2:9]=ref_regions[site_to_alignment$site[which(site_to_alignment$site!=0)],]
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
      # print(paste("k=",k))
      k= k+1
    }
  }
  
  ref_codon = matrix(0,dim(regions)[1],1)
  ref_codon[curr_codon_sites] = codons
  regions = cbind(regions,ref_codon)
  return(regions)
}