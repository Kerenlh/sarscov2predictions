find_site_patterns = function(begin_path,ref_seq,seqs,tree,testing_flag){
  start_time = Sys.time()
  if (testing_flag==0){
    print("Taking only sequences that appear in the tree, this may take a while...")
    aligned_seqs = return_only_tree_seqs(seqs,tree)
  }else{
    aligned_seqs = seqs
  }
  rm(seqs)
  seqs_names = names(aligned_seqs)
  num_sites = length(aligned_seqs[1][[1]])
  setwd(paste0(begin_path,"vars/"))
  save(seqs_names, file = "seqs_names")
  save(num_sites, file = "num_sites")
  # save(aligned_seqs,file = "all_aligned_seqs")
  
  #####################
  # Find site patterns:
  #####################
  dir.create("site_patterns")
  setwd(paste0(begin_path,"vars/site_patterns/",collapse = ""))
  num_seqs = length(aligned_seqs)
  seq_length = length(aligned_seqs[[1]])
  max_j = floor(length(aligned_seqs)/1000)
  pb <- txtProgressBar(min = 0, max = (max_j+1), style = 3, width = 50, char = "=")
  cat("\n")
  print(paste("Splitting sequences in groups of 1000, there are",(max_j+1),"groups"))
  for (j in 0:max_j){
    # print(paste("group number:",j))
    split_seqs = matrix("m",1000,length(aligned_seqs[[1]]))
    count = 0
    for (i in (j*1000+1):min(length(aligned_seqs),(j+1)*1000)){
      # print(i)
      count = count+1
      split_seqs[count,] = strsplit(as.character(aligned_seqs[[i]]),split = "")[[1]]
    }
    if (j==max_j){
      split_seqs = split_seqs[1:(num_seqs-j*1000),]
    }
    save(split_seqs,file = paste0("split_seqs_",(j*1000+1),"_",(j+1)*1000,collapse=""))
    setTxtProgressBar(pb, (j+1))
  }
  
  ##################
  # Unite split_seqs: 
  ##################
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
  dir.create("all")
  num_seqs = length(aligned_seqs)
  seq_length = length(aligned_seqs[[1]])
  rm(aligned_seqs)
  count = 0
  num_sites = 5500
  # print(file_names)
  cat("\n")
  print("Splitting sites in groups of 5500, there are 6 groups:")
  pb <- txtProgressBar(min = 0, max = 6, style = 3, width = 50, char = "=")
  for (j in 1:6){
    # print(paste("group number:",j))
    plcs = c((1+(j-1)*num_sites):min(seq_length,(j*num_sites)))
    all_split_seqs = matrix("M",num_seqs,length(plcs))
    count = 0
    for (i in 1:length(file_names)){
      # print(i)
      setwd(paste0(begin_path,"/vars/site_patterns/"))
      load(file_names[i])
      # print(count)
      # print(dim(split_seqs))
      all_split_seqs[(count+1):(dim(split_seqs)[1]+count),] = split_seqs[,plcs]
      count = count + dim(split_seqs)[1]
      rm(split_seqs)
    }
    split_seqs = all_split_seqs
    setwd(paste0(begin_path,"/vars/site_patterns/all"))
    save(split_seqs,file = paste0("all_split_seqs_",plcs[length(plcs)],collapse =""))
    setTxtProgressBar(pb, j)
  }
  ###########
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
  setwd(paste0(begin_path,"vars/site_patterns/",collapse = ""))
  dir.create("all2")
  cat("\n")
  print("Splitting sites in groups of 5500 - second stage, there are 6 groups:")
  for (j in 1:6){
    setwd(paste0(begin_path,"vars/site_patterns/all",collapse = ""))
    load(file_names[j])
    site_patterns = site_patterns_plcs = no_subs_plcs = NULL
    for (i in 1:dim(split_seqs)[2]){
      # print(i)
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
    setTxtProgressBar(pb, j)
  }
  
  ######33
  # Unite site_patterns:
  cat("\n")
  print("Uniting site patterns and removing temporary files")
  setwd(paste0(begin_path,"vars/site_patterns/all2",collapse = ""))
  file_names = dir()
  all_site_patterns = NULL
  for (i in 7:12){
    # print(i)
    load(file_names[i])
    all_site_patterns = c(all_site_patterns,site_patterns)
  }
  site_patterns = all_site_patterns
  
  all_site_patterns_plcs = NULL
  for (i in 13:18){
    # print(i)
    load(file_names[i])
    all_site_patterns_plcs = c(all_site_patterns_plcs,site_patterns_plcs)
  }
  site_patterns_plcs = all_site_patterns_plcs
  
  all_no_subs_plcs = NULL
  for (i in 1:6){
    # print(i)
    load(file_names[i])
    all_no_subs_plcs = c(all_no_subs_plcs,no_subs_plcs)
  }
  no_subs_plcs = all_no_subs_plcs
  
  setwd(paste0(begin_path,"vars/",collapse = ""))
  save(site_patterns,file = "all_site_patterns")
  save(site_patterns_plcs,file = "all_site_patterns_plcs")
  save(no_subs_plcs,file = "all_no_subs_plcs")
  site_patterns_numbers = site_patterns_plcs[which(duplicated(site_patterns)==FALSE)]
  save(site_patterns_numbers,file = "site_patterns_numbers")
  
  unlink("site_patterns", recursive = TRUE)
  print(Sys.time()-start_time)
}

#########################
# Testing functions: #
#########################
get_no_subs_plcs_with_indels = function(site_patterns_plcs,site_patterns,no_subs_plcs){
  cat("\n")
  print("Find sites in which there were no subsitutions when indels are disregarded")
  cat("\n")
  pb <- txtProgressBar(min = 0, max = length(site_patterns_plcs), style = 3, width = 50, char = "=")
  more_no_subs_plcs =bugs =  NULL
  for (i in 1:length(site_patterns_plcs)){
    # print(i)
    curr_alignment_site = site_patterns_plcs[i]
    split_site_pattern = strsplit(site_patterns[which(site_patterns_plcs==curr_alignment_site)],split = "")[[1]]
    table_site_pattern = table(split_site_pattern)
    # table_site_pattern = table_site_pattern[which(is.element(names(table_site_pattern),c("A","C","G","T")))]
    if (length(table_site_pattern)==1){
      # print("bug!")
      bugs = c(bugs,i)
      break
    }
    if (length(table_site_pattern)==2 & is.element("-",names(table_site_pattern))){
      # print(table_site_pattern)
      more_no_subs_plcs = c(more_no_subs_plcs,site_patterns_plcs[i])
    }
    setTxtProgressBar(pb, i)
  }
  
  no_subs_plcs_with_indels = c(no_subs_plcs,more_no_subs_plcs)
  return(no_subs_plcs_with_indels)
}

get_testing_subs = function(begin_path,no_subs_plcs,site_patterns_plcs,site_patterns,site_to_alignment,
                            new_seqs,codons_table,site_details,
                            site_to_alignment_testing,no_subs_plcs_testing,
                            site_patterns_plcs_testing,site_patterns_testing){
  
  no_subs_plcs_with_indels = get_no_subs_plcs_with_indels(site_patterns_plcs,site_patterns,no_subs_plcs)
  
  subs_ref_sites = unique(codons_table$ref_site[which(codons_table$y>0)])
  no_subs_plcs_with_indels_ref_sites = 
    site_to_alignment$site[no_subs_plcs_with_indels]
  codon_ref_sites = unique(regions$ref_site[which(regions$codon.pos>0)])
  no_subs_ref_sites = codon_ref_sites[which(is.element(codon_ref_sites,subs_ref_sites)==FALSE)]
  no_subs_ref_sites = no_subs_ref_sites[which(is.element(no_subs_ref_sites,no_subs_plcs_with_indels_ref_sites)==TRUE)]
  no_subs_plcs_testing_ref_sites = site_to_alignment_testing$site[no_subs_plcs_testing]
  rel_sites_no_subs = no_subs_ref_sites[which(is.element(no_subs_ref_sites,no_subs_plcs_testing_ref_sites))]
  rel_sites_subs = no_subs_ref_sites[which(is.element(no_subs_ref_sites,no_subs_plcs_testing_ref_sites)==FALSE)]
  
  subs = problematic_site_patterns = problematic_site_patterns_2_vals = NULL
  for (i in 1:length(rel_sites_subs)){
    # print(i)
    curr_alignment_site = site_to_alignment_testing$alignment_num[which(site_to_alignment_testing$site==rel_sites_subs[i])]
    split_site_pattern = strsplit(site_patterns_testing[which(site_patterns_plcs_testing==curr_alignment_site)],split = "")[[1]]
    table_site_pattern = table(split_site_pattern)
    table_site_pattern = table_site_pattern[which(is.element(names(table_site_pattern),c("A","C","G","T")))]
    if (length(table_site_pattern)==2 & min(table_site_pattern)==1){
      # print(i)
      diff_seqs_num = 1
      before = names(table_site_pattern)[which.max(table_site_pattern)]
      after = names(table_site_pattern)[which.min(table_site_pattern)]
      sub_seq = strsplit(as.character(new_seqs[which(split_site_pattern==after)]),split = "")[[1]]
      if (sub_seq[curr_alignment_site]!=after){
        # print("bug!")
      }
      subs = rbind(subs,c(before,sub_seq[(curr_alignment_site-2):(curr_alignment_site+2)],
                          rel_sites_subs[i],curr_alignment_site,diff_seqs_num))
    }else if((length(table_site_pattern)==2 & min(table_site_pattern)>1)){
      before = names(table_site_pattern)[which.max(table_site_pattern)]
      after = names(table_site_pattern)[which.min(table_site_pattern)]
      diff_seqs_num = table_site_pattern[after]
      sub_seqs = NULL
      for (k in 1:table_site_pattern[after]){
        curr_sub_seq = strsplit(as.character(new_seqs[which(split_site_pattern==after)[k]]),split = "")[[1]]
        sub_seqs = rbind(sub_seqs,curr_sub_seq[(curr_alignment_site-2):(curr_alignment_site+2)])
      }
      problem = 0
      for (k in 1:dim(sub_seqs)[2]){
        if(length(unique(sub_seqs[,k]))>1){
          problem = 1
        }
      }
      if (problem){
        tmp = tmp = c(0,0,0,0)
        names(tmp) = c("A","C","G","T")
        for (j in 1:length(table_site_pattern)){
          tmp[names(table_site_pattern)[j]] = table_site_pattern[j]
        }
        problematic_site_patterns_2_vals = rbind(problematic_site_patterns_2_vals,
                                                 c(rel_sites_subs[i],tmp))
      }else{
        sub_seq = sub_seqs[1,]
        subs = rbind(subs,c(before,sub_seq,
                            rel_sites_subs[i],curr_alignment_site,diff_seqs_num))
      }
    }else{
      tmp = c(0,0,0,0)
      names(tmp) = c("A","C","G","T")
      for (j in 1:length(table_site_pattern)){
        tmp[names(table_site_pattern)[j]] = table_site_pattern[j]
      }
      problematic_site_patterns = rbind(problematic_site_patterns,
                                        c(rel_sites_subs[i],length(table_site_pattern),tmp))
      # print("problematic site pattern:")
      # print(rel_sites_subs[i])
      # print(table_site_pattern)
    }
  }
  
  # problematic_site_patterns = data.frame(problematic_site_patterns)
  # colnames(problematic_site_patterns) = c("ref_site","num_vals","A","C","G","T")
  # problematic_site_patterns_2_vals = data.frame(problematic_site_patterns_2_vals)
  # # sites in which there is more than 1 different seq and the seqs are not the same in enviornment
  # colnames(problematic_site_patterns_2_vals) = c("ref_site","A","C","G","T")
  colnames(subs) = c("before","L2","L1","after","R1","R2","ref_site","aln_site","diff_seqs_num")
  subs = data.frame(subs)
  
  setwd(paste0(begin_path,"testing/vars/"))
  new_subs_testing = subs
  save(new_subs_testing,file = "new_subs_testing")
  # new_subs_testing = only sites with 2 values in all seqs.
  
  rel_ref_sites_no_subs = rel_sites_no_subs
  save(rel_ref_sites_no_subs,file = "rel_ref_sites_no_subs")
  
  # rel_ref_sites_no_subs = c(rel_sites_no_subs,
  #                           problematic_site_patterns$ref_site[which(problematic_site_patterns$num_vals==1)])
  
  # new_subs_sites = c(subs$ref_site,problematic_site_patterns$ref_site[which(problematic_site_patterns$num_vals>1)])
  # new_subs_sites = cbind(new_subs_sites,c(matrix(2,dim(subs)[1],1),
  #                                         problematic_site_patterns$num_vals[which(problematic_site_patterns$num_vals>1)]))
  # new_subs_sites = data.frame(new_subs_sites)
  # colnames(new_subs_sites) = c("ref_site","num_vals")
  # new_subs_sites$ref_site = as.numeric(new_subs_sites$ref_site)
  # new_subs_sites$num_vals = as.numeric(new_subs_sites$num_vals)
  # save(new_subs_sites,file = "new_subs_sites")
  # # new_subs_sites = sites with more than 1 value in all seqs
}
