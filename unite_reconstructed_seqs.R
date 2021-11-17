sort_file_names = function(file_names){
  tmp = unlist(strsplit(file_names,split = "\\."))
  nums = as.numeric(tmp[seq(2,length(tmp),2)])
  file_names = file_names[order(nums)]
  return(file_names)
}

unite_reconstructed_seqs = function(begin_path, site_patterns_plcs,no_subs_plcs,num_sites,
                                    init_seq_no_subs_plcs){
  setwd(trees_path)
  file_names = dir()
  nums = NULL
  for (i in 1:length(file_names)){
    if (strsplit(file_names[i],split = "[.]")[[1]][1]=="tree_list"){
      nums = c(nums,strsplit(file_names[i],split = "[.]")[[1]][2])
    }else{
      nums = c(nums,0)
    }
  }
  nums = as.numeric(nums)
  load(paste0("tree_list.",max(nums),collapse = ""))
  
  #####
  # Add duplicated site patterns:
  site_patterns_numbers = site_patterns_plcs[which(duplicated(site_patterns)==FALSE)]
  site_patterns_plcs = cbind(site_patterns_plcs,matrix(0,length(site_patterns_plcs),1))
  colnames(site_patterns_plcs) = c("site","tree_list_num")
  site_patterns_plcs = data.frame(site_patterns_plcs)
  cat("\n")
  print("Some sites have identical site patterns (identical letters in all sequences for the specific sites)")
  print("The ASR is done once for identical site patterns so here we add the ASR results for identical site patterns.")
  cat("\n")
  # pb <- txtProgressBar(min = 0, max = length(unique_site_patterns), style = 3, width = 50, char = "=")
  cat("\n")
  for (i in 1:length(unique_site_patterns)){
    site_patterns_plcs$tree_list_num[which(site_patterns==unique_site_patterns[i])] =
      site_patterns_numbers[i]
    # setTxtProgressBar(pb, i)
  }
  ##########
  pb <- txtProgressBar(min = 0, max = num_sites, style = 3, width = 50, char = "=")
  cat("\n")
  print(paste("Reconstructing sequences at internal nodes, going over all sites in groups of 1000, takes ~40 minutes"))
  cat("\n")
  # Initializing partial sequences with "M" at all sites. For
  for (j in 1:ceiling(num_sites/1000)){
    setwd(trees_path)
    start_time = Sys.time()
    if (j==ceiling(num_sites/1000)){
      part_seqs = data.frame(matrix("M",length(bases),(num_sites- 1000*floor(num_sites/1000)) ) )
    }else{
      part_seqs = data.frame(matrix("M",length(bases),1000))
    }
    rownames(part_seqs) = names(bases)
    rel_plcs = c(((j-1)*1000+1):(j*1000))
    rel_file_plcs = which(is.element(nums,rel_plcs))
    # loading and putting in place relevant tree_list.# files
    if (length(rel_file_plcs)>0){
      for (i in 1:length(rel_file_plcs)){
        # print(i)
        load(file_names[rel_file_plcs[i]])
        curr_num = as.numeric(strsplit(file_names[rel_file_plcs[i]],split = "[.]")[[1]][2])
        rel_sites = site_patterns_plcs$site[which(site_patterns_plcs$tree_list_num==curr_num)] 
        rel_sites = rel_sites[which((is.element(rel_sites,rel_plcs)))]
        part_seqs[names(bases),(rel_sites-((j-1)*1000))] = bases
        setTxtProgressBar(pb, ((j-1)*1000+i))
      }
    }
    # Adding the reference sequence letter to sites that were identical in all sequences
    rel_no_subs_plcs = which(is.element(no_subs_plcs,rel_plcs))
    if (length(rel_no_subs_plcs)>0){
      for(i in 1:length(rel_no_subs_plcs)){
        # print(i)
        part_seqs[,(no_subs_plcs[rel_no_subs_plcs[i]]-((j-1)*1000))] = 
          init_seq_no_subs_plcs[rel_no_subs_plcs[i]]
      }
    }
    # Adding tree_list.# files information for sites whose relevent site_pattern is in a different 1000 segment
    missing_sites = which(part_seqs[1,]=="M")
    if (length(missing_sites)>0){
      for(i in 1:length(missing_sites)){
        curr_num = site_patterns_plcs$tree_list_num[which(site_patterns_plcs$site==missing_sites[i]+(j-1)*1000)]
        load(paste0("tree_list.",curr_num))
        rel_sites = site_patterns_plcs$site[which(site_patterns_plcs$tree_list_num==curr_num)] 
        rel_sites = rel_sites[which((is.element(rel_sites,rel_plcs)))]
        part_seqs[names(bases),(rel_sites-((j-1)*1000))] = bases
      }
    }
    missing_sites = which(part_seqs[1,]=="M")
    if (length(missing_sites)>0){
      print("bug! there are missing sites!")
    }
    setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
    save(part_seqs,file = paste0("part_seqs.",j))
    # print(Sys.time()-start_time)
  }
  
  # setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  # file_names = dir()
  # file_names = sort_file_names(file_names)
  # cat("\n")
  # print("Make sure file names are sorted!") 
  # print(file_names) 
  # cat("\n")
  # pb <- txtProgressBar(min = 0, max = 30, style = 3, width = 50, char = "=")
  # cat("\n")
  # print("Some of the tree_list.# files were saved with an addition of .file_part")
  # print("Here we add results frome these files.")
  # for (i in 1:30){
  #   setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  #   # print(file_names[i])
  #   start_time = Sys.time()
  #   load(file_names[i])
  #   missing_plcs_local = which(part_seqs[1,]=="M")
  #   missing_plcs = missing_plcs_local +(i-1)*1000
  #   setwd(trees_path)
  #   # print(length(missing_plcs))
  #   for (j in 1:length(missing_plcs)){
  #     rel_tree_list_num = 
  #       site_patterns_plcs$tree_list_num[which(site_patterns_plcs$site==missing_plcs[j])]
  #     tree_name = paste0("tree_list.",rel_tree_list_num)
  #     if(file.exists(tree_name)){
  #       load(tree_name)
  #     }else{
  #       load(paste0(tree_name,".filepart"))
  #     }
  #     part_seqs[names(bases),missing_plcs_local[j]] = bases
  #     rm(bases)
  #   }
  #   setTxtProgressBar(pb, i)
  #   # print(Sys.time()-start_time)
  #   # print(length(which(part_seqs[1,]=="M")))
  #   setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  #   save(part_seqs,file = paste0("P",file_names[i]))
  # }
  
  setwd(paste0(begin_path,"vars/reconstructed_seqs/"))
  # It's important that this directory will contain only the files "part_seqs.(1-30)"
  file_names = dir()
  file_names = sort_file_names(file_names)
  all_seqs = NULL
  pb <- txtProgressBar(min = 0, max = length(file_names), style = 3, width = 50, char = "=")
  cat("\n")
  print("Uniting partial sequences")
  cat("\n")
  for (i in 1:length(file_names)){
    # print(i)
    load(file_names[i])
    tmp = apply(part_seqs,1,paste0,collapse = "")
    rm(part_seqs)
    all_seqs = apply(cbind(all_seqs,tmp),1,paste0,collapse = "")
    setTxtProgressBar(pb, i)
  }
  setwd(paste0(begin_path,"vars/"))
  save(all_seqs, file = "all_seqs")
  unlink("reconstructed_seqs", recursive = TRUE)
}
