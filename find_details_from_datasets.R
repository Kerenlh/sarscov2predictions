get_all_reg_hash = function(file_names,length_reg_hash){
  all_reg_hash = hash()
  cat("\n")
  pb <- txtProgressBar(min = 0, max = length(file_names), style = 3, width = 50, char = "=")
  for (i in 1:length(file_names)){
    if (paste0(strsplit(file_names[i],split = "")[[1]][1:9]
               ,collapse = "")=="reg_hash."){ # change to reg_hash.
      load(file_names[i])
      curr_hash_names = names(reg_hash)
      # print(i)
      for (j in 1:length(curr_hash_names)){
        all_reg_hash[[curr_hash_names[j]]] = reg_hash[[curr_hash_names[j]]][1:length_reg_hash]
      }
    }
    setTxtProgressBar(pb, i)
  }
  return(all_reg_hash)
}

find_details_from_datasets = function(curr_folder,folder_num){
  setwd(paste0(begin_path,"vars/unique_regs/",curr_folder))
  file_names = dir()
  load(file_names[1])
  length_reg_hash = length(reg_hash[[names(reg_hash)[1]]])
  cat("\n")
  print("Collects regressions hashes into one list")
  all_reg_hash = get_all_reg_hash(file_names = file_names,length_reg_hash = length_reg_hash)
  # save(all_reg_hash,file = "all_reg_hashes")
  # load(paste0("all_reg_hashes",collapse = ""))
  
  setwd(paste0(begin_path,"vars/reg_hash_plcs/",curr_folder,collapse =""))
  load(paste0("all_reg_hash_plcs_",curr_folder))
  
  setwd(paste0(begin_path,"vars/"))
  load(paste0("dataset_reg_name.",folder_num,collapse = ""))
  reg_output = matrix(0,length(dataset_reg_name),length_reg_hash)
  reg_hash_names = unique(dataset_reg_name)
  # print(length(reg_hash_names))
  cat("\n")
  print("Creates a matrix with the regressions results in rows that fit the dataset matrix.")
  cat("\n")
  pb <- txtProgressBar(min = 0, max = length(reg_hash_names), style = 3, width = 50, char = "=")
  for (j in 1:length(reg_hash_names)){
    # print(j)
    # k=1
    # while(is.element(reg_hash_names[j],eval(parse(text = paste0("reg_hash_plcs_names_",k))))==FALSE){
    #   k = k+1
    # }
    # plcs = eval(parse(text = paste0("all_reg_hash_plcs_",k)))[[reg_hash_names[j] ]]
    # if (length(plcs)>0){
    plcs = all_reg_hash_plcs[reg_hash_names[j]][[1]]
    reg_output[plcs,] = matrix(all_reg_hash[[reg_hash_names[j]]][1:length_reg_hash],
                               length(plcs),dim(reg_output)[2],byrow = TRUE)
    setTxtProgressBar(pb, j)
    # }
  }
  setwd(paste0(begin_path,"vars/unique_regs/"))
  save(reg_output,file = paste0("reg_output_",folder_num,collapse = ""))
  rm(reg_output)
  
  cat("\n")
  print("Saving the details from datasets along with the regression results into a
        matrix called details. This may take a while..")
  setwd(paste0(begin_path,"vars/"))
  load(paste0("dataset.",folder_num,collapse = ""))
  details = dataset[,1:11]
  rm(dataset)
  setwd(paste0(begin_path,"vars/unique_regs/"))
  load(paste0("reg_output_",folder_num,collapse = ""))
  details = cbind(details,matrix(0,dim(details)[1],1),reg_output)
  colnames(details) = c("L.neighbor","R.neighbor","codon",            
                        "codon.pos","gene","mat_peptide","stem_loop","CG","amino_acid",       
                        "base","output","ID","logLik","df","AIC","theta",
                        "NB.converged","rows_num","cols_num","P.GLR","P.test","P.logLik",
                        "P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum",
                        "added_log_lik","added_log_lik.P", "output_0_or_1_row") 
  details = data.frame(details)
  numeric_cols = c("rows_num","cols_num","logLik","P.logLik","df","P.df","theta","AIC","P.AIC")
  for (j in 1:length(numeric_cols)){
    details[,numeric_cols[j]] = as.numeric(details[,numeric_cols[j]])
  }
  setwd(paste0(begin_path,"vars/"))
  save(details,file = paste0("details.",folder_num,collapse = ""))
  
  # Check:
  # tmp = apply(details,1,paste0,collapse = ".")
  # tmp2 = which(duplicated(tmp)==FALSE)
  # details = details[tmp2,]
}


