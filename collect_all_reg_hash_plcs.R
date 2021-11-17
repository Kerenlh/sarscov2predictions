collect_all_reg_hash_plcs = function(begin_path,curr_folder){
  setwd(paste0(begin_path,"vars/reg_hash_plcs/",curr_folder,collapse = ""))
  file_names = dir()
  all_reg_hash_plcs = list()
  problematic_files = NULL
  bad_files = NULL
  for (i in 1:length(file_names)){
    if(strsplit(file_names[i],"[.]")[[1]][1]=="reg_hashes_plcs"){
      #   tmp = tryCatch(load(file_names[i]),
      #            error = function(e) {})
      #   if (is.null(tmp)){
      #     bad_files = c(bad_files,i)
      #     print(bad_files)
      #   }
      #   else if (tmp!="reg_hashes_plcs"){
      #     bad_files = c(bad_files,i)
      #     print(bad_files)
      #   }
      load(file_names[i])
      # print(i)
      curr_names = names(reg_hashes_plcs)
      count = 0
      for(j in 1:length(curr_names)){
        all_reg_hash_plcs[[curr_names[j]]] = reg_hashes_plcs[[curr_names[j]]]
        if (length(reg_hashes_plcs[[curr_names[j]]])==1){
          count = count+1
        }
      }
    }
  }
  
  file_name = paste0("all_reg_hash_plcs_",curr_folder,collapse = "")
  save(all_reg_hash_plcs,file = file_name)
}