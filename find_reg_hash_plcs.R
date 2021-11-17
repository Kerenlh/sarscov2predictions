find_reg_hash_plcs = function(dataset_reg_name,iter_num,reg_hash_num){
  reg_hashes_plcs = hash()
  for (j in ((reg_hash_num-1)*iter_num+1):min((reg_hash_num*iter_num),length(reg_hash_names))){
    # print(j)
    plcs = which(dataset_reg_name==reg_hash_names[j])
    reg_hashes_plcs[[ reg_hash_names[j] ]] = plcs
  }
  file_name = paste0("reg_hashes_plcs.",reg_hash_num,collapse = "")
  save(reg_hashes_plcs,file = file_name)
}

