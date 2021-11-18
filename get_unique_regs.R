get_unique_regs = function(iter_num,curr_num,dataset,iterate_vals,output,data_input,data_output,binary_data,
                           save_path){
  reg_hash = hash()
  for (i in ((curr_num-1)*iter_num+1):min((curr_num*iter_num),dim(dataset)[1])){
    # print(i)
    tmp = get_data(dataset[i,],iterate_vals = iterate_vals,output = output,
                   data_input = data_input,data_output = data_output,binary_data = binary_data)
    curr_data = tmp[[1]];
    curr_details = tmp[[2]]
    ys = curr_data[,dim(curr_data)[2]]
    added_log_lik = added_log_lik.P = output_0_or_1_row = 0
    curr_details = cbind(curr_details,added_log_lik,added_log_lik.P, output_0_or_1_row)
    # print(curr_details)
    if (is.na(curr_details[,"AIC"])){
      curr_lik = 1
      for (j in 1:length(ys)){
        ys_j = ys[j]
        curr_lik = curr_lik*dpois(ys_j,ys_j)
        if(curr_lik==0){
          break
        }
      }
      curr_details[,"logLik"] = log(curr_lik)
      curr_details[,"added_log_lik"] = 1
      curr_details[,"df"] = curr_details[,"rows_num"]
    }
    
    #Poisson correction:
    if (is.na(curr_details[,"P.AIC"]) | (curr_details[,"P.converged"]==0 &
                                         (curr_details[,"P.df"]==0 | (curr_details[,"P.df"]!=0 & curr_details[,"P.AIC"]>1e5)))){
      curr_lik = 1
      for (j in 1:length(ys)){
        ys_j = ys[j]
        curr_lik = curr_lik*dpois(ys_j,ys_j)
        if(curr_lik==0){
          break
        }
      }
      curr_details[,"P.logLik"] = log(curr_lik)
      curr_details[,"added_log_lik.P"] = 1
      curr_details[,"P.df"] = curr_details[,"rows_num"] 
    }
    
    #models with 1 row and output!=0
    if (is.na(curr_details[,"NB.converged"]) & is.na(curr_details[,"AIC"])==FALSE 
        & curr_details[,"output_sum"]!=0){
      curr_details[,"P.df"] = 1
      curr_details[,"df"] = 1
      curr_val = curr_details[,"output_sum"]
      curr_lik = exp(-curr_val)*(curr_val^curr_val)/factorial(curr_val)
      if (is.na(curr_details[,"NB.converged"]) & is.na(curr_details[,"AIC"])==FALSE
          & curr_details[,"output_sum"]==curr_val){
        curr_details[,c("logLik","P.logLik")] = log(curr_lik)
        curr_details[,c("added_log_lik","added_log_lik.P")] = 1
      }
      
    }
    if(is.na(curr_details[,"NB.converged"]) & is.na(curr_details[,"AIC"])==FALSE){
      curr_details[,"added_log_lik"] = 1
      curr_details[,"added_log_lik.P"] = 1
    }
    
    if(curr_details[,"rows_num"]==1 | curr_details[,"output_sum"]==0){
      curr_details[,"output_0_or_1_row"] = 1
    }
    # print(i)
    # print(curr_details[(length(iterate_vals)+3):length(curr_details)])
    reg_hash[[dataset$reg_name2[i]]] = 
      curr_details[(length(iterate_vals)+3):length(curr_details)]
    # print(reg_hash[[dataset$reg_name2[i]]])
  }
  
  file_name = paste0("reg_hash.",curr_num,collapse = "")
  save(reg_hash,file = file_name)
}
