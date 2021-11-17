get_dataset_size = function(path,output_num){
  setwd(path)
  file_names = dir()
  count = 0
  count2 = 0
  start_time = Sys.time()
  cat("\n")
  pb <- txtProgressBar(min = 0, max = length(file_names), style = 3, width = 50, char = "=")
  for (i in 1:length(file_names)){
    # print(i)
    if (paste0(strsplit(file_names[i],split = "")[[1]][1:14]
               ,collapse = "")=="short_details."){
      load(file_names[i])
      count2 = count2+1
      # print(count2)
      if (is.null(short_details)==FALSE){
        output_plcs = which(short_details$output==output_num)
        if (length(output_plcs)>0){
          rel_short_details = short_details[output_plcs,]
          # dataset[(count+1):(count+dim(rel_short_details)[1]),] = rel_short_details
          count = count + dim(rel_short_details)[1]
          rm(short_details)
          rm(rel_short_details)
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  return(count)
}

collect_datasets = function(path,output_num,total_count){
  setwd(path)
  file_names = dir()
  x = file.info(file_names)
  biggest_file_name = row.names(x[order(x$size),])[length(file_names)]
  load(biggest_file_name)
  dataset = data.frame(matrix("-",total_count, length(short_details)))
  colnames(dataset) = colnames(short_details)
  rm(short_details)
  count = 0
  count2 = 0
  start_time = Sys.time()
  cat("\n")
  pb <- txtProgressBar(min = 0, max = length(file_names), style = 3, width = 50, char = "=")
  for (i in 1:length(file_names)){
    # print(i)
    if (paste0(strsplit(file_names[i],split = "")[[1]][1:14]
               ,collapse = "")=="short_details."){
      load(file_names[i])
      count2 = count2+1
      # print(count2)
      if (is.null(short_details)==FALSE){
        output_plcs = which(short_details$output==output_num)
        if (length(output_plcs)>0){
          rel_short_details = short_details[output_plcs,]
          dataset[(count+1):(count+dim(rel_short_details)[1]),] = rel_short_details
          count = count + dim(rel_short_details)[1]
          rm(short_details)
          rm(rel_short_details)
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  print(Sys.time()-start_time)
  # print(dim(dataset))
  setwd(paste0(begin_path,"vars/"))
  save(dataset,file = paste0("dataset.",output_num))
  dataset_reg_name = dataset$reg_name2
  save(dataset_reg_name, file = paste0("dataset_reg_name.",output_num))
  plcs = which(duplicated(dataset$reg_name2)==FALSE)
  dataset = dataset[plcs,]
  save(dataset,file = paste0("unique_dataset.",output_num))
}