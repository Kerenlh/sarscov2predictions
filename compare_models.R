begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
# begin_path = "C:/Users/ifog/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/"))
load("codons_table")
ti_tv_correction = sum(codons_table$y)*log(2)-sum(lchoose(codons_table$y,codons_table$transitions))
syn_non_syn_correction = sum(codons_table$y)*log(2)-sum(lchoose(codons_table$y,codons_table$syn))
ACGT_correction = sum(codons_table$y)*log(2)-sum(lfactorial(codons_table$y))+sum(lfactorial(codons_table$A))+
  sum(lfactorial(codons_table$C))+sum(lfactorial(codons_table$G))+sum(lfactorial(codons_table$T))

add_likelihood_correction = function(unique_models,correction){
  table(2*(unique_models$df.NB -unique_models$logLik.NB)==unique_models$AIC.NB)
  unique_models$logLik.NB = unique_models$logLik.NB+correction
  unique_models$logLik.P = unique_models$logLik.P+correction
  unique_models$AIC.NB = 2*(unique_models$df.NB -unique_models$logLik.NB)
  unique_models$AIC.P = 2*(unique_models$df.P -unique_models$logLik.P)
  return(unique_models)
}

setwd(paste0(begin_path,"vars/models/"))
load("unique_models.12")
unique_models = add_likelihood_correction(unique_models = unique_models,correction = 
                                            syn_non_syn_correction)
save(unique_models,file = "unique_models_corrected.12")
load("unique_models.34")
unique_models = add_likelihood_correction(unique_models = unique_models,correction = 
                                            ti_tv_correction)
save(unique_models,file = "unique_models_corrected.34")
load("unique_models.6789")
unique_models = add_likelihood_correction(unique_models = unique_models,correction = 
                                            ACGT_correction)
save(unique_models,file = "unique_models_corrected.6789")

setwd(paste0(begin_path,"vars/models/"))
load("unique_models_corrected.12")
all_unique_models = unique_models
load("unique_models_corrected.34")
all_unique_models = rbind(all_unique_models,unique_models)
load("unique_models.5")
all_unique_models = rbind(all_unique_models,unique_models)
load("unique_models_corrected.6789")
all_unique_models = rbind(all_unique_models,unique_models)
unique_models = all_unique_models
rm(all_unique_models)
unique_models = unique_models[order(apply(unique_models[,c("AIC.NB","AIC.P")],1,min)),]
unique_models[1:20,]
save(unique_models,file = "all_unique_models_corrected")


# setwd(paste0(begin_path,"vars/models/"))
# load("unique_models_corrected.12")
# all_unique_models = unique_models[which(unique_models$base==3),]
# load("unique_models_corrected.34")
# all_unique_models = rbind(all_unique_models,unique_models[which(unique_models$base==3),])
# load("unique_models.5")
# all_unique_models = rbind(all_unique_models,unique_models[which(unique_models$base==3),])
# load("unique_models_corrected.6789")
# all_unique_models = rbind(all_unique_models,unique_models[which(unique_models$base==3),])
# unique_models = all_unique_models
# rm(all_unique_models)
# save(unique_models,file = "all_unique_models_corrected_no_base")

# unique_models = unique_models[order(apply(unique_models[,c("AIC.NB","AIC.P")],1,min)),]
# unique_models[1:20,]

# begin_path = "/Users/keren/Dropbox/covid/new2/"
# setwd(paste0(begin_path,"vars/models_old/"))
# load("all_unique_models_corrected_no_base")
# unique_models = unique_models[order(apply(unique_models[,c("AIC.NB","AIC.P")],1,min)),]
# unique_models[1:20,]


begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/models/"))
load("all_unique_models_corrected")
unique_models_Lanfear = unique_models
begin_path = "/Users/keren/Dropbox/covid/new2/ncbi_tree/"
setwd(paste0(begin_path,"vars/models/"))
load("all_unique_models_corrected")
unique_models_ncbi = unique_models

score = c(1:dim(unique_models)[1])
unique_models_Lanfear = cbind(unique_models_Lanfear,score)
unique_models_ncbi = cbind(unique_models_ncbi,score)

names(score) = unique_models_ncbi$ID
ncbi_score = score[as.character(unique_models_Lanfear$ID)]
unique_models_Lanfear = cbind(unique_models_Lanfear,ncbi_score)

names(score) = unique_models_Lanfear$ID
Lanfear_score = score[as.character(unique_models_ncbi$ID)]
unique_models_ncbi = cbind(unique_models_ncbi,Lanfear_score)

which.min(unique_models_Lanfear$score+unique_models_Lanfear$ncbi_score)

begin_path = "/Users/keren/Dropbox/covid/new2/"
setwd(paste0(begin_path,"vars/models/"))
save(unique_models_Lanfear,file = "Lanfear_models")
save(unique_models_ncbi,file = "ncbi_models")

setwd(paste0(begin_path,"vars/models/"))
load("Lanfear_models")
load("ncbi_models")

unique_models_Lanfear[order(unique_models_Lanfear$AIC.P),]
unique_models_ncbi[order(unique_models_ncbi$AIC.P),]