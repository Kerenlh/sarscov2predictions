# A main directory containing all relevant functions:
begin_path = "/Users/keren/Dropbox/covid/new2/github/" 

# Under the main directory create a directory called vars containing the following files:

# "ref_seq.fasta":  save the reference sequence in the vars directory

# "training_seqs.txt": save all aligned training sequences in the vars directory. 
# For aligning the sequences one can use the framework developed by Lanfear at 
# https://github.com/roblanf/sarscov2phylo . When using this frame work "global.fa" is the 
# desired training sequences alignment.

# "tree": save the phylogenentic training tree in the vars directory. 
# For creating a tree one can use the framework developed by Lanfear at 
# https://github.com/roblanf/sarscov2phylo . When using this framework "ft_SH.tree" is the 
# desired training tree.

# "site_to_alignment" and "site_to_alignment_testing": A matrix mapping the alignment site number
# to the site number at the reference sequence. 
# This should be done carefully according to the alignment method. 
# I provide here an example for creating this matrix when no sites from the reference sequence are
# omitted. In case some sites from the reference sequence are omitted this should be done
# with a tailored code. When using the framework developed by Lanfear for alignment, one should
# check the "problematic_sites_sarsCov2.mask.vcf" file to find which sites were omitted.

# "testing_seqs.aln": save all aligned testing sequences in the vars directory.

setwd(begin_path)
source("COVID_functions.R") # General functions
#####################################################
# Preprocessing:
setwd(paste0(begin_path,"vars/"))
ref_seq = readDNAStringSet(file = "ref_seq.fasta")
seqs = readDNAStringSet(file = "training_seqs.txt")
seqs = short_names_seqs(seqs)
tree = read.tree("tree.txt")
tree = short_names_tree(tree)
save(tree,file = "tree")
# Make sure that the tip labels of the tree are identical to the names(seqs) and if not change 
# the short_names_seqs/short_names_tree functions in "COVID_functions.R".

###################
setwd(begin_path)
source("find_site_patterns.R")
find_site_patterns(begin_path = begin_path, ref_seq = ref_seq, seqs = seqs, tree = tree, testing_flag = 0)
# This is a preprocessing stage for the ASR (Ancestral Sequence Reconstruction) part.
# Goes over all sites to find the unique site patterns and prepares the site patterns in the format 
# required for the next stage. Sites for which all sequences are identical
# do not require ASR. Sites that have identical site patterns should have ASR only once.
# Saves to the ./vars/ folder the files:  "all_no_subs_plcs", "all_site_patterns", "all_site_patterns_plcs",
# "site_patterns_numbers", "seqs_names".
# Takes about an hour
##################################################
# Create "site_to_alignment":
# Note! This code should be used only if no sites are omitted from the reference sequence. 
# If some of the sites were omitted please write a tailored code to generate this matrix as
# explained above. alignment_num is the site number in the alignment and site is the site number
# in the reference sequence.
warning("Improtant!!! Please use this code to generate site_to_alignment only if no sites were
        omitted from the reference sequence. The generated site_to_alignment matrix should always
        be inspected manually.")
create_site_to_aligment = function(num_sites,ref_seq){
  splitted_ref_seq = strsplit(as.character(ref_seq),split="")[[1]]
  ref_seq_plcs = which(splitted_ref_seq!="-")
  site_to_alignment = cbind(1:num_sites,matrix(0,num_sites,1))
  colnames(site_to_alignment) = c("alignment_num","site")
  site_to_alignment = data.frame(site_to_alignment)
  site_to_alignment$site[ref_seq_plcs] = c(1:length(ref_seq_plcs))
  return(site_to_alignment)
}
setwd(paste0(begin_path,"vars/"))
load("num_sites")
site_to_alignment = create_site_to_aligment(num_sites,seqs[["NC_045512"]])
save(site_to_alignment,file = "site_to_alignment")
rm(seqs)
###################
# Creating "tree_list.initial" and "node_depth"
setwd(begin_path)
source("create_initial_tree_and_nodes_depths.R")
create_initial_tree_list(tree,begin_path) #saves "tree_list.initial" to the ./vars/ folder 
setwd(paste0(begin_path,"vars/"))
load("tree_list.initial")
load("seqs_names")
create_node_depth(begin_path, tree, tree_list, seqs_names) # saves "node_depth" to the ./vars/ folder
##################################################

# Ancestral Sequence Reconstruction (ASR)
setwd(paste0(begin_path,"vars/",collapse = ""))
load("site_patterns_numbers")
load("node_depth")
load("all_site_patterns")
unique_site_patterns = site_patterns[which(duplicated(site_patterns)==FALSE)]
dir.create("trees")
dir.create("muts")
setwd(begin_path)
source("ASR.R")

#########################################################
# Important note!!! This part should be done in parallel.
#########################################################
# The ASR for each site takes about 5 minutes. I recommend dividing the unique_site_patterns and
# corresponding site_patterns_numbers into groups according to the available number of cores. 
# The function saves each reconstructed tree as tree_list.(relevant site number) at ./vars/trees.
# The size of a directory containing all tree_list.# files is 6.65GB for the tree reconstructed
# by Lanfear's method described in the paper. I recommend saving all files to a local directory 
# with enough memory and deleting the local copies from all cores. The path to this directory 
# should be inserted as trees_path at the next section.

ASR(begin_path = begin_path,site_patterns_numbers = site_patterns_numbers,node_depth = node_depth,
    unique_site_patterns = unique_site_patterns, tree_list = tree_list)
#########################################################
# Combine all tree_list.# files (that contain ASR for a single site pattern) to reconstruct 
# sequences at all internal nodes of the tree and save them in file called "all_seqs" at the ./vars folder.

setwd(paste0(begin_path,"vars/",collapse = ""))
dir.create("reconstructed_seqs")
load("all_site_patterns_plcs")
load("all_no_subs_plcs")
num_seqs = length(seqs_names)
init_seq_no_subs_plcs = strsplit(as.character(ref_seq),split ="")[[1]][no_subs_plcs]

# trees_path = "/Users/keren/Desktop/covid_files/new2/Lanfear trees/"
trees_path = paste0(begin_path, "vars/trees/")
# This is a path to a directory that contains all the tree_list.(relevant site number) created by
# the ASR function. 

setwd(begin_path)
source("unite_reconstructed_seqs.R")
unite_reconstructed_seqs(begin_path,site_patterns_plcs,no_subs_plcs,num_sites,
                         init_seq_no_subs_plcs)
#################
# Creates "site_details" (saved in ./vars folder) that contains regions - a matrix with details about each site in the 
# reference sequence (codon position, gene,...). Also saves codons_to_amino_acids in the ./vars folder.
setwd(begin_path)
source("site_details.R")
regions = create_site_details(num_sites,site_to_alignment)
setwd(paste0(begin_path,"vars/",collapse = ""))
save(regions,file = "site_details")
###########
# Creates codon_states and codon_outputs - matrices with details about all initial states
# along the tree, their exposure and substitutions.

setwd(paste0(begin_path,"vars/",collapse = ""))
load("all_seqs")
load("site_details")
load("node_depth")
setwd(begin_path)
source("create_codon_states_and_outputs.R")
create_codon_states_and_outputs(tree_list,regions,tree,all_seqs,node_depth,leaves_flag = 0)
create_codon_states_and_outputs(tree_list,regions,tree,all_seqs,node_depth,leaves_flag = 1)
#############
# Summarizes codon_states and codon_outputs into a matrix called codons_table containing all
# information about the different states along the tree and their substitutions (also includes
# site based information like gene, amino acid, amino acid num, details about the reference sequence...)
setwd(paste0(begin_path,"/vars/"))
load("codon_states")
load("codon_outputs")
load("site_details")
setwd(begin_path)
source("codons_table.R")
codons_table = create_codons_table(codon_states,codon_outputs,regions,testing_flag = 0)
setwd(paste0(begin_path,"vars/"))
save(codons_table,file = "codons_table")
#################
# This section produces all datasets for regression and gives them a hash number so that
# regression will be done only once in case of identical datasets (this is common as different
# combinations of partitioning factors and explaining factors produce identical datasets for
# regression).
#########################################################
# Important note!!! This part should be done in parallel.
#########################################################
# It takes a lot of time to go over all possible combinations of inclusion, partition and omission
# of all factors. The for loop below is a convenient way to run this part on multiple cores. 

setwd(begin_path)
source("datasets.R")
setwd(paste0(begin_path,"/vars/"))
load("codons_table")
datasets_preprocess(codons_table)

load("iterate_vals")
load("data_rows")
load("output")
load("selected_input")
load("curr_details")
load("data_input")
load("data_output")
load("binary_data")
load("input_iterate")

dir.create("datasets")

parallel_vals = expand.grid(66,23,1:6,1:5)
# parallel_vals = expand.grid(1:66,1:23,1:6,1:5)
# For each factor the number of possible values is the number of categories + 2:
# The first values refer to partition according to the current categorical value and the last
# 2 options refer to inclusion as an explnantory facto and omission.
# For example, codon.pos (=codon position) gets the values 1-3 for codon positions 1-3,
# 4 is for including the codon position as an explanatory factor and 5 is for omitting the
# the codon position from the regression. 
colnames(parallel_vals) = c("codon","amino_acid","base","codon.pos")
plcs = c(which(parallel_vals[,"codon"]<66 & parallel_vals[,"amino_acid"]!=23),
         which(parallel_vals[,"codon"]<66 & parallel_vals[,"codon.pos"]<5 &
                 parallel_vals[,"base"]<6))
if (length(plcs)>0){
  parallel_vals = parallel_vals[-plcs,]
}

library(hash)
library(digest)
library(MASS)
cat("\n")
pb <- txtProgressBar(min = 0, max = dim(parallel_vals)[1], style = 3, width = 50, char = "=")
start_time = Sys.time()
for (i in 1:dim(parallel_vals)[1]){
  q = parallel_vals[i,"codon"]
  w = parallel_vals[i,"amino_acid"]
  b = parallel_vals[i,"base"]
  j = parallel_vals[i,"codon.pos"]
  datasets(codons_table,q,w,b,j,i,iterate_vals,data_rows,output,selected_input,
           curr_details,data_input,data_output,binary_data,input_iterate)
  setTxtProgressBar(pb, i)
}
print(Sys.time()-start_time)
################################
# Collect datasets to seperate files, 1 for synonymous substitutions and 2 for non-synonymous
# substitutions. Saves the files dataset.#, dataset_reg_name.# and unique_dataset.# to the ./vars folder.
setwd(begin_path)
source("collect_datasets.R")
count1 = get_dataset_size(path = paste0(begin_path,"vars/datasets"),output_num = 1)
collect_datasets(path = paste0(begin_path,"vars/datasets"),total_count = count1,output_num = 1)
count2 = get_dataset_size(path = paste0(begin_path,"vars/datasets"),output_num = 2)
collect_datasets(path = paste0(begin_path,"vars/datasets"),total_count = count2,output_num = 2)

##############
# Run regressions for all unique datasets. The sysnonymous and non-synonymous regressions are
# done separately (hence the separation for unique dataset.1/2. This is more convenient but
# not mandatory).
#########################################################
# Important note!!! This part should be done in parallel.
#########################################################
# I recommend to go over the for loops in parallel when multiple cores are available.
setwd(begin_path)
source("regression_functions.R")
source("get_unique_regs.R")
library(hash)
library(digest)
setwd(paste0(begin_path,"vars/"))
load("codons_table")
load("iterate_vals")
load("output")
load("data_input")
load("data_output")
load("binary_data")
load("unique_dataset.1")
iter_num = 500
rel_nums = 1:ceiling(dim(dataset)[1]/iter_num)
dir.create("unique_regs")
setwd(paste0(begin_path,"vars/unique_regs/"))
dir.create("dataset_1")
dir.create("dataset_2")
setwd(paste0(begin_path,"vars/unique_regs/dataset_1"))
cat("\n")
pb <- txtProgressBar(min = 0, max = length(rel_nums), style = 3, width = 50, char = "=")
for (curr_num in rel_nums){
  get_unique_regs(iter_num,curr_num,dataset,iterate_vals,output,data_input,data_output,
                  binary_data)
  setTxtProgressBar(pb, curr_num)
}
rm(dataset)
setwd(paste0(begin_path,"vars/"))
load("unique_dataset.2")
rel_nums = 1:ceiling(dim(dataset)[1]/iter_num)
setwd(paste0(begin_path,"vars/unique_regs/dataset_2"))
cat("\n")
pb <- txtProgressBar(min = 0, max = length(rel_nums), style = 3, width = 50, char = "=")
for (curr_num in rel_nums){
  get_unique_regs(iter_num,curr_num,dataset,iterate_vals,output,data_input,data_output,
                  binary_data)
  setTxtProgressBar(pb, curr_num)
}
rm(dataset)
#################
# Creates all_reg_hash_plcs_dataset_1/2. This is a list that maps regression names to the
# relevant rows in dataset_1/2.
iter_num = 500
source(paste0(begin_path,"find_reg_hash_plcs.R"))
source(paste0(begin_path,"collect_all_reg_hash_plcs.R"))
setwd(paste0(begin_path,"vars/"))
load("dataset_reg_name.1")
dir.create("reg_hash_plcs")
setwd(paste0(begin_path,"vars/reg_hash_plcs/"))
dir.create("dataset_1")
dir.create("dataset_2")
setwd(paste0(begin_path,"vars/reg_hash_plcs/dataset_1/"))
reg_hash_names = unique(dataset_reg_name)
total_num_files = round(length(reg_hash_names)/iter_num)+1
cat("\n")
pb <- txtProgressBar(min = 0, max = total_num_files, style = 3, width = 50, char = "=")
for (k in 1:total_num_files){
  find_reg_hash_plcs(dataset_reg_name,iter_num,reg_hash_num = k)
  setTxtProgressBar(pb, k)
}
rm(dataset_reg_name)
collect_all_reg_hash_plcs(begin_path,curr_folder = "dataset_1")

setwd(paste0(begin_path,"vars/"))
load("dataset_reg_name.2")
setwd(paste0(begin_path,"vars/reg_hash_plcs/dataset_2/"))
reg_hash_names = unique(dataset_reg_name)
total_num_files = round(length(reg_hash_names)/iter_num)+1
cat("\n")
pb <- txtProgressBar(min = 0, max = total_num_files, style = 3, width = 50, char = "=")
for (k in 1:total_num_files){
  find_reg_hash_plcs(dataset_reg_name,iter_num,reg_hash_num = k)
  setTxtProgressBar(pb, k)
}
collect_all_reg_hash_plcs(begin_path,curr_folder = "dataset_2")

##########################
# Creates details.1/2 and saves it in the ./vars directory. details.1/2 contain the regression
# results for all datasets.

source(paste0(begin_path,"find_details_from_datasets.R"))
find_details_from_datasets(curr_folder = "dataset_1",folder_num = 1)
find_details_from_datasets(curr_folder = "dataset_2",folder_num = 2)

###############
# Each model is composed of several datasets and their regression results that appear in the details
# file. In this section, the results of the sub-models composing each model are combined to 
# create a unified score for the model. All models are then ranked by the minimum of their 
# Poisson and Negative Binomial (NB) AIC score and saved as "unique_models" at ./vars.
# "model_ids" is also produced and saved here (necessary later on in "predict_substitutions.R").
source(paste0(begin_path,"score_models.R"))

###############################################################################
#                                   Testing:                                  #
###############################################################################

setwd(begin_path)
dir.create("testing")
setwd(paste0(begin_path,"testing"))
dir.create("vars")
setwd(paste0(begin_path,"vars/"))
load("seqs_names")
new_seqs = readDNAStringSet(file = "testing_seqs.aln")
new_seqs = short_names_seqs(new_seqs)
new_seqs_ref_seq = new_seqs[["NC_045512"]]
new_seqs_plcs = which(is.element(names(new_seqs),seqs_names)==FALSE)[-1]
new_seqs = new_seqs[new_seqs_plcs]
save(new_seqs,file = "new_seqs_for_testing")
source(paste0(begin_path,"find_site_patterns.R"))
find_site_patterns(begin_path = paste0(begin_path,"testing/"), ref_seq = ref_seq,
                   seqs = new_seqs, tree = NULL, testing_flag = 1)
#########################################
# Create "site_to_alignment_testing":
# Note! This code should be used only if no sites are omitted from the reference sequence. 
# If some of the sites were omitted please write a tailored code to generate this matrix as
# explained above. alignment_num is the site number in the alignment and site is the site number
# in the reference sequence.
warning("Improtant!!! Please use this code to generate site_to_alignment only if no sites were
        omitted from the reference sequence. The generated site_to_alignment matrix should always
        be inspected manually.")
setwd(paste0(begin_path,"testing/vars/"))
load("num_sites")
site_to_alignment_testing = create_site_to_aligment(num_sites,new_seqs_ref_seq)
save(site_to_alignment_testing,file = "site_to_alignment_testing")
#########################################
# The function get_testing_subs in find_site_patterns.R creates "new_subs_testing" and 
# "rel_ref_sites_no_subs" and saves them at ./testing/vars/
setwd(paste0(begin_path,"testing/vars/"))
load("all_no_subs_plcs")
load("all_site_patterns_plcs")
load("all_site_patterns")
no_subs_plcs_testing = no_subs_plcs; site_patterns_plcs_testing = site_patterns_plcs
site_patterns_testing = site_patterns
setwd(paste0(begin_path,"vars/"))
load("all_no_subs_plcs")
load("all_site_patterns_plcs")
load("all_site_patterns")
load("site_details")
load("codons_table")
load("site_to_alignment")
get_testing_subs(begin_path,no_subs_plcs,site_patterns_plcs,site_patterns,site_to_alignment,
                 new_seqs,codons_table,regions,
                 site_to_alignment_testing,no_subs_plcs_testing,
                 site_patterns_plcs_testing,site_patterns_testing)
rm(new_seqs)

#############################
source(paste0(begin_path,"codons_table.R"))
setwd(paste0(begin_path,"testing/vars/"))
load("new_subs_testing")
load("rel_ref_sites_no_subs")
setwd(paste0(begin_path,"vars/"))
load("codon_states_leaves")
load("codon_outputs_leaves")
load("syn_tv_const")
load("non_syn_tv_const")
tmp = create_codons_table(codon_states = codon_states,
                          codon_outputs = codon_outputs,regions = regions,
                          testing_flag = 1)
codons_table_testing = create_codons_table_testing(tmp,new_subs_testing,rel_ref_sites_no_subs,
                            syn_tv_const,non_syn_tv_const)
save(codons_table_testing,file = "codons_table_testing")

##############################
# Predict substitutions for the testing data based on top models trained on the training data.
# num_first_models is the number of top models evaluated.
source(paste0(begin_path,"predict_substitutions.R"))
setwd(paste0(begin_path,"vars/"))
load("codons_table")
load("codons_table_testing")
load("iterate_vals")
load("unique_models")
load("details.12")
load("output")
load("data_input")
load("data_output")
load("binary_data")
load("model_ids")
predict_substitutions(begin_path,codons_table_training = codons_table,codons_table_testing,
                      iterate_vals,unique_models,details,model_ids,num_first_models=3,output,data_input,
                      data_output,binary_data)

#######################################################
# Get AUC results and lift plot for the examined model.
source(paste0(begin_path,"analyse_predictions.R"))
examined_model_num = 3 # The number is the rank of the model in the sorted unique_models list.
model_type = "p.non_syn" # The options are: "nb.syn","p.syn","nb.non_syn","p.non_syn"
# p. = Poisson, nb. = Negative Binomial, syn = synonymous, non_syn = non-synoynymous
spike_flag = 0 # 0 gives results for all genes and 1 gives results only for the spike gene.
tmp = analyse_predcitions(begin_path,spike_flag,examined_model_num,model_type)
results = tmp[[1]]; lift_plot = tmp[[2]]                    
