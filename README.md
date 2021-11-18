# sarscov2predictions
Statistical modeling of SARS-CoV-2 substitution processes: predicting the next variant

We hereby provide the code used in the paper https://www.researchsquare.com/article/rs-654547/v1.
The workflow described here enables to evaluate multiple substitution models for the training data, score them and test their predictions using the test data.

The file main.R contains the entire workflow:
A preprocessing section is followed by Ancestral Sequence Reconstruction using a heuristic we developed, inspired by Fitch's algorithm. 

# Input variables
In order to apply our method, the following input variables are required:

"ref_seq.fasta"  - The reference sequence in .fasta format. We used accession NC_045512.2 from the NCBI website as a reference sequence  (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049).

"training_seqs.txt" - Aligned training sequences.
For aligning the sequences one can use the framework developed by Lanfear at https://github.com/roblanf/sarscov2phylo. 
When using this frame work "global.fa" is the desired training sequences alignment file.

"tree.txt" - A phylogenentic tree of the training sequences in Newick format.  
For creating a tree one can use the framework developed by Lanfear at https://github.com/roblanf/sarscov2phylo. When using this framework "ft_SH.tree" is the 
desired training tree.

"site_to_alignment" and "site_to_alignment_testing": A matrix mapping the alignment site number to the site number at the reference sequence. 
This should be done carefully according to the alignment method. 
I provide here (in the file "main.R") an example for creating this matrix when no sites from the reference sequence are omitted. In case some sites from the reference sequence are omitted this should be done with a tailored code. When using the framework developed by Lanfear for alignment, one should check the "problematic_sites_sarsCov2.mask.vcf" file to find which sites were omitted.

"testing_seqs.aln": Aligned test sequences.

* All input variables above should be saved in a directory called "vars" under the main directory. 
* We provide a very small example of training and test sequences (200 sequences each) along with a phylogenetic tree inferred for the training sequences, downloaded from the NCBI website. https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049
* For this small example, "site_to_alignment" and "site_to_alignment_testing" are created automatically.
* Even though the example contains only 200 training and testing sequences, the running time for the modeling part is very long since all models are considered. To make the running time reasonable it is possible to consider only a small part of the models by uncommenting line 194 in the "main.R" code (making parallel_vals = expand.grid(66,23,6,5)). 

# Getting started
In order to use our workflow, please clone the repository (git clone https://github.com/Kerenlh/sarscov2predictions.git).
Replace the input variables in the ./vars directory with your desired reference sequence, aligned training and test sequences and a phylogenetic tree for the training sequecnes. Also create "site_to_alignment" and "site_to_alignment_testing" as explained above and add them to the ./vars directory.

# Important note!
Some parts of the "main.R" code take a very long time and should be done in parallel. I wrote some of them as the inner part of for loops in the "main.R" code so that it will be easier to run them on multiple cores.
The time consuming functions are:
"ASR.R" (Ancestral sequence reconstruction per site).
The function "datasets" in "datasets.R".
The function "get_unique_regs" in "get_unique_regs.R".




# How to cite
Keren Levinstein Hallak, Saharon Rosset. Modeling SARS-CoV-2 substitution processes: predicting the next variant, 11 August 2021, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-654547/v1]
