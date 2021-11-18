# sarscov2predictions
"Statistical modeling of SARS-CoV-2 substitution processes: predicting the next variant"

This is the code used in the paper https://www.researchsquare.com/article/rs-654547/v1.

We describe a workflow for evaluating multiple substitution models, scoring them, and testing their predictions against the test data.

The file "main.R" contains detailed comments and the entire workflow:
1. A preprocessing section is followed by Ancestral Sequence Reconstruction using a heuristic we developed, inspired by Fitch's algorithm. 
2. We create tables of observed states (defined as site + codon + left and right neighboring nucleotides) along with their corresponding substitutions and exposure. 
3. Information from these tables ("codon_states" and "codon_outputs") is combined with site-based information (appears in "site_details") to create "codons_table" - a matrix containing all states, their substitutions and site-based information.
4. Modeling: All possible models (all possible combinations of inclusion as an explanatory factor, partition, and omission of all factors mentioned in the paper) are considered. Each model is composed of sub-models requiring subsets of the codons_table. The function "datasets.R" goes over all subsets and assigns them a hash number (many different sub-models produce identical subsets of the "codons_table", assigning a hash number allows the time-consuming regression part to be done once). Then, regression is performed for all unique subsets. 
5. The regression results for all sub-models are combined to find the overall AIC score of each model, and models are then ranked by their AIC score (minimum of the Poisson and Negative Binomial AIC scores). All ranked models are saved as "unique_models".
6. Testing: Substitutions are inferred only for sites with no substitutions along the phylogenetic training tree, that their test sequences had at most one base different from the reference sequence base. A table containing all test states, their inferred substitutions, and site-based information is created ("codons_table_testing"). The top 3 models (the number can be changed by changing "num_first_models" in the "main.R" file) evaluated from the training data are used to predict substitutions for the leaves of the training phylogenetic tree. The results are then compared to the substitutions inferred from the test sequences. The code produces AUC and 3% Lift results for the top models and for the two null models mentioned in the paper, as well as a Lift plot for a selected model (saved as "results" and "lift_plot", please see the comments in the "main.R" file to control the plot).


# Input variables
In order to apply our method, the following input variables are required:

1. "ref_seq.fasta" - The reference sequence in .fasta format. We used accession NC_045512.2 from the NCBI website as a reference sequence (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049).

2. "training_seqs.txt" - Aligned training sequences. For aligning the sequences one can use the framework developed by Lanfear at https://github.com/roblanf/sarscov2phylo. When using this frame work "global.fa" is the desired training sequences alignment file.

3. "tree.txt" - A phylogenetic tree of the training sequences in Newick format.
For creating a tree one can use the framework developed by Lanfear at https://github.com/roblanf/sarscov2phylo. When using this framework "ft_SH.tree" is the desired training tree.

4. "site_to_alignment" and "site_to_alignment_testing": A matrix mapping the alignment site number to the site number at the reference sequence. This should be done carefully according to the alignment method.

We provide here (in the file "main.R") an example for creating this matrix when no sites from the reference sequence are omitted. In case some sites from the reference sequence are omitted this should be done with a tailored code. When using the framework developed by Lanfear for alignment, one should check the "problematic_sites_sarsCov2.mask.vcf" file to find which sites were omitted.

5. "testing_seqs.aln": Aligned test sequences.

* All input variables above should be saved in a directory called "vars" under the main directory.
* We provide a very small example of training and test sequences (200 sequences each) along with a phylogenetic tree inferred for the training sequences, downloaded from the NCBI website. https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049
* For this small example, "site_to_alignment" and "site_to_alignment_testing" are created automatically.
* Even though the example contains only 200 training and testing sequences, the running time for the modeling part is very long since all models are considered. To make the running time reasonable it is possible to consider only a small part of the models by uncommenting line 194 in the "main.R" code (making parallel_vals = expand.grid(66,23,6,5)).



# Getting started
In order to use our workflow, clone the repository (git clone https://github.com/Kerenlh/sarscov2predictions.git).
Replace the input variables in the ./vars directory with your desired reference sequence, aligned training and test sequences and a phylogenetic tree for the training sequecnes. Also create "site_to_alignment" and "site_to_alignment_testing" as explained above and add them to the ./vars directory.
Next, go to the main directory and source the "main.R" file. 
For higher efficiency, some parts of the code are needed to be done in parallel using multiple cores. Please see the next section for details. 

# Running time note
Some parts of the "main.R" code take a very long time and should be done in parallel. We wrote some of them as the inner part of for loops in the "main.R" code so that it will be easier to run them on multiple cores.
The time consuming functions are:
"ASR.R" (Ancestral sequence reconstruction per site).
The function "datasets" in "datasets.R".
The function "get_unique_regs" in "get_unique_regs.R".

# How to cite
Keren Levinstein Hallak, Saharon Rosset. Modeling SARS-CoV-2 substitution processes: predicting the next variant, 11 August 2021, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-654547/v1]
