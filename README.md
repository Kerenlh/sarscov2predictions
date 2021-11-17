# sarscov2predictions
Statistical modeling of SARS-CoV-2 substitution processes: predicting the next variant

We hereby provide the code used in the paper.
The workflow described here enables to evaluate multiple substitution models for the training data, score them and test their predictions using the test data.

The file main.R contains the entire workflow:

# Input variables
In order to apply our method, the following input variables are required:

"ref_seq.fasta"  - The reference sequence in .fasta format. We used accession NC_045512.2 from the NCBI website as a reference sequence  (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049).

"training_seqs.txt" - Aligned training sequences.
For aligning the sequences one can use the framework developed by Lanfear at https://github.com/roblanf/sarscov2phylo. 
When using this frame work "global.fa" is the desired training sequences alignment file.

"tree.txt" - A phylogenentic tree of the training sequences in Newick format.  
For creating a tree one can use the framework developed by Lanfear at https://github.com/roblanf/sarscov2phylo. When using this framework "ft_SH.tree" is the 
desired training tree.

"site_to_alignment" and "site_to_alignment_testing": A matrix mapping the alignment site numberto the site number at the reference sequence. 
This should be done carefully according to the alignment method. 
I provide here (in the file "main.R") an example for creating this matrix when no sites from the reference sequence are omitted. In case some sites from the reference sequence are omitted this should be done with a tailored code. When using the framework developed by Lanfear for alignment, one should check the "problematic_sites_sarsCov2.mask.vcf" file to find which sites were omitted.

"testing_seqs.aln": Aligned test sequences.

* All input variables above should be saved in a directory called "vars" under the main directory.


# How to cite
Keren Levinstein Hallak, Saharon Rosset. Modeling SARS-CoV-2 substitution processes: predicting the next variant, 11 August 2021, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-654547/v1]
