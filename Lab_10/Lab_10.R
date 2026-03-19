#install packages
library(BiocManager)
BiocManager::install("GenomicAlignments")
install.packages("UniprotR")
library(UniprotR)
install.packages("protti")
library(protti)
library(devtools)
devtools::install_github("swsoyee/r3dmol") 
library("r3mol")
library("Biostrings")

#I had to re-make the protein file. read in the sequence 
#file, translated it to protein, made the new .fasta
read.csv("seq_1.fasta")
seq <- readDNAStringSet("seq_1.fasta")
protein_seq <- translate(seq)
print(protein_seq)
writeXStringSet(protein_seq, file="seq1_protein.fasta")
#ran a blast
#made the txt file, now reading in 
accessions <- read.csv("accession_matches.txt", 
                       header = FALSE)
acc_string <- paste(accessions$V1)
print(acc_string)
#getting GO terms
GO_terms <- GetProteinGOInfo(acc_string)
print(GO_terms)
PlotGoInfo((GO_terms))
#the format is weird
PlotGOAll(
  GOObj = GO_terms,
  Top = 10,
  directorypath = "~/Desktop/Git_Bioinformatics/Lab_10",  
  width = 8, 
  height = 5
)
#that worked!
#interesting GO terms respiratory chain comlplex 1 I suppose  
#finding stuff now
Pathology_object <- GetPathology_Biotech(acc_string)
Get.diseases(Pathology_object)
#NULL pops up
# #10 and 11
fetch_uniprot(acc_string)
fetch_pdb(acc_string)
fetch_pdb("1ZMR")
fetch_pdb("2HWG")
fetch_alphafold_prediction(acc_string) 
GetpdbStructure(acc_string)
#I did it through the website

