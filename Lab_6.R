setwd("~/Desktop/Git_Bioinformatics")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

BiocManager::install("msa")

library("biostrings")
library("BiocManager")
library("seqinr")
library("msa")

# Read in each sequence file
seq1 <- readDNAStringSet("seq_1.fasta") 
seq2 <- readDNAStringSet("seq_2.fasta") 
seq3 <- readDNAStringSet("seq_3.fasta") 
seq4 <- readDNAStringSet("seq_4.fasta") 
seq5 <- readDNAStringSet("seq_5.fasta")

# Combine them into one DNAStringSet object
mySequences <- c(seq1, seq2, seq3, seq4, seq5)

# Rename them so you know which is which
names(mySequences) <- c("seq1", "seq2", "seq3", "seq4", "seq5")

#Check it
mySequences

#Run MUSCLE alignment
muscleAlignment <- msa(mySequences, method="Muscle")
print(muscleAlignment, show = "complete")

#Counting the gaps
alignedSeqs <- as.matrix(muscleAlignment)
totalGaps <- sum(alignedSeqs == "-")
totalGaps

#Determining alignment length
ncol(alignedSeqs)

#Determining GC amount
GC(alignedSeqs)

# msa to seqinr
muscleAlignment_seqinr <- msaConvert(muscleAlignment, type = "seqinr::alignment")

#Distance matrix
distMatrix <- dist.alignment(muscleAlignment_seqinr, matrix = "identity")
distMatrix

#Into Amino Acids
seq1_aa <- translate(seq1)
seq1_aa

#Phangorn
install.packages("phangorn")
library("phangorn")
Alignment_phyDat <- msaConvert(muscleAlignment, type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")
write.phyDat
