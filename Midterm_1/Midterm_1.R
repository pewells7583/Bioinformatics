#Set working environment
setwd("/Users/peytonwells/Desktop/Git_Bioinformatics/Midterm_1/")

#Loading in necessary packages

library(BiocManager)
library(Biostrings)
library(msa)
library(seqinr)

#Importing DNA
sequences <- "sequences.fasta"

#Reading sequence file
seqs <- readDNAStringSet(sequences)
print(seqs)

#Aligning with MUSCLE
alignment <- msa(seqs, method="Muscle")
print (alignment)

#Converting to seqinr format for later analysis
seqinr_alignment <- msaConvert(alignment, type = "seqinr::alignment")

#Measuring alignment quality: 
#Pairwise identity matrix
dist_id <- dist.alignment(seqinr_alignment, matrix = "identity")
identity_matrix <- round ((1 - as.matrix(dist_id)^2) * 100, 2)

cat("\nPairwise percent identity matrix (5):\n")
print(identity_matrix )

#Consensus sequence
#consensus matrix and consensus string from biostrings
alignment_biostrings <- as(alignment, "DNAMultipleAlignment")
consensus_seq <- consensusString(alignment_biostrings)
cat("\nConcensus sequence:\n")
cat(consensus_seq, "\n")

#GC Content
aligned_seqs <- as.matrix(alignment)
GC(aligned_seqs)

#Finding the most different individual
print(identity_matrix)
mean_identity <- rowMeans(identity_matrix)
cat("\nMean percent identity to all others:\n")
print(round(sort(mean_identity), 2))

#Homo_Sapiens_6 is the most different individual with a 
#99.08% identity. Homo_sapiens 10 and 4 also have 
#PID lower than 99.93, so there might be 
#substitutions in their sequences

#Exporting sequences for BLAST
writeXStringSet(seqs, "for_blast.fasta")
cat("Exported for_blast.fasta\n")
#Homo sapiens hbb gene for beta globin, 
#Accession number: LC121775.1

#Translating the most different individual's sequence to protein
different_seq <- seqinr_alignment$seq[[6]]
print(different_seq)
different_protein <- translate(s2c(different_seq))
cat("\nProtein sequence for Homo_sapiens_6:\n")
print(different_protein)

#to fasta
different_protein_string <- paste(different_protein, collapse = "")
different_protein_set <- AAStringSet(different_protein_string)
names(different_protein_set) <- "Homo_sapiens_6_protein"
writeXStringSet(different_protein_set, "Homo_sapiens_6_protein.fasta")
cat("Protein FASTA written\n")

# #8- after being put in BlastP, 
# given hemoglobin subunit beta isoform X1 [Theropithecus gelada]
#Sequence ID: XP_025213810.1
#it is associated with sickle cell disease
#and beta thalassemia
#yes, it is likely homo sapien 6 has the disease
