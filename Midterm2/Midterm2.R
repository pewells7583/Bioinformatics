library("BiocManager")
library("ape")
#reading in the file
gene_alignment <- read.dna("metazoa_alignment.gene.fasta", format = "fasta")
#looking at the names in the file
rownames(as.character(gene_alignment))
#translating H. sapiens to protein
sapiens_alignment <- as.character(gene_alignment)["Homo_sapiens", ]
sapiens_nogaps <- sapiens_alignment[sapiens_alignment != "-"]
sapiens_dna <- as.DNAbin(matrix(sapiens_nogaps, nrow = 1,
                           dimnames = list("Homo_sapiens", NULL)))

sapiens_protein <- trans(sapiens_dna, codonstart = 1)

#writing as .fasta file
write.dna(sapiens_protein,
          file   = "homo_sapiens_protein.fasta",
          format = "fasta")
#copy for blast
as.character(sapiens_protein)
protein_chars <- as.character(sapiens_protein)
cat(paste(protein_chars, collapse=""))

#uniprot/GO load ins
library("UniprotR")

#getting GO terms
polg_data <- GetProteinGOInfo("P54098")

#viewing it
print(polg_data)

#plotting GO info and saving as pdf
pdf("POLG_GO_terms_plot.pdf", width = 8, height = 6)


