raxml-ng --msa metazoa_alignment.5k.fasta --model GTR+G --prefix metazoa_5k_out
raxml-ng --msa metazoa_alignment.gene.fasta --model GTR+G --prefix metazoa_gene_out 
raxml-ng --all --msa metazoa_alignment.5k.fasta --model GTR --prefix metazoa_5k_out --bs-trees 100 --redo_lh_impr

