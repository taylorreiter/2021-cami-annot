library(rtracklayer)

gff <- import(snakemake@input[['gff']])
export(gff, snakemake@output[['no_fasta']])
