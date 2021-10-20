library(readr)
library(dplyr)
library(tidyr)

gs_read_mapping <- read_tsv("inputs/CAMI_low/gs_read_mapping.binning.gz", 
                            comment = "@",
                            col_names = c("SEQUENCEID", "BINID", "TAXID", "READID"))

genomes_and_contigs <- gs_read_mapping %>%
  separate(col = READID, into = c("contig", "unk"), sep = "-") %>%
  select(BINID, contig) %>%
  distinct() %>%
  mutate(genome_and_contig = paste0(BINID, "-", contig))

write_tsv(genomes_and_contigs, "inputs/CAMI_low_genomes_and_contigs.tsv")