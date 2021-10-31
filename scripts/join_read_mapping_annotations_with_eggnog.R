library(dplyr)
library(purrr)
library(readr)

files <- unlist(snakemake@input[['eggnog']])
#files = Sys.glob("outputs/eggnog_source_genomes/*.annotations")

eggnog <- files %>%
  set_names() %>%
  map_dfr(read_tsv, col_names = c('id', 'seed_ortholog', 'evalue', 'score', 
                                  'eggNOG_OGs', 'max_annot_lvl', 'COG_category', 
                                  'Description', 'Preferred_name', 'GOs', 'EC',
                                  'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                  'KEGG_Reaction', 'KEGG_rclass', 'BRITE',
                                  'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'),
          comment = "#", na = "-", .id = "source_genome") %>%
  mutate(source_genome = gsub("outputs/eggnog_source_genomes/", "", source_genome)) %>%
  mutate(source_genome = gsub("\\.emapper\\.annotations", "", source_genome))

gs_read_annot <- read_tsv(snakemake@input[['gs_read_annot']])
full_gs_read_annots <- left_join(gs_read_annot, eggnog, by = c("source_genome", "id"))

write_tsv(full_gs_read_annots, snakemake@output[['tsv']])
