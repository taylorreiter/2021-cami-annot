library(dplyr)
library(readr)
library(purrr)
library(tidyr)

files <- unlist(snakemake@input[['tsv']])
#files = Sys.glob("outputs/gs_read_annotations/*tsv")

gs_read_annots <- files %>%
  set_names() %>%
  map_dfr(read_tsv, col_names = c("scaffold", "source", "feature", "start", "end",
                                  "annot_score", "strand", "frame", "attribute", "scaffold2", 
                                  "read_start", "read_end", "read_name", "map_score", 
                                  "read_strand"), .id = "source_genome") %>%
  mutate(source_genome = gsub("outputs\\/gs_read_annotations\\/", "", source_genome)) %>%
  separate(source_genome, into = c("source_genome", "scaffold3"), sep = "-") %>%
  select(-scaffold2, -scaffold3) %>%
  filter(feature != "region") %>%
  separate(attribute, into = "id", sep = ";", remove = F) %>%
  mutate(id = gsub("ID=", "", id))

write_tsv(gs_read_annots, snakemake@output[['tsv']])
