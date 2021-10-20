library(readr)
library(dplyr)
library(tidyr)

gs_read_mapping <- read_tsv("inputs/CAMI_low/gs_read_mapping.binning.gz", 
                            #n_max = 10000,
                            comment = "@",
                            col_names = c("SEQUENCEID", "BINID", "TAXID", "READID"))

# add paired end information. Arrange so read names are in interleaved format.
gs_read_mapping_r1 <- gs_read_mapping %>%
  mutate(SEQUENCEID = paste0(SEQUENCEID, "/1"))
gs_read_mapping_r2 <- gs_read_mapping %>%
  mutate(SEQUENCEID = paste0(SEQUENCEID, "/2"))
gs_read_mapping_r1_r2 <- bind_rows(gs_read_mapping_r1, gs_read_mapping_r2) %>%
  arrange(SEQUENCEID)

# write all reads that correspond to a single genome:contig to a file
gs_read_mapping_r1_r2 %>%
  separate(col = READID, into = c("contig", "unk"), sep = "-") %>%
  select(-TAXID, -unk) %>%
  group_by(BINID, contig) %>%
  group_walk(~ write_tsv(., paste0("outputs/test/", .y$BINID, "-", .y$contig, ".txt"), col_names = F))
