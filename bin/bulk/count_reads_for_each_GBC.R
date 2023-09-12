library(tidyverse)

args <- commandArgs(TRUE)

## load corrected GBC
GBC_corrected <- read_tsv(args[1], col_names = FALSE)

## calculate read count by GBC
GBC_count_corrected <- GBC_corrected %>%
  group_by(X1) %>%
  summarize(read_count = n()) %>%
  arrange(-read_count)

## write read counts per GBC to TSV file
GBC_count_corrected %>%
rename(GBC=X1) %>%
write_tsv('read_count_by_GBC_corrected.tsv')
