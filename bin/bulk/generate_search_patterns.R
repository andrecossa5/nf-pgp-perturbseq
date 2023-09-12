library(tidyverse)

args <- commandArgs(TRUE)

sequence <- args[1]

array <- tibble(pattern = rep(NA, nchar(sequence)))

for ( i in seq(nchar(sequence)) ) {
  temp_sequence <- strsplit(sequence, split = '')[[1]]
  temp_sequence[i] <- '.'
  temp_sequence <- paste0(temp_sequence, collapse = '')
  array$pattern[i] <- temp_sequence
}

array %>%
mutate(pattern = paste0(pattern, '.*')) %>%
write_tsv('search_patterns.tsv', col_names = FALSE)
