library(tidyverse)

args <- commandArgs(TRUE)

## read data
t <- read_tsv(args[1], col_names = FALSE) %>%
  mutate(length = nchar(X1))
t

## check how many GBC are shorter than 18 bp
message('Number of GBC which are 18bp long:')
t %>%
mutate(correct_length = length == 18) %>%
pull(correct_length) %>%
table()

## check format
t %>% arrange(length)

## filter for GBC which are 18 bp in length and write to new TSV file
t %>%
filter(length == 18) %>%
select(-length) %>%
write_tsv('GBC_not_corrected_18bp.tsv.gz', col_names = FALSE)
