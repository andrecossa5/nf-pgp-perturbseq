library(tidyverse)

args <- commandArgs(TRUE)

## read in previously created whitelist
whitelist <- read_tsv(args[1], col_names = c('correct','degenerated'))

## bring 'degenerated' column in longer format (one row for every degenerated
## version of every correct GBC)
temp_correct_barcodes <- whitelist %>%
  dplyr::select(correct) %>%
  dplyr::mutate(degenerated = correct) # Correct / correct
temp_degenerated_barcodes <- whitelist %>% 
  dplyr::filter(!is.na(degenerated)) %>%
  dplyr::select(correct,degenerated) %>%
  tidyr::separate_rows(degenerated, sep = ',') # Correct / degenerated
temp_barcodes_final <- bind_rows(temp_correct_barcodes, temp_degenerated_barcodes) %>%
  dplyr::select(degenerated, correct) %>% 
  dplyr::distinct() # Invert degenereated / correct

## check number of unique degenerated (original) GBC and corrected GBC
message(
  paste0(
    'Number of unique degenerated GBC: ', # ALL barcodes found. All of them
    temp_barcodes_final %>%
      dplyr::select(degenerated, correct) %>%
      dplyr::pull(degenerated) %>%
      unique() %>%
      length() %>%
      format(big.mark = ',')
  )
)

message(
  paste0(
    'Number of unique correct GBC: ', # ALL umi-tools correct barcodes. All of them
    temp_barcodes_final %>%
      dplyr::select(degenerated, correct) %>%
      dplyr::pull(correct) %>%
      unique() %>%
      length() %>%
      format(big.mark = ',')
  )
)

## function to calculate Hamming distance (number of mismatches) between two strings
calculateHammingDistance <- function(s1,s2) {
  mapply(function(c1,c2) sum(c1!=c2), strsplit(s1,''), strsplit(s2,''))
}

## Apply 
temp_barcodes_final <- temp_barcodes_final %>%
  mutate(Hamming_distance = calculateHammingDistance(correct, degenerated))

## check how many corrections we can do with only a single mismatch
message(
  paste0(
    'Number of GBC that can be corrected with a single mismatch: ',
    temp_barcodes_final %>%
    filter(Hamming_distance <= 1) %>%
    pull(degenerated) %>%
    unique() %>%
    length() %>%
    format(big.mark = ',')
  )
)

## filter single mismatch corrections and write to TSV file for error correction
temp_barcodes_final %>% 
filter(Hamming_distance <= 1) %>%
select(-Hamming_distance) %>%
write_tsv('GBC_whitelist_proper_format.tsv', col_names = FALSE)
