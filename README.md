# New version, bulk_is_the_key

A Nextflow pipeline for lentiviral-barcoding (Dixit et al., 2016) + scRNA-seq experiments data pre-processing.

Major steps:

BULK

- Find GBC reads (<=1 mismatch from DNA anchor sequence)
- UMI-tools correction. All read counts from degenerate sequences are passed to the "correct" ones.
- (Optional) clone calling step from bulk DNA. Can use spikeis information or Gaussian Mixture Models to do this.
  MOI of lentiviral infection must be <=.1 to retrieve clones from bulk DNA-sequencing, to ensure one GBC --> one clone. 

SC

- Solo alignment of 10x reads
- Extract GBC reads element from lentiviral-barcoding library (<= 2 mismatch from RNA anchor sequence), using only reads from putative cells (Solo CBs)
- Bulk-reference-based GBC filtering, correction and CB-GBC UMI table computation. SC CB-UMI-GBC combinations with GBC matching 
  a unique "correct" bulk GBC (1-2 max mismatches allowed) are filtered and (if necessary) corrected. Outliers CB 
  (too high/low GBCs counts) are removed. Highly supported CB-GBC combinations (>=30 reads, >=5 UMIs and >=10 coverage) are filtered, 
  and used to compute a CB-GBC UMI table. 
- Clone calling. Unique GBCs-combinations in the CB-GBC UMI table are computed, and individual CBs are assigned to them. GBCs-combinations with one or
  more of their GBCs shared across other GBCs-combos are filtered out (un-realistic). The remaining GBCs-combinations are considered clones.
- Cell assignment. CB assigned to the final pool of GBCs-combinations are assigned the respective clonal identity. Clonal frequencies are computed
  after this step, using assigned cell counts.     
