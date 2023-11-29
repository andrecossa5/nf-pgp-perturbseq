# Beta version: bulk_is_the_key

A Nextflow pipeline for lentiviral-barcoding single-cell lineage tracing (scLT) experiments data pre-processing.

Major steps:

BULK DNA libraries

- Extract GBC reads (<=1 mismatch from DNA anchor sequence) observed >1 read.
- Graph-based correction. Unique GBC sequences are clustered and read counts from degenerate
  sequences are passed to the "correct" ones.
- Reference GBC pool construction. A reference GBC whitelist is obtained from all samples.
- (Optional) clone calling step from bulk DNA. Two methods implemented: a. uses spike-ins
  information, and b. does not, relying on distributional properties of the GBC-sequences
  read counts only. 

Note: MOI of lentiviral infection must be very low (e.g., <=.1) to retrieve clones from bulk DNA-sequencing, to ensure one GBC --> one clone. 

SC libraries (10X and GBC reads)

- Solo alignment of 10x reads
- GBC reads element extraction (<= 2 mismatch from RNA anchor sequence), using only reads from
  putative cells (Solo-corrected CBs).
- Bulk-reference-based GBC filtering, correction and CB-GBC UMI table computation. SC CB-UMI-GBC
  combinations with GBC matching a unique, "correct" bulk GBC (1-2 max mismatches allowed) are filtered and (if necessary) corrected. Outliers CB (too high/low GBCs counts) are removed. Highly supported CB-GBC combinations (>=30 reads, >=5 UMIs and >=10 coverage) are filtered, 
  and used to compute a CB-GBC UMI table. 
- Clone calling. Unique GBCs-combinations in the CB-GBC UMI table are computed, and individual CBs
  are assigned to them. GBCs-combinations with one or more of their GBCs shared by other GBCs-combinations are filtered out (un-realistic combinations from multiple identical viral particles infecting independet cells). The remaining GBCs-combinations are considered clones.
- Cell assignment. CB assigned to the final pool of GBCs-combinations are assigned the respective
  clonal identity. Clonal frequencies are computed after this step, using assigned cell counts.     
