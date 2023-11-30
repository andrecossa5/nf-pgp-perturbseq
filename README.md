# Beta version: bulk_is_the_key

A Nextflow pipeline for lentiviral-barcoding single-cell lineage tracing (scLT) experiments data pre-processing.

Major steps:

BULK DNA libraries

- Extract GBC reads (<=1 mismatch from DNA anchor sequence).
- Graph-based correction. Unique GBC sequences (>1 read) are clustered, and read counts from "degenerate"
  sequences are passed to "correct" ones. Spikeins are removed, and remaining GBC sequences are filtered 
  to retain only "correct" and "uniquely observed" sequences, with read_counts > bulk_min_n_reads.
- Reference GBC pool construction. A reference GBC whitelist is obtained from all samples.
- (Optional) clone calling step from bulk DNA. Two methods implemented: a. uses spike-ins
  information, and b. does not, relying only distributional properties of the GBC sequences
  read counts. 

Note: MOI of lentiviral infection must be very low (e.g., <= 0.1) to retrieve clones from bulk DNA-sequencing, to ensure a GBCs-clones one-to-one mapping. 

SC libraries (10X and GBC reads)

- Solo alignment of 10x reads
- GBC reads elements (i.e., CBC, GBC and UMI) extraction (<= 1 mismatch from cDNA anchor sequence) extraction.
  Only reads from putative cells (Solo-corrected CBC, passing Solo cell calling) are retained in this step.
- Bulk-reference-based GBC correction, filtering and CBC-GBC table computation. Single-cell CBC-UMI-GBC
  combinations with GBC matching a unique, "correct" bulk GBC (1 max mismatches allowed) are filtered and (if necessary) corrected to the bulk reference sequence. Highly supported CB-GBC combinations (>=15 reads, >=5 UMIs, >=3 coverage
  and umi counts ratio with most abundant GBC for a certain CBC >=.3) are filtered, and used to compute a CB-GBC
  UMI table. 
- Clone calling. Unique GBCs-combinations are from computed the CB-GBC UMI table, and individual CBs
  are assigned to them.
- Cell assignment. CB assigned to the final pool of GBCs-combinations are assigned the respective
  clonal identity. Clonal frequencies are computed after this step, using assigned cell counts.     
