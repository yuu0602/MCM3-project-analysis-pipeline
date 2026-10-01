# MCM3 Project Analysis Pipeline

## Analysis Conditions

- Bulk RNA-seq reads are quantified with Salmon and tested with limma-voom.
  DEGs require BH-adjusted p-value <=0.05 and absolute log2 fold change >=0.28.
- CUT&RUN mouse and yeast reads are aligned separately with Bowtie2. Properly
  paired primary mouse alignments with MAPQ >=30 are retained.
- Mouse paired-fragment coverage is yeast-normalized to 10,000 retained yeast
  read pairs per library.
- Pooled matched-IgG MACS3 peaks require q-value <=0.05 and fold enrichment
  >=3.
- Promoters are GENCODE M25 TSS +/-1,000 bp. Promoter binding requires a retained
  peak overlap >=250 bp and yeast-normalized factor/IgG coverage ratio >=2 in
  both biological replicates.
