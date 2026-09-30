# MCM3 Project Analysis Pipeline

## Analysis Conditions

- Bulk RNA-seq reads are quantified with Salmon and tested with limma-voom.
  DEGs require BH-adjusted p-value <=0.05 and absolute log2 fold change >=0.28.
- CUT&RUN mouse and yeast reads are aligned separately with Bowtie2. Properly
  paired primary mouse alignments with MAPQ >=30 are retained.
- Mouse paired-fragment coverage is yeast-normalized to 10,000 retained yeast
  read pairs per library.
- Pooled matched-IgG MACS3 peaks require q-value <=0.05 and fold enrichment
  >=3; no separate p-value threshold is applied.
- Promoters are GENCODE M25 TSS +/-1,000 bp. Promoter binding requires a retained
  peak overlap >=250 bp and yeast-normalized factor/IgG coverage ratio >=2 in
  both biological replicates.

## CUT&RUN Figure Universes

Step 06 keeps three analyses separate:

- `promoters_visuals`: genes satisfying the promoter-binding rule;
- `all_genes_visuals`: protein-coding genes assigned to all retained peaks;
- `all_peaks_visuals`: retained genomic peak loci without gene filtering.

All-gene assignment gives priority to every promoter overlapping a peak by
>=250 bp, then uses nearest-TSS fallback only when no promoter qualifies.
Protein-coding filtering is applied after assignment. Genes are counted by
Ensembl ID; peaks are counted by canonical genomic locus.
Shared-gene classes are defined by overlapping genomic loci before assignment
to genes. `Venn_PeakAssociatedGenes.png` is a schematic three-circle summary
of locus-defined gene categories, not a conventional gene-set Venn: a gene can
belong to several locus classes, so region counts are not additive unique-gene
totals. Shared-gene profiles use these same
classes; "only" describes the locus, not exclusive binding across a gene.

## Regulatory Integration

Step 07 uses the same current DEG tables for three separate analyses:

- `regulatory_work/visuals/promoters_visuals`: promoter-bound genes AND DEGs;
- `regulatory_work/visuals/all_genes_visuals`: protein-coding genes assigned to
  three-factor shared peak loci AND each knockdown's DEGs.
- `regulatory_work/visuals/all_peaks_visuals`: three-factor shared peak loci
  with at least one assigned protein-coding DEG per knockdown; counted once
  per locus. Group A-D profiles align to peak centers; the MCM3 target panel
  instead profiles unique associated Up/Down DEG gene bodies from TSS to TES,
  excluding mixed-direction peak assignments from that panel only.

The all-gene regulatory branch uses the same shared-peak gene universe for all
knockdowns; Group A is DEG in all three, and Groups B-D describe pairwise DEG
overlap. Its input gene universe matches the triple-shared class of the CUT&RUN
locus-defined gene-group summary. Regulatory all-gene candidates
include nearest-TSS assignments outside promoters and should not be called
established direct regulation.

In the peak branch, different assigned genes may support the same peak in
different KD comparisons. Peak counts are not gene counts or differential
binding calls. Up/Down/Mix summarizes the directions of linked RNA-seq DEGs.
Run Step 07 with `--branch all-peaks` for only this branch; the default `all`
runs all three, while `both` retains the two gene branches.

`--publication-figures` adds text-free counterparts within each branch's
`publication_figures/` folder.

See [CUT&RUN](scripts/CUT&RUN/README.md) and
[Regulatory Target Analyses](scripts/regulatory_target/README.md).
