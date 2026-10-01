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

Step 06 keeps three analyses separate. The two protein-coding branches are the
default publication outputs; the promoter branch is optional:

- `promoters_visuals` (optional): genes satisfying the promoter-binding rule;
- `all_protein_coding_genes_visuals`: protein-coding genes assigned to all retained peaks;
- `all_protein_coding_peaks_visuals`: retained genomic peak loci with at least one protein-coding gene assignment.

All-gene assignment gives priority to every promoter overlapping a peak by
>=250 bp, then uses nearest-TSS fallback only when no promoter qualifies.
Protein-coding filtering is applied after assignment. Genes are counted by
Ensembl ID; peaks are counted by canonical genomic locus.
`Venn_PeakAssociatedGenes.png` is a conventional Venn of independent per-factor
gene lists, with mutually exclusive regions. Its current triple intersection
contains 6,156 genes. Associated peaks need not overlap between proteins.
All-gene profiles use the same gene-list intersections. Same-locus analyses
remain separate in `all_protein_coding_peaks_visuals`.

## Regulatory Integration

The peak branch has `all_protein_coding_peaks_visuals/shared_locus/` and
`all_protein_coding_peaks_visuals/non-shared_locus/` in both CUT&RUN and regulatory outputs.
Steps 06/07 accept `--peak-membership shared_locus`, `non-shared_locus`, or
`both` (default). CUT&RUN peak figures remain identical: peaks and their
genomic-overlap definition are unchanged. Regulatory membership differs:
the shared-locus branch starts with triple-shared peaks, while the independent
branch requires each factor's retained peak plus its corresponding KD DEG
association without a triple-sharing prefilter. Gene-level overlap never
substitutes for genomic peak overlap in the peak Venn.
Both peak branches exclude loci lacking a protein-coding assignment from all
visual and regulatory backgrounds. The plotted totals are MCM3 19,588,
NONO 15,880, PSPC1 29,013, and IgG 273; the three-factor intersection is 3,036.
Existing noncoding assignments are excluded, not redirected to coding genes.

Both CUT&RUN and regulatory `all_protein_coding_genes_visuals` contain `shared_locus/` and
`non-shared_locus/`. The former retains genes assigned after genomic peak
intersection (3,049 triple-associated genes); the latter intersects independent
factor gene lists (6,156). `--gene-membership both` is the default for Steps
06 and 07. Select either membership name to render only that version. The
shared-locus circle diagram is schematic with non-disjoint gene categories;
the independent-gene diagram is a conventional gene Venn.

Step 07 uses the same current DEG tables for three separate analyses:

- `regulatory_work/visuals/promoters_visuals`: promoter-bound genes AND DEGs;
- `regulatory_work/visuals/all_protein_coding_genes_visuals`: protein-coding genes assigned to
  all three factors independently AND each knockdown's DEGs.
- `regulatory_work/visuals/all_protein_coding_peaks_visuals`: protein-coding-associated peak loci
  with at least one assigned protein-coding DEG per knockdown; counted once
  per locus. Group A-D profiles align to peak centers; the MCM3 target panel
  instead profiles unique associated Up/Down DEG gene bodies from TSS to TES,
  excluding mixed-direction peak assignments from that panel only.

Each all-gene regulatory version uses its own common associated-gene universe for all
knockdowns; Group A is DEG in all three, and Groups B-D describe pairwise DEG
overlap. Its input gene universe matches the triple-shared class of the CUT&RUN
gene Venn. Regulatory all-gene candidates
include nearest-TSS assignments outside promoters and should not be called
established direct regulation.

In the peak branch, different assigned genes may support the same peak in
different KD comparisons. Peak counts are not gene counts or differential
binding calls. Up/Down/Mix summarizes the directions of linked RNA-seq DEGs.
Run Step 07 with `--branch all-peaks` for only this branch. The default
`--branch all` runs the two protein-coding branches; `--branch all-with-promoters`
also renders the optional promoter branch, while `--branch both` retains the two
gene branches only.

`--publication-figures` adds text-free counterparts within each branch's
`publication_figures/` folder.

See [CUT&RUN](scripts/CUT&RUN/README.md) and
[Regulatory Target Analyses](scripts/regulatory_target/README.md).
