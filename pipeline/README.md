# MCM3 Project Analysis Pipeline
## Analysis Conditions

- bulk-RNAseq reads are quantified with Salmon and tested with limma-voom. DEGs
  require **BH-adjusted p-value <= 0.05** and **absolute log2 fold change >= 0.28**.
- CUT&RUN mouse and yeast reads are aligned separately with Bowtie2. Properly
  paired primary mouse alignments with **MAPQ >= 30** are retained.
- Mouse paired-fragment coverage is yeast-normalized to **10,000 retained yeast
  read pairs** per library. Factor and IgG libraries use the same alignment,
  filtering, and normalization procedure.
- For each factor, two biological-replicate BAMPE libraries are pooled and compared
  with pooled matched-IgG libraries using MACS3. Downstream peaks require a direct
  **q-value <= 0.05** and **fold enrichment >= 3**; no separate p-value threshold is
  applied.
- Promoters are **GENCODE M25 TSS +/- 1,000 bp**. A retained pooled peak must overlap
  a promoter by **at least 250 bp**. A promoter-bound gene additionally requires a
  yeast-normalized factor/IgG coverage ratio **>= 2** in **both biological replicates**.
- Regulatory targets are promoter-bound genes that are also DEGs after
  knockdown of the corresponding factor.

## Peak-Associated Gene Analysis

Whole-genome peak/gene Venns and pies use protein-coding-associated canonical
loci. Assign peaks to all GENCODE M25 gene promoters with >=250 bp overlap
(TSS +/-1 kb). Only when no promoter of any biotype qualifies, use the original
nearest-TSS assignment without a distance limit. Filter assignments to
protein-coding genes afterward; do not redirect excluded noncoding assignments.
Count distinct Ensembl IDs for genes and unique canonical loci for peaks.

Step 07 also intersects genes associated with the three-factor shared peaks
with the existing KD DEG lists. These candidate associations have data under
`regulatory_work/visuals/shared_peaks_visuals/data/` and figures under
`regulatory_work/visuals/shared_peaks_visuals/`, not substituted for the
promoter-bound regulatory targets. This branch includes non-promoter assignments;
nearest-TSS proximity alone does not demonstrate regulation.

Both branches are part of the standard pipeline. `--publication-figures` adds
text-free counterparts in a `publication_figures/` subfolder of each figure folder.
The [CUT&RUN](scripts/CUT&RUN/README.md) and
[regulatory](scripts/regulatory_target/README.md) READMEs describe the inputs,
assignment audits, and figure meanings.
