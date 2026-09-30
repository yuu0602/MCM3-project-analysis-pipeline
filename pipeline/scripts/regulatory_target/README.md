# Regulatory Target Analyses

| Stage | Script | Purpose |
| --- | --- | --- |
| 07 | `07_define_regulatory_targets.py` | Integrate current directional DEGs with promoter-bound genes, shared-peak-associated genes, and shared genomic peaks. |

Step 07 uses the existing bulk RNA-seq DEG tables and does not rerun
differential expression. It generates three separate branches under
`regulatory_work/visuals/`.

## Promoter Branch

`promoters_visuals/` intersects each factor's promoter-bound genes with DEGs
from knockdown of the same factor. This is the stricter direct-target analysis.
The promoter-bound rule is defined in CUT&RUN Step 05.

## All-Gene Branch

`all_genes_visuals/` starts with the three-factor shared genomic loci in the
CUT&RUN peak Venn (`region_mask=7`), then uses Step 06's protein-coding gene
assignments to those loci. All promoters overlapping a peak by >=250 bp at
TSS +/-1 kb are assigned; nearest-TSS fallback applies only when no promoter
qualifies. Biotype filtering happens after assignment.

The same shared-peak gene universe is intersected with each knockdown's DEGs.
Ensembl IDs without version suffixes are used throughout; symbols are labels.
This gene universe matches the triple-shared class in CUT&RUN's
`Venn_PeakAssociatedGenes.png`: membership requires association with a common
genomic peak locus, not separate peaks assigned to the same gene.
Nearest-TSS associations remain candidates rather than proof of direct regulation.
`data/SharedPeakLoci.tsv`, `SharedPeakGeneAssignments.tsv`,
`SharedPeakAssociatedGenes.tsv`, and `SharedPeakParameters.json` record the
source loci, assignments, gene universe, and input checksums.

Both gene branches contain:

- directional CUT&RUN-gene versus Up/Down DEG panels;
- a three-factor target Venn;
- Group A-D direction pies and combined pie panel;
- NONO and PSPC1 target-association pies;
- Group A-D and MCM3-target gene-body metaprofiles;
- a `data/` folder with gene lists, counts, BEDs, matrices, and parameters;
- text-free counterparts in `publication_figures/` when requested.

Groups A-D describe the DEG intersections: A, all three; B, MCM3+PSPC1 only;
C, NONO+PSPC1 only; D, MCM3+NONO only. In the all-gene branch, every input gene
has three-factor shared-peak association, so "only" refers to knockdown DEG
status, never exclusive CUT&RUN binding. Up/Down/Mix reports RNA-seq direction
agreement among the indicated knockdowns. Metaprofiles average the selected
gene bodies, not the shared peak centers.

## All-Peak Branch

`all_peaks_visuals/` uses all three-factor shared canonical loci (`region_mask=7`)
as its background. Each peak enters a KD-associated set when at least one of
its assigned protein-coding genes is a DEG for that KD. Existing Step 06 gene
assignments and current RNA-seq thresholds are unchanged. A peak is counted
once per KD regardless of how many genes it maps to. Different genes assigned
to the same peak may support its membership in different KD sets.

The target Venn and Groups A-D therefore count **peaks**, not genes. All loci
already overlap all three CUT&RUN factors; the Venn groups describe linked-gene
DEG status, not factor-exclusive binding or differential peak occupancy.
Within one KD, a peak is Up or Down if all its linked DEGs agree, and Mix if it
links to both directions. A group's Up/Down classification requires agreement
across every indicated KD; all other combinations are Mix. "No DEG" means no
assigned significant DEG, including loci without coding assignments or genes
not tested. It is not proof of no biological effect.

This branch includes:

- `Venn_target.png`: overlap of KD-associated shared-peak sets;
- `Pie_DEGAssociation_<factor>.png`: Up/Down/Mix among DEG-linked shared
  peaks only; No-DEG peaks are excluded from these pies but remain in the
  complete source tables;
- `peak_vs_DEG_direction_<factor>.png` and `peak_vs_DEG_direction_all.png`:
  counterparts of the gene-branch directional circle panels, in the same
  design. The center is the three-factor shared peak set. Side sets contain
  all retained peaks of the indicated factor linked only to Up or only to Down
  DEGs. Mixed-direction links are excluded from both sides; the center-only
  count includes No DEG and Mix, with their breakdown in the data table.
  Every count is a peak count, not a DEG-gene count;
- Group A-D pies and `Pie_Groups.png`, using the existing gene-branch design;
- `Pie_TargetBinding_NONO.png` and `Pie_TargetBinding_PSPC1.png`: fractions of
  all shared peaks with versus without a DEG link, not fractions of DEG genes;
- `Metaprofile_GroupA.png` through `Metaprofile_GroupD.png`: arithmetic mean coverage aligned to peak
  centers over +/-3 kb, 25-bp bins, existing mean bigWigs, missing coverage zero;
- `Metaprofile_MCM3_target.png`: strand-aware TSS-to-TES profiles of unique
  MCM3 Up/Down DEG genes assigned to the respective Up-only/Down-only shared
  peaks. Gene bodies are scaled to 3 kb, with 3-kb flanks and 25-bp bins;
  existing mean bigWigs, equal gene weights, and missing coverage zero;
- `data/PeakGeneDEGEvidence.tsv`, `Peak_DEG_directions.tsv`, group peak lists,
  BEDs, per-peak matrices, counts, and `AnalysisParameters.json` with input hashes;
- `_noTexts` versions in `publication_figures/` when requested.

The MCM3 target profile compares DEG genes linked to Up-only and Down-only peaks.
Each gene is counted once regardless of its number of linked peaks. Non-DEG
genes are not assigned their neighboring DEG's direction. Mixed-direction
peak assignments are excluded from this profile and listed in its exclusion table;
they remain in regulatory membership, group figures, and complete source tables.
These are RNA-seq directions of linked genes, never CUT&RUN gain/loss calls.
The `Metaprofile_*` filenames are consistent across branches, but peak-branch
x axes use the peak center except for the MCM3 target gene-body profile.
Its gene lists, peak-to-gene assignments, stranded BEDs, and per-gene matrix
are saved under `data/Metaprofile_MCM3_target*`. Old `PeakProfile_*` names are no
longer generated. The two `Pie_TargetBinding_*` figures still use all shared
peaks as their denominator because they assess the fraction with a DEG link.

```bash
python pipeline/scripts/regulatory_target/07_define_regulatory_targets.py --publication-figures
```

Use `--branch promoters`, `--branch all-genes`, or `--branch all-peaks` to render
one branch. The default `--branch all` renders all three; `--branch both` remains
available for the two gene branches only.

```bash
python pipeline/scripts/regulatory_target/07_define_regulatory_targets.py --branch all-peaks --publication-figures
```
