# CUT&RUN

| Stage | Script | Purpose |
| --- | --- | --- |
| 04 | `04_align_and_normalize.py` | Align mouse and yeast reads, filter BAMs, and generate normalized bigWigs. |
| 05 | `05_call_peaks_and_define_binding.py` | Call and filter peaks, then define promoter-bound genes. |
| 06 | `06_render_figures.py` | Render CUT&RUN figures and generate promoter-priority peak-to-gene assignments. |

## Peak and Associated-Gene Figures

`Venn_Peaks`, `Venn_PeakAssociatedGenes`, `Pie_PeakDistribution`,
`Pie_PeakAssociatedGenes`, and the individual factor/IgG pies show
protein-coding-associated subsets. Peak-calling thresholds, original canonical
locus coordinates, normalization, and promoter-binding criteria are unchanged.

Assignment starts from the original all-biotype canonical loci, not a previously
filtered protein-coding subset:

1. Assign each peak to **all** GENCODE M25 gene promoters overlapping it by
   at least 250 bp. Promoters are gene TSS +/-1,000 bp, strand-aware.
2. Only if no promoter from any biotype qualifies, use the original nearest-TSS
   assignment across all biotypes, with no distance cap.
3. Filter assignment records to `protein_coding`. Do not redirect a noncoding
   fallback assignment to a more distant coding gene.

A peak with at least one retained coding assignment is counted once per assay.
Genes are counted by distinct Ensembl gene IDs with version suffixes removed;
symbols are labels only. A peak may be associated with several genes. For a gene
pie, select its assigned peak with the smallest midpoint-to-TSS distance, breaking
ties by chromosome, start and end. Peak pies use each unique retained locus.
The existing midpoint-based genomic-feature categories and plot designs remain
unchanged; these categories are distinct from the >=250-bp assignment rule.

Gene Venn intersections mean common gene associations, not necessarily a common
binding site. Peak Venn intersections use canonical locus identities.
Nearest-TSS association alone does not establish direct regulation.

Source data are under
`cutrun_work/data/figure_inputs/protein_coding_peak_associations/`:

- `PeakGeneAssignments_all_biotypes_audit.tsv`: all assignment records before filtering.
- `PeakGeneAssignments.tsv` and factor-specific counterparts: retained many-to-many links.
- `Peaks_<factor>.tsv`: unique retained loci.
- `PeakAssociatedGenes_<factor>.tsv`: representative records per Ensembl gene ID.
- `AssignmentSummary.tsv`, Venn membership, and pie counts: figure inputs and totals.

`Pie_SharedPeaks_Distribution.png` uses the three-factor peak intersection.
`MCM3_Peak_profiles` and `MCM3_NONO_PSPC1_peak_profiles` use the same filtered
Venn loci, restricted to those overlapping a promoter by >=250 bp. The existing
maximum-overlap rule chooses one promoter per locus. Each peak contributes one
TSS/TES-scaled gene-body row; these are not whole-genome peak-centered averages.

## Regeneration

From the repository root, using existing peak calls and bigWigs:

```bash
python 'pipeline/scripts/CUT&RUN/06_render_figures.py' --publication-figures
```

Labeled PNGs are in `cutrun_work/visuals/`; text-free counterparts are in
`visuals/publication_figures/`. Omit `--publication-figures` for labeled figures
only. Use `--peak-associations-only` for Venns and pies, then
`--peak-profiles-only` for the two affected peak-associated profiles.

Step 07 consumes these assignments for the shared-peak-associated DEG branch
with data under `regulatory_work/visuals/shared_peaks_visuals/data/` and figures under
`regulatory_work/visuals/shared_peaks_visuals/`. The original promoter-bound
regulatory analysis remains separate under `regulatory_work/data` and `visuals`.
