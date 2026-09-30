# CUT&RUN

| Stage | Script | Purpose |
| --- | --- | --- |
| 04 | `04_align_and_normalize.py` | Align mouse and yeast reads, filter BAMs, and generate normalized bigWigs. |
| 05 | `05_call_peaks_and_define_binding.py` | Call/filter peaks and define promoter-bound genes. |
| 06 | `06_render_figures.py` | Render promoter-gene, all-gene, and all-peak figure branches. |

## Figure Branches

Step 06 writes labeled figures directly to three folders under
`cutrun_work/visuals/`. Each folder contains a `publication_figures/`
subfolder when `--publication-figures` is used.

- `promoters_visuals/`: promoter-bound gene Venn, gene-body profiles,
  metaprofile, and RPKM panels. These use the Step 05 promoter-bound sets.
- `all_genes_visuals/`: protein-coding genes associated with retained peaks,
  including a locus-defined three-circle summary, feature pies, gene-body profiles,
  metaprofile, and RPKM panels.
- `all_peaks_visuals/`: whole-genome retained peak loci, including Venns,
  feature pies, the shared-peak distribution, and peak-locus profile panels.

## Peak-To-Gene Assignment

All-gene assignment starts from the complete canonical peak loci:

1. Assign a peak to every GENCODE M25 promoter (TSS +/-1,000 bp) that it
   overlaps by at least 250 bp.
2. Only when no promoter of any biotype qualifies, assign the nearest TSS with
   no distance cap.
3. Retain `protein_coding` assignment records without redirecting excluded
   noncoding assignments.

Peak figures count unique genomic loci. Gene figures count distinct Ensembl
gene IDs after removing version suffixes; symbols are labels only. One peak may
map to several genes, so a gene intersection can exceed a peak intersection.
Shared-gene membership is defined from overlapping genomic loci first, then
gene assignments. `Venn_PeakAssociatedGenes.png` keeps the original circle
design as a schematic of locus-defined gene categories, not a conventional
gene-set Venn. Each region counts distinct genes associated with one exact
peak-overlap class. A gene may belong to several classes through different
loci, so the counts are not mutually exclusive and must not be summed as
unique genes. Circle labels report distinct per-factor gene totals, and circle
areas are schematic. "Only" refers
to the proteins present at the locus, not exclusive binding across the gene.

The shared-gene profiles and RPKM panels use these same locus-defined groups.
The triple group is also the gene universe used in the all-gene regulatory
analysis. Profiles align gene bodies at TSS/TES; they are not peak-centered
profiles. Feature pies summarize each factor's complete unique gene
associations; the all-gene metaprofile summarizes their union. Nearest-TSS
assignment is an association, not proof of direct regulation.

Source data and assignment audits are under
`cutrun_work/data/figure_inputs/protein_coding_peak_associations/`.
`LocusDefinedPeakGeneAssignments.tsv` links every class membership to its locus;
`LocusDefinedGeneGroups.tsv` and `LocusDefinedGeneGroups_counts.tsv` provide
the plotted gene lists and counts.

## Regeneration

```bash
python 'pipeline/scripts/CUT&RUN/06_render_figures.py' --publication-figures
```

Use `--branch promoters`, `--branch all-genes`, or `--branch all-peaks` to
render one branch. The default `--branch all` renders all three. Text-free
figures reuse existing deepTools matrices during the second pass.
