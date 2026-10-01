# CUT&RUN

| Stage | Script | Purpose |
| --- | --- | --- |
| 04 | `04_align_and_normalize.py` | Align mouse and yeast reads, filter BAMs, and generate normalized bigWigs. |
| 05 | `05_call_peaks_and_define_binding.py` | Call/filter peaks and define promoter-bound genes. |
| 06 | `06_render_figures.py` | Render the main protein-coding gene and peak figure branches; promoter figures are optional. |

## Figure Branches

Step 06 writes labeled figures directly to folders under
`cutrun_work/visuals/`. Each folder contains a `publication_figures/`
subfolder when `--publication-figures` is used. The default renders the two
protein-coding branches; the promoter branch is requested explicitly.

- `promoters_visuals/` (optional): promoter-bound gene Venn, gene-body profiles,
  metaprofile, and RPKM panels. These use the Step 05 promoter-bound sets.
- `all_protein_coding_genes_visuals/`: protein-coding genes associated with retained peaks,
  including a conventional gene Venn, feature pies, gene-body profiles,
  metaprofile, and RPKM panels.
- `all_protein_coding_peaks_visuals/`: protein-coding-associated whole-genome peak loci, including Venns,
  feature pies, the shared-peak distribution, and peak-locus profile panels.

## Peak-To-Gene Assignment

The all-gene branch contains two separate result folders:
`all_protein_coding_genes_visuals/shared_locus/` (3,049 triple shared-locus-associated genes)
and `all_protein_coding_genes_visuals/non-shared_locus/` (6,156 genes in all three independent
factor lists). Both retain the same figure design and peak-assignment rules.
Use `--gene-membership shared_locus` or `--gene-membership non-shared_locus`
to render one; the default `both` renders both.

In `shared_locus`, the three-circle figure is a schematic of locus-defined
gene categories, not a conventional gene Venn. A gene can belong to several
classes through different loci; region counts are not additive. Profiles use
these locus-defined gene groups. `non-shared_locus` uses disjoint gene-list
intersections as described below. Their membership tables and matrices are
stored separately under the corresponding membership-named data subfolders.

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
In `non-shared_locus`, `Venn_PeakAssociatedGenes.png` intersects the independent per-factor gene lists.
Its seven regions are mutually exclusive gene categories. A shared gene can
have different associated peak locations for different proteins; overlapping
genomic loci are not required. The current triple intersection is 6,156 genes.

The shared-gene profiles and RPKM panels use these same gene-list intersections.
The triple group is also the gene universe used in the all-gene regulatory
analysis. Profiles align gene bodies at TSS/TES; they are not peak-centered
profiles. Feature pies summarize each factor's complete unique gene
associations; the all-gene metaprofile summarizes their union. Nearest-TSS
assignment is an association, not proof of direct regulation.

Source data and assignment audits are under
`cutrun_work/data/figure_inputs/protein_coding_peak_associations/`.
Each membership subfolder's `Venn_PeakAssociatedGenes_members.tsv` and `Venn_PeakAssociatedGenes_counts.tsv`
provide the plotted gene lists and counts. Separate `LocusDefined*` tables
retain the stricter locus-first associations for audit; they are not the
independent-gene Venn or profile membership. They define the `shared_locus`
version. The all-peak branch remains locus-first.

## Regeneration

The peak branch also keeps two folders:
`all_protein_coding_peaks_visuals/shared_locus/` and `all_protein_coding_peaks_visuals/non-shared_locus/`.
Use `--peak-membership shared_locus`, `--peak-membership non-shared_locus`,
or `--peak-membership both` (default).

Both CUT&RUN folders show the same protein-coding-associated peak sets and genomic-overlap Venn.
Only loci with at least one existing `protein_coding` assignment are retained;
excluded noncoding assignments are never redirected to another gene. The current
totals are MCM3 19,588, NONO 15,880, PSPC1 29,013, and IgG 273, with 3,036
three-factor shared loci. These are canonical overlap segments, not original
MACS peak-call counts. Unfiltered Step 05 loci remain available as source data.
The upstream thresholds, peak identities, peak-to-gene assignments, and figure
design do not change. The distinction is the downstream regulatory eligibility:
shared-locus integration starts with triple-shared loci; independent integration
starts with each protein's own retained loci. Different loci assigned to the
same gene are never relabeled as the same shared peak.

Each peak result has its own `data/` folder containing counts, canonical loci,
membership records, provenance, and profile matrices. The existing peak profile
design aligns protein-coding promoter-associated gene bodies to TSS/TES; this
profile subset does not change the whole-genome Venn or peak-pie denominators.

```bash
python 'pipeline/scripts/CUT&RUN/06_render_figures.py' --publication-figures
```

Use `--branch promoters`, `--branch all-genes`, or `--branch all-peaks` to
render one branch. The default `--branch all` renders the two main
protein-coding branches. Use `--branch all-with-promoters` to include the
optional promoter branch. Text-free figures reuse existing deepTools matrices
during the second pass.
