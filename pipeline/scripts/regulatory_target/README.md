# Regulatory Analyses

| Stage | Script | Purpose |
| --- | --- | --- |
| 07 | `07_define_regulatory_targets.py` | Generate promoter-bound regulatory targets and shared-peak-associated DEG figures. |

Run after bulk RNA-seq and CUT&RUN Steps 05-06. Both branches use the existing
directional DEG lists; neither reruns differential expression.

## Promoter-Bound Targets

Existing factor-specific promoter-bound gene sets are intersected with the
corresponding KD DEG lists. Outputs remain in `regulatory_work/data/` and
`regulatory_work/visuals/`. Their thresholds and membership are unchanged by
the peak-to-gene assignment update.

## Shared-Peak-Associated DEGs

Data are in `regulatory_work/visuals/shared_peaks_visuals/data/`.
Figures are in `regulatory_work/visuals/shared_peaks_visuals/`, with text-free
counterparts in its `publication_figures/` subfolder.
Start with the three-factor intersection of the coding-associated canonical
peaks from Step 06. Use all retained promoter-priority assignments, including
nearest-TSS fallback only for peaks with no qualifying promoter. Intersect the
resulting distinct Ensembl gene IDs with each KD DEG list (BH-adjusted p <=0.05,
absolute log2FC >=0.28). These are candidate gene associations, not established
promoter-bound direct targets.

- A: DEG in all three KDs.
- B: DEG in MCM3 and PSPC1 only.
- C: DEG in NONO and PSPC1 only.
- D: DEG in MCM3 and NONO only.

"Only" describes DEG status, not factor-exclusive binding: all branch input
genes are associated with three-factor shared peaks. Up/Down/Mix describes
agreement of the indicated RNA-seq effects. The branch includes directional
overlap panels, target Venns, group pies, DEG-association pies, and gene-body
metaprofiles in the established designs. Profiles use existing mean bigWigs,
25-bp bins, 3-kb flanks, and a 3-kb scaled gene body.

```bash
python pipeline/scripts/regulatory_target/07_define_regulatory_targets.py --publication-figures
```

Use `--shared-peaks-only` to update only the new branch. Text-free counterparts
are generated in `publication_figures/` within each branch's figure folder when
`--publication-figures` is supplied. Omit that option for labeled figures only.
No testing-directory inputs are required by either branch.
