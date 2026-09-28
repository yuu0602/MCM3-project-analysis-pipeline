# Regulatory Analyses

| Stage | Script | Purpose |
| --- | --- | --- |
| 07 | `07_define_regulatory_targets.py` | Generate promoter-bound targets; optionally generate shared-peak-associated DEG figures. |

The default run generates only promoter-bound regulatory outputs. To also write
shared-peak figures and source data under
`regulatory_work/visuals/shared_peaks_visuals/`, run:

```bash
python pipeline/scripts/regulatory_target/07_define_regulatory_targets.py \
  --shared-peaks --publication-figures
```

Use `--shared-peaks-only` to update only that optional branch. Omitting both
shared-peak options leaves existing shared-peak outputs untouched.
