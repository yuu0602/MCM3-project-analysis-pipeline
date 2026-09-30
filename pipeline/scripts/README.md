# Pipeline Scripts

`run_pipeline.py` executes the numbered stages in dependency order.

| Stage | Script | Result |
| --- | --- | --- |
| 01 | `bulkRNAseq/01_prepare_inputs.py` | References, manifests, and accepted/raw input staging |
| 02 | `bulkRNAseq/02_quantify_and_test.py` | Salmon quantification and DEG tables |
| 03 | `bulkRNAseq/03_render_figures.py` | RNA-seq figures and optional IGV tracks |
| 04 | `CUT&RUN/04_align_and_normalize.py` | Filtered BAMs and yeast-normalized tracks |
| 05 | `CUT&RUN/05_call_peaks_and_define_binding.py` | Retained peaks and promoter-bound genes |
| 06 | `CUT&RUN/06_render_figures.py` | Promoter-gene, all-gene, and all-peak CUT&RUN branches |
| 07 | `regulatory_target/07_define_regulatory_targets.py` | Promoter-bound, shared-locus gene, and shared-locus peak regulatory branches |

```bash
python pipeline/scripts/run_pipeline.py --step all --source raw --publication-figures
```

Use `--source raw` to rebuild from raw sequencing files. Step 06 writes to
`cutrun_work/visuals/{promoters_visuals,all_genes_visuals,all_peaks_visuals}`.
Step 07 writes to
`regulatory_work/visuals/{promoters_visuals,all_genes_visuals,all_peaks_visuals}`. Each branch
contains its own data or figure provenance and optional text-free figures.

See [CUT&RUN](CUT&RUN/README.md) and
[Regulatory Target Analyses](regulatory_target/README.md) for branch meanings.
