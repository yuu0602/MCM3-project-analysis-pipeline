# Pipeline Scripts

`run_pipeline.py` executes the numbered stages in dependency order.

| Stage | Script | Result |
| --- | --- | --- |
| 01 | `bulkRNAseq/01_prepare_inputs.py` | References, manifests, and accepted/raw input staging |
| 02 | `bulkRNAseq/02_quantify_and_test.py` | Salmon quantification and DEG tables |
| 03 | `bulkRNAseq/03_render_figures.py` | RNA-seq figures and optional IGV tracks |
| 04 | `CUT&RUN/04_align_and_normalize.py` | Filtered BAMs and yeast-normalized tracks |
| 05 | `CUT&RUN/05_call_peaks_and_define_binding.py` | Retained peaks and promoter-bound genes |
| 06 | `CUT&RUN/06_render_figures.py` | Main protein-coding gene and peak CUT&RUN branches; optional promoter branch |
| 07 | `regulatory_target/07_define_regulatory_targets.py` | Main protein-coding gene and peak regulatory branches; optional promoter branch |

Steps 06 and 07 preserve both `all_protein_coding_genes_visuals/shared_locus/` and
`all_protein_coding_genes_visuals/non-shared_locus/`. Use `--gene-membership shared_locus`
or `--gene-membership non-shared_locus` to render only one. The default `both`
renders both separately. They share the same peak and RNA-seq thresholds;
only gene-set membership differs.

Steps 06 and 07 also preserve `all_protein_coding_peaks_visuals/shared_locus/` and
`all_protein_coding_peaks_visuals/non-shared_locus/`. Select either with `--peak-membership`,
or use `both` (default). The CUT&RUN peak figures use identical protein-coding-associated genomic peak
sets; regulatory integration either requires triple-shared loci up front or
uses each factor's own retained peaks with its corresponding KD DEGs.

```bash
python pipeline/scripts/run_pipeline.py --step all --source raw --publication-figures
```

Use `--source raw` to rebuild from raw sequencing files. By default, Steps 06
and 07 write the main results to `all_protein_coding_genes_visuals` and
`all_protein_coding_peaks_visuals` under their respective `visuals/` folders.
Use `--branch all-with-promoters` to also write the optional
`promoters_visuals/` branch. Each branch contains its own data or figure
provenance and optional text-free figures.

See [CUT&RUN](CUT&RUN/README.md) and
[Regulatory Target Analyses](regulatory_target/README.md) for branch meanings.
