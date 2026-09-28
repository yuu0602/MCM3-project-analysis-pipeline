# Pipeline

`run_pipeline.py` executes the numbered stages in dependency order.

| Stage | Script | Result |
| --- | --- | --- |
| 01 | `bulkRNAseq/01_prepare_inputs.py` | References, manifests, and local inputs |
| 02 | `bulkRNAseq/02_quantify_and_test.py` | DEG tables and bulk-RNAseq figures |
| 03 | `bulkRNAseq/03_render_figures.py` | Optional RNA-seq IGV tracks (`--igv-tracks`): local-SSD HISAT2 index build, MAPQ >=30/proper-pair filtering, CPM bigWigs, and group means |
| 04 | `CUT&RUN/04_align_and_normalize.py` | Alignments, filtered BAMs, and normalized tracks |
| 05 | `CUT&RUN/05_call_peaks_and_define_binding.py` | Peaks and promoter-bound genes |
| 06 | `CUT&RUN/06_render_figures.py` | CUT&RUN figures and promoter-priority peak-to-gene assignments |
| 07 | `regulatory_target/07_define_regulatory_targets.py` | Promoter-bound targets; shared-peak DEG figures only with `--shared-peaks` |

Step 06 assigns each canonical peak to all gene promoters overlapping by >=250 bp
within TSS +/-1 kb; only peaks with no qualifying promoter use nearest-TSS
fallback. Protein-coding filtering follows assignment. See [CUT&RUN](CUT&RUN/README.md).

Step 07 retains the promoter-bound target analysis and, with `--shared-peaks`, writes the
shared-peak-associated DEG data under `regulatory_work/visuals/shared_peaks_visuals/data/`
and figures under `regulatory_work/visuals/shared_peaks_visuals/`.
See [Regulatory Analyses](regulatory_target/README.md) for their distinct meanings.

The main runner accepts `--shared-peaks` and forwards it only to Step 07.
Combine it with `--publication-figures` for text-free counterparts. Existing
shared-peak results are not removed when this option is omitted.
