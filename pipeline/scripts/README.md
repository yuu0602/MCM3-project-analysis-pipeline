# Pipeline Scripts

| Stage | Script | Result |
| --- | --- | --- |
| 01 | `bulkRNAseq/01_prepare_inputs.py` | References, manifests, and accepted/raw input staging |
| 02 | `bulkRNAseq/02_quantify_and_test.py` | Salmon quantification and DEG tables |
| 03 | `bulkRNAseq/03_render_figures.py` | RNA-seq figures and optional IGV tracks |
| 04 | `CUT&RUN/04_align_and_normalize.py` | Filtered BAMs and yeast-normalized tracks |
| 05 | `CUT&RUN/05_call_peaks_and_define_binding.py` | Retained peaks and promoter-bound genes |
| 06 | `CUT&RUN/06_render_figures.py` | Main protein-coding gene and peak CUT&RUN branches; optional promoter analysis |
| 07 | `regulatory_target/07_define_regulatory_targets.py` | Main protein-coding gene and peak regulatory branches; optional promoter analysis |
