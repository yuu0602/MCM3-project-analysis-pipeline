<h1 style="font-family: Arial, Helvetica, sans-serif; font-size: 2.2em; line-height: 1.2; color: #17365d; margin-bottom: 0.15em;">Scripts directory</h1>

<p style="font-family: Arial, Helvetica, sans-serif; font-size: 1.05em; line-height: 1.6; color: #4b5563;"><code>&lt;PROJECT_ROOT&gt;/pipeline/scripts/</code> — executable stages, shared configuration, and assay-specific rendering helpers for the MCM3–NONO–PSPC1 workflow.</p>

The numbered filenames indicate the intended execution order; the top-level [`run_pipeline.py`](run_pipeline.py) provides dependency-aware orchestration.

> Replace `<PROJECT_ROOT>` with the local path containing the repository clone. Run commands from `<PROJECT_ROOT>/pipeline` unless a command explicitly states otherwise.

## Structure

```text
scripts/
├── config.py                        # Shared paths, factors, thresholds, and constants
├── run_pipeline.py                  # Dependency-aware workflow launcher
├── bulkRNAseq/
│   ├── 01_prepare_inputs.py          # References, manifests, and input validation
│   ├── 02_quantify_and_test.py       # Salmon quantification and limma-voom DEG analysis
│   └── 03_render_figures.py          # RNA-seq figures and optional IGV tracks
├── CUT&RUN/
│   ├── 04_align_and_normalize.py     # Alignment, filtering, and BigWig construction
│   ├── 05_call_peaks_and_define_binding.py
│   ├── 06_render_figures.py
│   └── figures_rendering/            # CUT&RUN plotting and annotation helpers
└── regulatory_target/
    ├── 07_define_regulatory_targets.py
    └── figures_rendering/            # Regulatory direction, Venn, and profile renderers
```

## Shared entry points

| File | Responsibility | Accomplishes |
| --- | --- | --- |
| `config.py` | Central configuration | Defines project-relative workspaces, references, factor names, colors, and analysis thresholds. Altering this file changes the declared analytical specification. |
| `run_pipeline.py` | Workflow orchestration | Resolves stage order for `prepare`, `rna`, `cutrun`, `regulatory`, `figures`, or `all`; passes raw/accepted source mode and thread settings to compatible stages. |

## Bulk RNA-seq scripts

| Script | Inputs | Processing and products |
| --- | --- | --- |
| `01_prepare_inputs.py` | Raw FASTQ directories or packaged accepted intermediates. | Builds sample manifests, prepares references and indexes, and validates the upstream material required by later stages. |
| `02_quantify_and_test.py` | RNA-seq manifest, Salmon index, and transcript/gene annotation. | Quantifies expression with Salmon and generates factor-specific differential-expression tables using limma-voom. |
| `03_render_figures.py` | DEG tables and, optionally, alignments. | Produces RNA-seq plots and optional individual/mean browser tracks. |

## CUT&RUN scripts

### `04_align_and_normalize.py`

Consumes paired-end CUT&RUN libraries described in `cutrun_work/metadata/Samples.tsv`. It aligns reads independently to mouse and yeast with Bowtie2, creates filtered indexed BAMs, counts retained yeast read pairs, and scales mouse paired-fragment coverage to 10,000 yeast pairs. It then creates one normalized BigWig per library, a two-replicate mean track for each factor and IgG, and display-scaled tracks for browsers and figures.

During mean-track creation, a `pyBigWig` chromosome query with no stored intervals is explicitly treated as an empty interval list. This prevents empty chromosomes or contigs from terminating the analysis while preserving zero signal in the calculated mean.

### `05_call_peaks_and_define_binding.py`

Uses Step 04 filtered BAMs for pooled factor-versus-matched-IgG MACS3 calls, then retains peaks that satisfy the configured direct FDR and fold-enrichment criteria. It intersects those peaks with GENCODE M25 promoter windows and queries the Step 04 replicate BigWigs to calculate factor-to-IgG coverage ratios. A gene is designated promoter-bound only when it meets the peak/promoter criterion and ratio requirement in both biological replicates.

### `06_render_figures.py`

Consumes Step 04 display-scaled mean tracks and Step 05 figure-input tables. It configures deepTools `computeMatrix` and rendering helpers to produce gene-body and peak metaprofiles, coverage figures, peak-overlap Venn diagrams, and genomic peak-distribution panels. Results are written to `cutrun_work/visuals/` with source matrices and tables retained under `cutrun_work/data/figure_inputs/`.

## Regulatory-target scripts

### `07_define_regulatory_targets.py`

Combines Step 05 promoter-bound gene sets with factor-matched RNA-seq DEG tables. It produces direct-target lists and directional classifications, creates the three-factor regulatory-target Venn diagram, and calls its R renderers to generate group-level composition charts and CUT&RUN metaprofiles. The profile component uses Step 04 display-scaled mean BigWigs.

## Running scripts safely

- Use `run_pipeline.py --source raw` for a dependency-aware raw reconstruction.
- Invoke Steps 04 and 05 with `--from-raw` only when rebuilding CUT&RUN intermediates from FASTQs.
- Do not submit Steps 04–07 concurrently: each depends on files created by prior stages.
- Run `--dry-run` where supported to inspect planned commands without changing data.
- Preserve `environment.yml`, the exact command line, and script revision with each reproducibility record.

## Related documentation

- [Main workflow overview and installation](../../README.md).
- [CUT&RUN workspace products](../cutrun_work/README.md).
- [Regulatory-target workspace products](../regulatory_work/README.md).
