<h1 style="font-family: Arial, Helvetica, sans-serif; font-size: 2.2em; line-height: 1.2; color: #17365d; margin-bottom: 0.15em;">MCM3–NONO–PSPC1 Multi-omics Pipeline</h1>

<p style="font-family: Arial, Helvetica, sans-serif; font-size: 1.05em; line-height: 1.6; color: #4b5563;"><em>A workflow integrating CUT&amp;RUN chromatin profiling with bulk-RNAseq differential expression to define direct regulatory targets.</em></p>

[![Conda environment](https://img.shields.io/badge/Conda-mcm3--pipeline-44A833?logo=anaconda&logoColor=white)](pipeline/environment.yml)
[![Python](https://img.shields.io/badge/Python-3.11-3776AB?logo=python&logoColor=white)](pipeline/environment.yml)
[![Workflow](https://img.shields.io/badge/workflow-CUT%26RUN%20%2B%20RNA--seq-17365D)](#visual-abstract)

## Overview

This workflow characterizes the chromatin occupancy and transcriptional consequences of the MCM3, NONO, and PSPC1 complex.  Our pipeline utilizes paired-end CUT&RUN data with matched IgG controls and yeast spike-in normalization, then integrates factor-specific co-bound genes with bulk-RNAseq differential expression results. The direct regulatory targets are the genes at the intersection of the co-bound and differentially expressed genes.

## Visual abstract
