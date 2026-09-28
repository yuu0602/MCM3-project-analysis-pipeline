#!/usr/bin/env python3
"""Render promoter-bound targets and a separate shared-peak-associated DEG analysis."""

#Before you run this script, please replace directory placeholders with your own directory

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import CUTRUN_ROOT, REGULATORY_ROOT, RUN_ROOT


HELPERS = Path(__file__).resolve().parent / "figures_rendering"
VISUALS = REGULATORY_ROOT / "visuals"
PUBLICATION_FIGURES = VISUALS / "publication_figures"


def run_r(script: str, publication_figures: bool) -> None:
    environment = os.environ.copy()
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
    executable = shutil.which("Rscript", path=environment["PATH"])
    if executable is None:
        raise FileNotFoundError("Rscript")
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[variable] = "1"
    cache = REGULATORY_ROOT / "data" / ".matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(cache)
    command = [str(executable), str(HELPERS / script), str(RUN_ROOT)]
    if publication_figures:
        command.append(str(PUBLICATION_FIGURES))
    subprocess.run(command, check=True, env=environment)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--shared-peaks", action="store_true", help="Also generate regulatory_work/visuals/shared_peaks_visuals and its data")
    parser.add_argument("--publication-figures", action="store_true", help="Also render text-free PNGs in regulatory_work/visuals/publication_figures")
    parser.add_argument("--shared-peaks-only", action="store_true", help="Update only the shared-peak-associated DEG branch")
    args = parser.parse_args()
    include_shared = args.shared_peaks or args.shared_peaks_only
    if args.dry_run:
        print(f"[DRY-RUN] Promoter-bound targets enabled: {not args.shared_peaks_only}; shared-peak-associated DEGs enabled: {include_shared}")
        return
    if include_shared:
        assignments = CUTRUN_ROOT / "data/figure_inputs/protein_coding_peak_associations/PeakGeneAssignments.tsv"
        if not assignments.is_file():
            raise FileNotFoundError("Run CUT&RUN Step 06 before shared-peak integration")
    promoter_input = CUTRUN_ROOT / "data" / "figure_inputs" / "promoter_gene_venn"
    if not promoter_input.is_dir():
        raise FileNotFoundError("Run CUT&RUN Step 05 before regulatory integration")
    if args.publication_figures:
        PUBLICATION_FIGURES.mkdir(parents=True, exist_ok=True)
    if not args.shared_peaks_only:
        run_r("render_direction.R", args.publication_figures)
        run_r("render_targets.R", args.publication_figures)
    if include_shared:
        command = [sys.executable, str(HELPERS / "render_shared_peak_genes.py"), "--pipeline-root", str(RUN_ROOT)]
        if args.publication_figures:
            command.append("--publication-figures")
        subprocess.run(command, check=True)
    if args.publication_figures and not args.shared_peaks_only:
        print(f"[DONE] Regulatory-target publication figures: {PUBLICATION_FIGURES}")
    if not args.shared_peaks_only:
        print(f"[DONE] Regulatory-target figures: {VISUALS}")
    if include_shared:
        print(f"[DONE] Shared-peak-associated DEG figures: {VISUALS / 'shared_peaks_visuals'}")


if __name__ == "__main__":
    main()
