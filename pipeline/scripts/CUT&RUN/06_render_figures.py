#!/usr/bin/env python3
"""Render and publish the final CUT&RUN figures."""

#Before you run this script, please replace directory placeholders with your own directory

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import CUTRUN_ROOT, FACTORS, REFERENCE_DIR, RUN_ROOT


TAG = "q5e2_fe3_min2of2"
HELPERS = Path(__file__).resolve().parent / "figures_rendering"
FIGURE_DATA = CUTRUN_ROOT / "data" / "figure_inputs"
VISUALS = CUTRUN_ROOT / "visuals"
PROMOTER_VISUALS = VISUALS / "promoters_visuals"
GENE_VISUALS = VISUALS / "all_genes_visuals"
PEAK_VISUALS = VISUALS / "all_peaks_visuals"
TRACK_ROOT = CUTRUN_ROOT / "03_bigwig" / "IGV_representation"
TRACKS = {factor: TRACK_ROOT / f"{factor}_mean.bw" for factor in FACTORS}
IGG = TRACK_ROOT / "IgG_mean.bw"


def load_module(name: str):
    path = HELPERS / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def executable(name: str, fallback: str) -> str:
    path = shutil.which(name)
    if path:
        return path
    if Path(fallback).is_file():
        return fallback
    raise FileNotFoundError(f"Required executable not found: {name}")


def prepare_inputs() -> None:
    FIGURE_DATA.mkdir(parents=True, exist_ok=True)
    counts = CUTRUN_ROOT / "data" / "Venn_Peaks_counts.tsv"
    if not counts.is_file():
        raise FileNotFoundError(counts)
    shutil.copy2(counts, FIGURE_DATA / "Venn_Peaks_counts.tsv")


def enable_text_free_rendering() -> None:
    """Suppress all Matplotlib text and append `_noTexts` during figure export."""
    from matplotlib.figure import Figure
    from matplotlib.text import Text

    original_savefig = Figure.savefig

    def savefig_without_text(self, fname, *args, **kwargs):
        path = Path(fname)
        if path.suffix.lower() == ".png":
            fname = path.with_name(f"{path.stem}_noTexts{path.suffix}")
        return original_savefig(self, fname, *args, **kwargs)

    Figure.savefig = savefig_without_text
    Text.draw = lambda self, renderer: None


def configure_primary():
    module = load_module("render_profiles.py")
    module.PROJECT = HELPERS
    module.RUN = RUN_ROOT
    module.CUTRUN = CUTRUN_ROOT
    module.DATA = FIGURE_DATA
    module.VISUALS = PROMOTER_VISUALS
    module.PROMOTER = FIGURE_DATA / "promoter_gene_venn"
    module.TAG = TAG
    module.GENE_BED = CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6"
    module.TRACKS = TRACKS
    module.IGG = IGG
    module.COMPUTE_MATRIX = Path(executable("computeMatrix", "/opt/anaconda3/envs/cutrun_env/bin/computeMatrix"))
    module.RSCRIPT = executable("Rscript", "Rscript")
    return module


def configure_additional():
    module = load_module("render_peak_profiles.py")
    module.PROJECT = HELPERS
    module.RUN = RUN_ROOT
    module.CUTRUN = CUTRUN_ROOT
    module.DATA = FIGURE_DATA
    module.VISUALS = VISUALS
    module.TAG = TAG
    module.GTF = REFERENCE_DIR / "gencode.vM25.annotation.gtf"
    module.GENE_BED = CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6"
    module.TRACKS = TRACKS
    module.PEAK_LOCI = CUTRUN_ROOT / "data" / "PeakLoci" / "Venn_Peaks_loci.tsv"
    module.PROMOTERS = CUTRUN_ROOT / "data" / "Promoters_M25_TSSplusminus1kb.bed"
    module.OUT = FIGURE_DATA / "additional_panels"
    module.COMPUTE_MATRIX = Path(executable("computeMatrix", "/opt/anaconda3/envs/cutrun_env/bin/computeMatrix"))
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--branch", choices=("all", "promoters", "all-genes", "all-peaks"), default="all")
    parser.add_argument("--publication-figures", action="store_true", help="Also render text-free PNGs within each visual branch")
    parser.add_argument("--no-text", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.no_text:
        enable_text_free_rendering()
    if args.dry_run:
        print("[DRY-RUN] Would render the registered CUT&RUN figures.")
        return
    prepare_inputs()
    selected = {
        "promoters": args.branch in ("all", "promoters"),
        "all-genes": args.branch in ("all", "all-genes"),
        "all-peaks": args.branch in ("all", "all-peaks"),
    }
    branch_dirs = {"promoters": PROMOTER_VISUALS, "all-genes": GENE_VISUALS, "all-peaks": PEAK_VISUALS}
    output_dirs = {
        key: path / "publication_figures" if args.no_text else path
        for key, path in branch_dirs.items()
    }
    for key, enabled in selected.items():
        if enabled:
            output_dirs[key].mkdir(parents=True, exist_ok=True)
    required = [REFERENCE_DIR / "gencode.vM25.annotation.gtf", CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6", IGG, *TRACKS.values()]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing CUT&RUN figure inputs:\n" + "\n".join(missing))
    if selected["all-genes"] or selected["all-peaks"]:
        render_peak_associations(args.no_text, output_dirs["all-genes"], output_dirs["all-peaks"],
                                 selected["all-genes"], selected["all-peaks"])
    primary = configure_primary()
    primary.NO_TEXT_VISUALS = None
    primary.TEXT_FREE = args.no_text
    additional = configure_additional()
    additional.NO_TEXT_VISUALS = None
    additional.TEXT_FREE = args.no_text
    if selected["promoters"]:
        primary.render_promoter_branch(output_dirs["promoters"])
        additional.VISUALS = output_dirs["promoters"]
        additional.render_promoter_gene_profile()
    if selected["all-genes"]:
        primary.render_peak_associated_branch(output_dirs["all-genes"])
        additional.VISUALS = output_dirs["all-genes"]
        additional.render_peak_associated_gene_profile()
    if selected["all-peaks"]:
        additional.VISUALS = output_dirs["all-peaks"]
        additional.render_peak_associated_gene_profiles()
    if args.publication_figures and not args.no_text:
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--no-text", "--branch", args.branch], check=True)
    print(f"[DONE] CUT&RUN figures: {VISUALS}")


def render_peak_associations(no_text: bool, gene_visuals: Path, peak_visuals: Path,
                             render_genes: bool, render_peaks: bool) -> None:
    module = load_module("render_peak_associated_genes.py")
    module.CUTRUN = CUTRUN_ROOT
    module.DATA = FIGURE_DATA
    module.GTF = REFERENCE_DIR / "gencode.vM25.annotation.gtf"
    module.GENE_VISUALS = gene_visuals
    module.PEAK_VISUALS = peak_visuals
    module.TEXT_FREE = no_text
    module.RENDER_GENES = render_genes
    module.RENDER_PEAKS = render_peaks
    module.main()


if __name__ == "__main__":
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[variable] = "1"
    os.environ.setdefault("MPLCONFIGDIR", str(CUTRUN_ROOT / "data" / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    main()
