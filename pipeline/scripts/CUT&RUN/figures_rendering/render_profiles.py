#!/usr/bin/env python3
"""Render promoter-overlap, promoter-profile, and RPKM figures."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parent
RUN = PROJECT.parents[2]
CUTRUN = RUN / "cutrun_work"
DATA = CUTRUN / "data" / "figure_inputs"
VISUALS = CUTRUN / "visuals"
NO_TEXT_VISUALS: Path | None = None
TEXT_FREE = False
GENE_MEMBERSHIP = "independent_genes"
PROMOTER = DATA / "promoter_gene_venn"
TAG = "q5e2_fe3_min2of2"
FACTORS = ("MCM3", "NONO", "PSPC1")
COLORS = {"MCM3": "#A80F14", "NONO": "#F39C12", "PSPC1": "#1F4AA8", "IgG": "#6B7280"}
GENE_BED = CUTRUN / "data" / "GeneBodies_M25.bed6"
TRACKS = {
    factor: CUTRUN / "03_bigwig" / "IGV_representation" / f"{factor}_mean.bw"
    for factor in FACTORS
}
IGG = CUTRUN / "03_bigwig" / "IGV_representation" / "IgG_mean.bw"
COMPUTE_MATRIX = Path("computeMatrix")
RSCRIPT = "Rscript"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_genes(path: Path) -> set[str]:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return set(table["gene"].str.strip()) - {""}


def promoter_sets() -> dict[str, set[str]]:
    return {
        factor: read_genes(PROMOTER / f"Venn_PromoterGenes_{factor}.tsv")
        for factor in FACTORS
    }


def venn_regions(sets: dict[str, set[str]]) -> dict[str, int]:
    a, b, c = (sets[factor] for factor in FACTORS)
    return {
        "111_MCM3_NONO_PSPC1": len(a & b & c),
        "110_MCM3_NONO": len((a & b) - c),
        "101_MCM3_PSPC1": len((a & c) - b),
        "100_only_MCM3": len(a - b - c),
    }


def render_promoter_venn(sets: dict[str, set[str]], visual_dir: Path) -> None:
    module = load_module(
        PROJECT / "promoter_venn.py",
        "promoter_venn",
    )
    module.PLOT_TITLE = "CUT&RUN promoter-gene overlap"
    module.plot_venn(
        sets,
        visual_dir / "Venn_PromoterGenes.png",
        show_numbers=True,
    )
    if NO_TEXT_VISUALS is not None:
        module.plot_venn(
            sets,
            NO_TEXT_VISUALS / "Venn_PromoterGenes_noTexts.png",
            show_numbers=False,
        )


def write_region_table(sets: dict[str, set[str]], path: Path) -> None:
    a, b, c = (sets[factor] for factor in FACTORS)
    regions = {
        "111_MCM3_NONO_PSPC1": a & b & c,
        "110_MCM3_NONO": (a & b) - c,
        "101_MCM3_PSPC1": (a & c) - b,
        "100_only_MCM3": a - b - c,
        "010_only_NONO": b - a - c,
        "001_only_PSPC1": c - a - b,
        "011_NONO_PSPC1": (b & c) - a,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{"region": region, "gene": gene} for region, genes in regions.items() for gene in sorted(genes)],
        columns=["region", "gene"],
    ).to_csv(path, sep="\t", index=False)


def peak_associated_inputs() -> tuple[dict[str, set[str]], Path, Path]:
    root = DATA / "protein_coding_peak_associations"
    sets = {
        factor: set(pd.read_csv(root / f"PeakAssociatedGenes_{factor}.tsv", sep="\t", dtype=str)["gene_id"].dropna())
        for factor in FACTORS
    }
    input_dir = DATA / "peak_associated_gene_profiles" / GENE_MEMBERSHIP
    input_dir.mkdir(parents=True, exist_ok=True)
    region_table = input_dir / "Venn_PeakAssociatedGenes_regions.tsv"
    if GENE_MEMBERSHIP == "shared_locus":
        groups = pd.read_csv(root / "LocusDefinedGeneGroups.tsv", sep="\t")
        groups[["region", "gene_id"]].rename(columns={"gene_id": "gene"}).to_csv(region_table, sep="\t", index=False)
    else:
        write_region_table(sets, region_table)
    gene_bed = input_dir / "GeneBodies_M25_byID.bed6"
    rows = []
    with (RUN / "reference" / "gencode.vM25.annotation.gtf").open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = re.search(r'gene_id "([^"]+)"', fields[8])
            if match is None:
                continue
            gene_id = match.group(1).split(".")[0]
            rows.append((fields[0], int(fields[3]) - 1, int(fields[4]), gene_id, 0, fields[6]))
    pd.DataFrame(rows).drop_duplicates(3).to_csv(gene_bed, sep="\t", header=False, index=False)
    return sets, region_table, gene_bed


def render_gene_profiles(
    sets: dict[str, set[str]], visual_dir: Path, region_table: Path,
    gene_bed: Path, data_namespace: str, gene_label: str, locus_defined: bool = False,
) -> None:
    module = load_module(
        PROJECT / "promoter_profiles.py",
        "promoter_profiles",
    )
    counts = pd.read_csv(region_table, sep="\t").groupby("region").gene.nunique().to_dict()
    module.OUTDIR = visual_dir
    module.TMPDIR = DATA / data_namespace / "intermediate"
    module.MATRIX_TMPDIR = module.TMPDIR / "matrices"
    module.REGION_TABLE = region_table
    module.GENE_BED = gene_bed
    module.BIGWIGS = TRACKS
    module.COMPUTE_MATRIX = COMPUTE_MATRIX
    module.DISPLAY_SCALE = 1.0
    module.SIGNAL_LABEL = "Coverage"
    module.REGIONS = {
        "111": {
            "data_stem": "Profile_MCM3_NONO_PSPC1",
            "table_region": "111_MCM3_NONO_PSPC1",
            "label": f"MCM3-PSPC1-NONO co-bound {gene_label}",
            "display_title": "MCM3-PSPC1-NONO Co-localized\ngenes",
            "tracks": FACTORS,
            "expected_n": counts["111_MCM3_NONO_PSPC1"],
        },
        "110": {
            "data_stem": "Profile_MCM3_NONO",
            "table_region": "110_MCM3_NONO",
            "label": f"MCM3-NONO-only {gene_label}",
            "display_title": f"MCM3-NONO-only\n{gene_label}",
            "tracks": ("MCM3", "NONO"),
            "expected_n": counts["110_MCM3_NONO"],
        },
        "101": {
            "data_stem": "Profile_MCM3_PSPC1",
            "table_region": "101_MCM3_PSPC1",
            "label": f"MCM3-PSPC1-only {gene_label}",
            "display_title": f"MCM3-PSPC1-only\n{gene_label}",
            "tracks": ("MCM3", "PSPC1"),
            "expected_n": counts["101_MCM3_PSPC1"],
        },
    }
    if locus_defined:
        for key, label in (("111", "MCM3-PSPC1-NONO shared loci"),
                           ("110", "MCM3-NONO-only loci"), ("101", "MCM3-PSPC1-only loci")):
            module.REGIONS[key]["label"] = label + ": associated genes"
            module.REGIONS[key]["display_title"] = label + "\nAssociated genes"
    elif data_namespace.startswith("peak_associated_gene_profiles"):
        for key, label in (("111", "MCM3-PSPC1-NONO"), ("110", "MCM3-NONO-only"), ("101", "MCM3-PSPC1-only")):
            module.REGIONS[key]["label"] = label + " peak-associated genes"
            module.REGIONS[key]["display_title"] = label + "\nPeak-associated genes"
    module.main(reuse_matrices=TEXT_FREE)
    source = visual_dir / "Profile_PromoterGenes_summary.tsv"
    if source.is_file():
        source.replace(DATA / data_namespace / "Profile_Genes_summary.tsv")


def render_rpkm_profiles(
    sets: dict[str, set[str]], visual_dir: Path, region_table: Path,
    gene_bed: Path, data_namespace: str, locus_defined: bool = False,
) -> None:
    module = load_module(PROJECT / "rpkm_profiles.py", "rpkm_profiles")
    counts = pd.read_csv(region_table, sep="\t").groupby("region").gene.nunique().to_dict()
    outdir = DATA / data_namespace / "rpkm_intermediate"
    module.CUTRUN_ROOT = CUTRUN
    module.OUTDIR = outdir
    module.TMPDIR = outdir / "regions"
    module.MATRIX_TMPDIR = outdir / "matrices"
    module.REGION_TABLE = region_table
    module.GENE_BED = gene_bed
    module.BIGWIGS = TRACKS
    module.COMPUTE_MATRIX = COMPUTE_MATRIX
    module.DISPLAY_SCALE = 1.0
    module.SIGNAL_LABEL = "Log2(RPKM)"
    module.OUTPUT = outdir / "MCM3_RPKM.png"
    module.OUTPUT_NONO = outdir / "NONO_RPKM.png"
    module.OUTPUT_PSPC1 = outdir / "PSPC1_RPKM.png"
    module.SUMMARY = outdir / "RPKM_profiles_summary.tsv"
    module.REGION_SPECS = (
        ("111", "111_MCM3_NONO_PSPC1", "MCM3-PSPC1-NONO co-binding", counts["111_MCM3_NONO_PSPC1"]),
        ("110", "110_MCM3_NONO", "MCM3-NONO co-binding", counts["110_MCM3_NONO"]),
        ("101", "101_MCM3_PSPC1", "MCM3-PSPC1 co-binding", counts["101_MCM3_PSPC1"]),
        ("100", "100_only_MCM3", "MCM3-specific binding", counts["100_only_MCM3"]),
    )
    if locus_defined:
        labels = ("Genes at MCM3-PSPC1-NONO loci", "Genes at MCM3-NONO-only loci",
                  "Genes at MCM3-PSPC1-only loci", "Genes at MCM3-only loci")
        module.REGION_SPECS = tuple((key, region, label, n) for (key, region, _, n), label
                                    in zip(module.REGION_SPECS, labels))
    elif data_namespace.startswith("peak_associated_gene_profiles"):
        labels = ("MCM3-PSPC1-NONO-associated genes", "MCM3-NONO-only-associated genes",
                  "MCM3-PSPC1-only-associated genes", "MCM3-only-associated genes")
        module.REGION_SPECS = tuple((key, region, label, n) for (key, region, _, n), label
                                    in zip(module.REGION_SPECS, labels))
    module.FACTOR_SPECS = (
        ("MCM3", 0, module.OUTPUT),
        ("NONO", 1, module.OUTPUT_NONO),
        ("PSPC1", 2, module.OUTPUT_PSPC1),
    )
    def fitted_r_plot(values_tsv: Path, output_png: Path) -> None:
        values = pd.read_csv(values_tsv, sep="\t").groupby("group", observed=True)["value"]
        whiskers = []
        for _, group_values in values:
            finite = group_values[np.isfinite(group_values)].to_numpy(dtype=float)
            q1, q3 = np.quantile(finite, [0.25, 0.75])
            whiskers.append(float(finite[finite <= q3 + 1.5 * (q3 - q1)].max()))
        target = max(1.0, max(whiskers) * 1.08)
        magnitude = 10.0 ** math.floor(math.log10(target))
        module.Y_MAX = next(value * magnitude for value in (1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10) if value * magnitude >= target)
        labels = [f"{label} (n={n:,})" for _, _, label, n in module.REGION_SPECS]
        levels = [label for _, _, label, _ in module.REGION_SPECS]
        r_script = r'''
# Keep implicit layout devices off disk; explicit PNG/PDF devices are unchanged.
local({
  .plot_devices_before <- grDevices::dev.list()
  .plot_device_option <- getOption("device")
  .rplots_path <- file.path(getwd(), "Rplots.pdf")
  on.exit({
    for (.plot_device in setdiff(grDevices::dev.list(), .plot_devices_before)) {
      try(grDevices::dev.off(.plot_device), silent = TRUE)
    }
    options(device = .plot_device_option)
    # Remove only the default-device artifact in this execution directory.
    unlink(c(.rplots_path, file.path(dirname(.rplots_path), "._Rplots.pdf")))
  }, add = TRUE)
  options(device = function(...) grDevices::pdf(file = NULL))
  grDevices::pdf(file = NULL)
args <- commandArgs(trailingOnly = TRUE)
df <- read.delim(args[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
levels <- strsplit(args[4], "\\|", fixed = FALSE)[[1]]
labels <- strsplit(args[5], "\\|", fixed = FALSE)[[1]]
df$group <- factor(df$group, levels = levels)
cols <- c("#cf111f", "#f2a010", "#4f7db1", "#777777")
no_text <- as.logical(args[7])
png(args[2], width = 2800, height = 2119, res = 300, bg = "white")
if (no_text) par(mar = c(.5, .5, .5, .5)) else par(mar = c(13, 10, 2, 5) + 0.1, cex.axis = 1.25, cex.lab = 1.45, font.axis = 2, font.lab = 2)
boxplot(value ~ group, data = df, col = cols, border = "#222222", lwd = 2.0,
        ylab = if (no_text) "" else args[3], xlab = "", xaxt = "n", yaxt = "n", outline = FALSE,
        ylim = c(0, as.numeric(args[6])), whisklty = 1, staplelty = 1)
if (!no_text) {
  axis(1, at = seq_along(levels), labels = FALSE, tick = FALSE)
  text(x = seq_along(levels), y = rep(par("usr")[3] - 0.39, length(levels)), labels = labels,
       srt = 35, xpd = NA, adj = 1, cex = 0.74, font = 2)
  axis(2, at = seq(0, as.numeric(args[6]), by = 1), labels = seq(0, as.numeric(args[6]), by = 1), las = 1, cex.axis = 1.25, font = 2)
}
dev.off()
})
'''
        process = subprocess.run(
            [RSCRIPT, "-", str(values_tsv), str(output_png), module.SIGNAL_LABEL, "|".join(levels), "|".join(labels), str(module.Y_MAX), str(TEXT_FREE).upper()],
            input=r_script,
            text=True,
            capture_output=True,
            check=False,
        )
        if process.returncode:
            raise RuntimeError(process.stderr or process.stdout)

    module.r_plot = fitted_r_plot
    module.main()
    for name in ("MCM3_RPKM.png", "NONO_RPKM.png", "PSPC1_RPKM.png"):
        output_name = f"{Path(name).stem}_noTexts.png" if TEXT_FREE else name
        (outdir / name).replace(visual_dir / output_name)
    if module.SUMMARY.is_file():
        module.SUMMARY.replace(DATA / data_namespace / module.SUMMARY.name)


def render_promoter_branch(visual_dir: Path) -> None:
    visual_dir.mkdir(parents=True, exist_ok=True)
    sets = promoter_sets()
    render_promoter_venn(sets, visual_dir)
    render_gene_profiles(
        sets, visual_dir, PROMOTER / "Venn_PromoterGenes_regions.tsv", GENE_BED,
        "promoter_profiles", "promoter genes",
    )
    render_rpkm_profiles(
        sets, visual_dir, PROMOTER / "Venn_PromoterGenes_regions.tsv", GENE_BED,
        "promoter_profiles",
    )


def render_peak_associated_branch(visual_dir: Path) -> None:
    visual_dir.mkdir(parents=True, exist_ok=True)
    sets, region_table, gene_bed = peak_associated_inputs()
    render_gene_profiles(
        sets, visual_dir, region_table, gene_bed,
        f"peak_associated_gene_profiles/{GENE_MEMBERSHIP}", "peak-associated genes",
        locus_defined=GENE_MEMBERSHIP == "shared_locus",
    )
    render_rpkm_profiles(
        sets, visual_dir, region_table, gene_bed,
        f"peak_associated_gene_profiles/{GENE_MEMBERSHIP}",
        locus_defined=GENE_MEMBERSHIP == "shared_locus",
    )


def main() -> None:
    required = [GENE_BED, *TRACKS.values(), IGG]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing rendering inputs:\n" + "\n".join(missing))
    for directory in (VISUALS, DATA / "promoter_profiles", DATA / "peak_associated_gene_profiles"):
        directory.mkdir(parents=True, exist_ok=True)
    render_promoter_branch(VISUALS)
    print(f"[DONE] CUT&RUN promoter and RPKM panels: {VISUALS}")


if __name__ == "__main__":
    main()
