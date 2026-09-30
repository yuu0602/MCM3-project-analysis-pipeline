#!/usr/bin/env python3
"""Render promoter, shared-locus gene, and shared-locus peak regulatory analyses."""

#Before you run this script, please replace directory placeholders with your own directory

from __future__ import annotations

import argparse
import hashlib
import json
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
BRANCHES = {
    "promoters": VISUALS / "promoters_visuals",
    "all_genes": VISUALS / "all_genes_visuals",
    "all_peaks": VISUALS / "all_peaks_visuals",
}


def prepare_shared_peak_genes(output: Path) -> None:
    """Assign the three-factor peak intersection using Step 06's existing links."""
    import pandas as pd

    loci_path = CUTRUN_ROOT / "data/PeakLoci/Venn_Peaks_loci.tsv"
    links_path = CUTRUN_ROOT / "data/figure_inputs/protein_coding_peak_associations/PeakGeneAssignments.tsv"
    loci = pd.read_csv(loci_path, sep="\t")
    shared = loci.loc[loci.region_mask.eq(7)].copy()
    if not shared.locus_id.is_unique or shared.empty:
        raise ValueError("Expected unique three-factor shared peak loci from Step 06")
    links = pd.read_csv(links_path, sep="\t")
    links = links.loc[links.peak_id.isin(shared.locus_id)].copy()
    if not links.gene_type.eq("protein_coding").all():
        raise ValueError("Step 06 assignments must contain only protein-coding genes")
    links["gene_id"] = links.gene_id.str.replace(r"\..*$", "", regex=True)
    links = links.drop_duplicates(["peak_id", "gene_id"])
    assigned = shared.loc[shared.locus_id.isin(links.peak_id)]
    genes = links.groupby("gene_id", as_index=False).agg(
        gene=("gene", "first"), n_shared_peaks=("peak_id", "nunique")
    )
    if genes.empty:
        raise ValueError("No protein-coding genes assigned to the shared peak loci")
    data = output / "data"
    data.mkdir(parents=True, exist_ok=True)
    shared.to_csv(data / "SharedPeakLoci.tsv", sep="\t", index=False)
    assigned.to_csv(data / "SharedPeaks.tsv", sep="\t", index=False)
    assigned[["chrom", "start", "end"]].to_csv(data / "SharedPeaks.bed", sep="\t", header=False, index=False)
    links.to_csv(data / "SharedPeakGeneAssignments.tsv", sep="\t", index=False)
    genes.to_csv(data / "SharedPeakAssociatedGenes.tsv", sep="\t", index=False)
    parameters = {
        "n_three_factor_shared_peak_loci": len(shared),
        "n_shared_peaks_with_protein_coding_assignments": len(assigned),
        "n_shared_peak_associated_genes": len(genes),
        "peak_rule": "Three-factor intersection (region_mask=7) of canonical Step 06 peak loci",
        "gene_assignment": "All promoters at TSS +/-1000 bp overlapping >=250 bp; nearest-TSS fallback only if no promoter qualifies; protein-coding filter after assignment",
        "regulatory_rule": "The same shared-peak-associated gene universe intersected with each knockdown DEG set by Ensembl ID",
        "group_meaning": "A-D describe knockdown DEG overlap, not factor-exclusive binding",
        "inputs_sha256": {str(p.relative_to(RUN_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (loci_path, links_path)},
    }
    (data / "SharedPeakParameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    print(f"[SHARED PEAKS] {len(shared)} loci; {len(assigned)} gene-assigned loci; {len(genes)} protein-coding genes", flush=True)


def run_r(script: str, mode: str, output: Path, publication_figures: bool) -> None:
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
    command = [str(executable), str(HELPERS / script), str(RUN_ROOT), mode, str(output)]
    if publication_figures:
        publication = output / "publication_figures"
        publication.mkdir(parents=True, exist_ok=True)
        command.append(str(publication))
    subprocess.run(command, check=True, env=environment)


def prepare_peak_targets(output: Path, publication_figures: bool) -> None:
    """Count each shared locus once per KD when any assigned coding gene is a DEG."""
    import pandas as pd

    prepare_shared_peak_genes(output)
    data = output / "data"
    loci = pd.read_csv(data / "SharedPeakLoci.tsv", sep="\t")
    evidence = pd.read_csv(data / "SharedPeakGeneAssignments.tsv", sep="\t")
    table = loci.rename(columns={"locus_id": "peak_id"}).copy()
    table["n_assigned_genes"] = table.peak_id.map(evidence.groupby("peak_id").gene_id.nunique()).fillna(0).astype(int)
    factors = ("MCM3", "NONO", "PSPC1")
    sets, summaries, hashes = {}, [], {}
    for factor in factors:
        gene_sets = {}
        for direction in ("UP", "DOWN"):
            path = RUN_ROOT / f"deg_work/data/VennDiagram_{direction}_{factor}_genes.tsv"
            gene_sets[direction] = set(pd.read_csv(path, sep="\t").gene_id.str.replace(r"\..*$", "", regex=True))
            hashes[str(path.relative_to(RUN_ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        if gene_sets["UP"] & gene_sets["DOWN"]:
            raise ValueError(f"Genes occur in both DEG directions for {factor}")
        evidence[factor] = evidence.gene_id.map(
            {**dict.fromkeys(gene_sets["UP"], "Up"), **dict.fromkeys(gene_sets["DOWN"], "Down")}
        ).fillna("No DEG")
        for direction in ("Up", "Down"):
            linked = evidence.loc[evidence[factor].eq(direction)].groupby("peak_id").gene_id.nunique()
            table[f"{factor}_n_{direction.lower()}_genes"] = table.peak_id.map(linked).fillna(0).astype(int)
        up = table[f"{factor}_n_up_genes"].gt(0)
        down = table[f"{factor}_n_down_genes"].gt(0)
        table[factor] = "No DEG"
        table.loc[up & ~down, factor] = "Up"
        table.loc[down & ~up, factor] = "Down"
        table.loc[up & down, factor] = "Mix"
        sets[factor] = set(table.loc[up | down, "peak_id"])
        table.loc[up | down].to_csv(data / f"Venn_target_{factor}_peaks.tsv", sep="\t", index=False)
        for direction in ("Up", "Down", "Mix", "No DEG"):
            summaries.append({"factor": factor, "direction": direction,
                              "n_peaks": int(table[factor].eq(direction).sum()), "total_shared_peaks": len(table)})
    evidence.to_csv(data / "PeakGeneDEGEvidence.tsv", sep="\t", index=False)
    table.to_csv(data / "Peak_DEG_directions.tsv", sep="\t", index=False)
    direction_counts = pd.DataFrame(summaries)
    direction_counts.to_csv(data / "Peak_DEG_direction_counts.tsv", sep="\t", index=False)
    pie_counts = direction_counts.loc[direction_counts.direction.ne("No DEG")].copy()
    pie_counts["DEG_linked_peaks"] = pie_counts.factor.map({f:len(sets[f]) for f in factors})
    pie_counts.to_csv(data / "Pie_DEGAssociation_counts.tsv", sep="\t", index=False)
    a, b, c = (sets[f] for f in factors)
    regions = {"100": a-b-c, "010": b-a-c, "001": c-a-b,
               "110": (a & b)-c, "101": (a & c)-b, "011": (b & c)-a, "111": a & b & c}
    names = {"100":"MCM3_only", "010":"NONO_only", "001":"PSPC1_only",
             "110":"MCM3_NONO_only", "101":"MCM3_PSPC1_only", "011":"NONO_PSPC1_only", "111":"MCM3_NONO_PSPC1"}
    pd.DataFrame([{"region":names[mask], "peak_id":peak} for mask, members in regions.items()
                  for peak in sorted(members)], columns=["region", "peak_id"]).to_csv(
        data / "Venn_target_peaks.tsv", sep="\t", index=False)
    pd.DataFrame([{"region":key,"n_peaks":len(value)} for key,value in regions.items()]).to_csv(
        data / "Venn_target_counts.tsv", sep="\t", index=False)
    metadata = json.loads((data / "SharedPeakParameters.json").read_text())
    metadata.update({
        "analysis_mode":"all_peaks", "counting_unit":"unique genomic peak locus",
        "regulatory_rule":"At least one assigned protein-coding gene is a DEG for that KD; each peak counted once per KD",
        "group_meaning":"A-D are intersections of KD DEG-linked peak sets, not exclusive CUT&RUN binding; different assigned genes can support different KDs",
        "direction_rule":"Up or Down when all assigned DEG genes agree within a KD; Mix when both directions occur. Group Up/Down requires agreement across all indicated KDs; otherwise Mix",
        "no_DEG_rule":"No assigned significant DEG; includes peaks lacking coding assignments or genes not tested, and does not imply no biological effect",
        "RNAseq_threshold":"BH-adjusted p <=0.05 and absolute log2FC >=0.28; unchanged directional DEG tables",
        "profile_rule":"Group A-D: peak-center aligned +/-3000 bp, 25-bp bins, unstranded; one row per locus; arithmetic mean, missing coverage zero",
        "MCM3_target_profile_rule":"Unique MCM3 Up/Down DEG genes assigned to Up-only/Down-only shared peaks; exclude mixed-direction peak assignments; strand-aware GENCODE M25 gene bodies scaled to 3000 bp with 3000-bp flanks and 25-bp bins. One row per gene, equal gene weights. Peak membership and other figures unchanged",
        "normalization":"Existing yeast-normalized IGV_representation mean bigWigs; no new scaling or normalization",
        "direction_pie_denominator":"DEG-linked shared peaks only; No DEG excluded from pies, retained in complete source tables",
        "n_DEG_linked_peaks":{f:len(sets[f]) for f in factors},
    })
    metadata["inputs_sha256"].update(hashes)
    (data / "AnalysisParameters.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (data / "AnalysisParameters.tsv").write_text("\n".join(f"{k}\t{json.dumps(v) if isinstance(v,dict) else v}" for k,v in metadata.items()) + "\n")
    for no_text in (False, True) if publication_figures else (False,):
        dest = output / "publication_figures" if no_text else output
        dest.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(HELPERS / "venn.py"), "--out", str(dest / ("Venn_target_noTexts.png" if no_text else "Venn_target.png")),
                   "--title", "DEG-linked shared peaks", "--subtitle", "Overlap denotes KD DEG association, not differential binding"]
        for letter,factor in zip("abc", factors):
            command += [f"--{letter}-name",factor,f"--{letter}-total",str(len(sets[factor]))]
        for mask,members in regions.items():
            command += [f"--n{mask}",str(len(members))]
        if no_text:
            command.append("--hide-numbers")
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(data / ".matplotlib")
        subprocess.run(command, check=True, env=environment)
    render_peak_direction_panels(output, table, publication_figures)
    print("[PEAK TARGETS]", {f:len(sets[f]) for f in factors}, "A-D:",
          [len(regions[key]) for key in ("111","101","011","110")], flush=True)


def render_peak_direction_panels(output: Path, shared, publication_figures: bool) -> None:
    """Use the reference circle design with exclusively peak-counted sets."""
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = output / "data"
    loci = pd.read_csv(CUTRUN_ROOT / "data/PeakLoci/Venn_Peaks_loci.tsv", sep="\t")
    links = pd.read_csv(CUTRUN_ROOT / "data/figure_inputs/protein_coding_peak_associations/PeakGeneAssignments.tsv", sep="\t")
    links["gene_id"] = links.gene_id.str.replace(r"\..*$", "", regex=True)
    shared_ids = set(shared.peak_id)
    records, backgrounds = [], []
    factors = ("MCM3", "NONO", "PSPC1")
    for i, (factor, color, bit) in enumerate(zip(factors, ("#A80F14", "#F39C12", "#1F4AA8"), (1,2,4))):
        factor_loci = loci.loc[(loci.region_mask & bit).ne(0)].copy()
        factor_links = links.loc[links.peak_id.isin(factor_loci.locus_id)]
        directional = {}
        for name in ("UP", "DOWN"):
            ids = set(pd.read_csv(RUN_ROOT / f"deg_work/data/VennDiagram_{name}_{factor}_genes.tsv", sep="\t").gene_id.str.replace(r"\..*$", "", regex=True))
            directional[name] = set(factor_links.loc[factor_links.gene_id.isin(ids), "peak_id"])
        mix = directional["UP"] & directional["DOWN"]
        up, down = directional["UP"] - mix, directional["DOWN"] - mix
        factor_loci["factor"] = factor
        factor_loci["DEG_direction"] = factor_loci.locus_id.map({**dict.fromkeys(up,"Up"), **dict.fromkeys(down,"Down"), **dict.fromkeys(mix,"Mix")}).fillna("No DEG")
        factor_loci["three_factor_shared"] = factor_loci.locus_id.isin(shared_ids)
        backgrounds.append(factor_loci)
        overlap_up, overlap_down = shared_ids & up, shared_ids & down
        assert overlap_up == set(shared.loc[shared[factor].eq("Up"), "peak_id"])
        assert overlap_down == set(shared.loc[shared[factor].eq("Down"), "peak_id"])
        numbers = {"n_center":len(shared_ids), "n_up":len(up), "n_down":len(down),
                   "n_only_center":len(shared_ids - up - down), "n_only_up":len(up - shared_ids),
                   "n_only_down":len(down - shared_ids), "n_overlap_up":len(overlap_up), "n_overlap_down":len(overlap_down)}
        detail = {"shared_No_DEG":int(shared[factor].eq("No DEG").sum()), "shared_Mix":len(shared_ids & mix)}
        records.append({"factor":factor, **numbers, **detail})
        for no_text in (False, True) if publication_figures else (False,):
            dest = output / "publication_figures" if no_text else output
            name = f"peak_vs_DEG_direction_{factor}" + ("_noTexts" if no_text else "") + ".png"
            opts = {"out":dest/name, "panel-letter":"ABC"[i], "factor":factor,
                    "center-label":f"{factor} CUT&RUN shared peaks", "center-col":color,
                    "left-label":"Up-only-linked peaks", "right-label":"Down-only-linked peaks",
                    "col-up":"#D87070", "col-down":"#6FA37A", "col-overlap":"#BFBFBF",
                    "stroke-col":"#222222", "stroke-lwd":2, "r-side":.170, "r-center":.205,
                    "cx":.5,"cy":.53,"lx":.305,"ly":.53,"rx":.695,"ry":.53,
                    "width-px":3200,"height-px":2400,"dpi":300,"fs-panel-letter":20,
                    "fs-title":24,"fs-subtitle":20,"fs-num-center":22,"fs-num-overlap":17,
                    "fs-bottom-lab":20,"fs-bottom-n":18}
            opts.update({k.replace("_","-"):v for k,v in numbers.items()})
            cmd = [sys.executable,str(HELPERS / "direction_panel.py")]
            for key,value in opts.items():
                cmd.extend([f"--{key}",str(value)])
            if no_text:
                cmd.append("--no-text")
            subprocess.run(cmd, check=True, env={**os.environ,"MPLCONFIGDIR":str(data / ".matplotlib")})
    pd.DataFrame(records).to_csv(data / "peak_vs_DEG_direction_counts.tsv", sep="\t", index=False)
    pd.concat(backgrounds, ignore_index=True).to_csv(data / "peak_vs_DEG_direction_background.tsv", sep="\t", index=False)
    (data / "peak_vs_DEG_direction_parameters.json").write_text(json.dumps({
        "unit":"unique genomic locus",
        "center":"All three-factor shared loci",
        "side_sets":"All retained loci of the indicated factor linked exclusively to Up or exclusively to Down DEGs; both-direction links are Mix",
        "center_only":"Shared loci without an Up-only or Down-only link: No DEGs plus Mix; breakdown retained in the counts table",
        "interpretation":"A schematic peak-set diagram, not an intersection between peak counts and gene counts; no promoter-only restriction",
    },indent=2)+"\n")
    for no_text in (False, True) if publication_figures else (False,):
        dest = output / "publication_figures" if no_text else output
        suffix = "_noTexts" if no_text else ""
        fig, axes = plt.subplots(3,1,figsize=(2200/300,5200/300),dpi=300)
        for ax,factor in zip(axes,factors):
            ax.imshow(plt.imread(dest / f"peak_vs_DEG_direction_{factor}{suffix}.png"))
            ax.set_axis_off()
        fig.subplots_adjust(left=0,right=1,bottom=0,top=1,hspace=0)
        fig.savefig(dest / f"peak_vs_DEG_direction_all{suffix}.png",dpi=300,facecolor="white")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--branch", choices=("all", "both", "promoters", "all-genes", "all-peaks"), default="all",
                        help="Default: all three branches; both retains the two gene branches")
    parser.add_argument("--publication-figures", action="store_true", help="Also render text-free PNGs within each branch")
    args = parser.parse_args()
    if args.dry_run:
        print(f"[DRY-RUN] Regulatory target branch: {args.branch}")
        return
    promoter_input = CUTRUN_ROOT / "data" / "figure_inputs" / "promoter_gene_venn"
    if args.branch in ("all", "both", "promoters") and not promoter_input.is_dir():
        raise FileNotFoundError("Run CUT&RUN Step 05 before regulatory integration")
    associated_input = CUTRUN_ROOT / "data/figure_inputs/protein_coding_peak_associations/PeakGeneAssignments.tsv"
    if args.branch in ("all", "both", "all-genes", "all-peaks") and not associated_input.is_file():
        raise FileNotFoundError("Run CUT&RUN Step 06 before all-gene regulatory integration")
    selected = []
    if args.branch in ("all", "both", "promoters"):
        selected.append("promoters")
    if args.branch in ("all", "both", "all-genes"):
        selected.append("all_genes")
    if args.branch in ("all", "all-peaks"):
        selected.append("all_peaks")
    for mode in selected:
        output = BRANCHES[mode]
        output.mkdir(parents=True, exist_ok=True)
        if mode == "all_genes":
            prepare_shared_peak_genes(output)
        if mode == "all_peaks":
            prepare_peak_targets(output, args.publication_figures)
        else:
            run_r("render_direction.R", mode, output, args.publication_figures)
        run_r("render_targets.R", mode, output, args.publication_figures)
        print(f"[DONE] Regulatory-target figures ({mode}): {output}")


if __name__ == "__main__":
    main()
