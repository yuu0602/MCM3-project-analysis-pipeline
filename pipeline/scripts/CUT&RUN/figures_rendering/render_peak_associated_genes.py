#!/usr/bin/env python3
"""Render peak Venns, independent peak-associated gene sets, and distribution pies."""

from bisect import bisect_left
from collections import defaultdict
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import pandas as pd

PROJECT = Path(__file__).resolve().parent
RUN = PROJECT.parents[2]
CUTRUN = RUN / "cutrun_work"
DATA = CUTRUN / "data" / "figure_inputs"
GENE_VISUALS = CUTRUN / "visuals" / "all_protein_coding_genes_visuals"
PEAK_VISUALS = CUTRUN / "visuals" / "all_protein_coding_peaks_visuals"
NO_TEXT_VISUALS: Path | None = None
TEXT_FREE = False
RENDER_GENES = True
RENDER_PEAKS = True
GENE_MEMBERSHIP = "independent_genes"
GTF = RUN / "reference" / "gencode.vM25.annotation.gtf"
FACTORS = ("MCM3", "NONO", "PSPC1")
ASSAYS = (*FACTORS, "IgG")
COORDS = ["chrom", "start", "end"]
LOCUS_REGIONS = {
    1: "100_only_MCM3", 2: "010_only_NONO", 3: "110_MCM3_NONO",
    4: "001_only_PSPC1", 5: "101_MCM3_PSPC1", 6: "011_NONO_PSPC1",
    7: "111_MCM3_NONO_PSPC1",
}


def locus_defined_gene_groups(output: Path) -> pd.DataFrame:
    """Keep the peak's overlap class when assigning genes; never union factor bits across loci."""
    loci = pd.read_csv(CUTRUN / "data/PeakLoci/Venn_Peaks_loci.tsv", sep="\t")
    assignments = pd.read_csv(output / "PeakGeneAssignments.tsv", sep="\t")
    links = assignments.merge(
        loci[["locus_id", "region_mask"]], left_on="peak_id", right_on="locus_id",
        how="inner", validate="many_to_one",
    )
    if not links.gene_type.eq("protein_coding").all():
        raise ValueError("Locus-defined groups require protein-coding assignments")
    links["region"] = links.region_mask.map(LOCUS_REGIONS)
    groups = links[["region", "gene_id", "gene"]].drop_duplicates().sort_values(["region", "gene_id"])
    links.to_csv(output / "LocusDefinedPeakGeneAssignments.tsv", sep="\t", index=False)
    groups.to_csv(output / "LocusDefinedGeneGroups.tsv", sep="\t", index=False)
    counts = links.groupby(["region_mask", "region"]).agg(
        n_genes=("gene_id", "nunique"), n_associated_peaks=("peak_id", "nunique")
    ).reset_index()
    counts.to_csv(output / "LocusDefinedGeneGroups_counts.tsv", sep="\t", index=False)
    shared = links.loc[links.region_mask.eq(7)].groupby("gene_id", as_index=False).agg(
        gene=("gene", "first"), n_shared_peaks=("peak_id", "nunique")
    )
    shared.to_csv(output / "SharedPeakAssociatedGenes.tsv", sep="\t", index=False)
    (output / "LocusDefinedGeneGroups_parameters.json").write_text(json.dumps({
        "source_loci": "cutrun_work/data/PeakLoci/Venn_Peaks_loci.tsv",
        "membership": "Exact overlapping-factor class of each canonical genomic locus, assigned to genes afterward",
        "assignment": "Existing Step 06 all qualifying promoter links (TSS +/-1000 bp; overlap >=250 bp), otherwise nearest-TSS fallback; protein-coding filter afterward",
        "counting_unit": "Unique Ensembl gene IDs within each locus class",
        "groups_are_disjoint": False,
        "role": "Membership source for shared_locus figures; not for independent_genes figures",
        "interpretation": "A gene may occur in several locus classes; only describes factor membership at a locus, not exclusive binding across the gene",
        "profiles": "shared_locus profiles use these groups; independent_genes profiles use independent per-factor gene-set intersections",
    }, indent=2) + "\n")
    return groups


def gene_tss() -> dict[str, list[tuple[int, str]]]:
    rows = []
    with GTF.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = fields[8]
            if 'gene_name "' not in attributes:
                continue
            symbol = attributes.split('gene_name "')[1].split('"')[0]
            start, end = int(fields[3]) - 1, int(fields[4])
            rows.append((fields[0], start if fields[6] == "+" else end, symbol))
    table = pd.DataFrame(rows, columns=["chrom", "tss", "gene"]).drop_duplicates()
    return {chrom: list(group[["tss", "gene"]].itertuples(index=False, name=None)) for chrom, group in table.sort_values(["chrom", "tss"]).groupby("chrom", sort=False)}


def nearest_gene(chrom: str, start: int, end: int, index: dict[str, list[tuple[int, str]]]) -> tuple[str, int] | None:
    entries = index.get(chrom)
    if not entries:
        return None
    center = (start + end) // 2
    positions = [entry[0] for entry in entries]
    point = bisect_left(positions, center)
    candidates = entries[max(0, point - 1):point + 1]
    tss, gene = min(candidates, key=lambda entry: (abs(entry[0] - center), entry[1]))
    return gene, abs(tss - center)


def peak_sets() -> dict[str, pd.DataFrame]:
    loci = pd.read_csv(CUTRUN / "data" / "PeakLoci" / "Venn_Peaks_loci.tsv", sep="\t")
    sets = {}
    for factor, bit in zip(FACTORS, (1, 2, 4)):
        sets[factor] = loci.loc[loci.region_mask.astype(int).map(lambda mask: bool(mask & bit)), ["chrom", "start", "end"]].copy()
    igg = pd.read_csv(CUTRUN / "data" / "PeakSets" / "Peaks_IgG.bed", sep="\t", header=None, names=["chrom", "start", "end"])
    sets["IgG"] = igg
    return sets


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gene_annotation_renderer(output: Path):
    distribution = load_module(PROJECT / "render_peak_distribution.py", "peak_distribution")
    renderer = load_module(PROJECT / "peak_annotation.py", "peak_annotation")
    renderer.WORK_ROOT = CUTRUN
    renderer.OUTDIR = output / "annotation"
    renderer.TMPDIR = renderer.OUTDIR / "intermediate"
    renderer.PAPER_TMPDIR = renderer.OUTDIR / "paperfig_intermediate"
    renderer.GTF_CANDIDATES = [GTF]
    renderer.ASSAYS = list(FACTORS)
    renderer.PEAK_UNIT_LABEL = "peak-associated genes"
    renderer.USE_NT_KD_UNION = False
    renderer.BEDTOOLS = Path(shutil.which("bedtools") or "/opt/anaconda3/envs/cutrun_env/bin/bedtools")
    renderer.ensure_dirs()
    return renderer, distribution.category_beds(renderer, CUTRUN / "data" / "PeakLoci" / "Venn_Peaks_MCM3.bed")


def annotate_unique_genes(renderer, categories: dict[str, Path], table: pd.DataFrame, factor: str) -> dict[str, int]:
    bed = renderer.TMPDIR / f"{factor}_unique_peak_associated_genes.bed"
    table[["chrom", "start", "end"]].to_csv(bed, sep="\t", header=False, index=False)
    midpoint = renderer.TMPDIR / f"{factor}_unique_peak_associated_genes_mid.bed"
    renderer.bed_to_peak_midpoints(bed, midpoint)
    return renderer.assign_categories(midpoint, categories)


def protein_coding_assignments(output: Path):
    """Assign all qualifying promoters, fall back to nearest TSS, then filter biotype."""
    data = output
    sources=peak_sets()
    peak_members={}
    for factor,table in sources.items():
        peak_members[factor]={f'{r.chrom}:{int(r.start)}-{int(r.end)}' for r in table.itertuples(index=False)}
    peaks=pd.concat(sources.values(),ignore_index=True).drop_duplicates(COORDS).sort_values(COORDS).reset_index(drop=True)
    peaks['peak_id']=[f'{r.chrom}:{int(r.start)}-{int(r.end)}' for r in peaks.itertuples(index=False)]
    peaks[COORDS+['peak_id']].to_csv(data/'OriginalPeaks.bed',sep='\t',index=False,header=False)
    rows=[]
    with GTF.open() as h:
        for line in h:
            fields=line.rstrip().split('\t')
            if len(fields)!=9 or fields[2]!='gene':continue
            a=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
            tss=int(fields[3])-1 if fields[6]=='+' else int(fields[4])
            rows.append({'chrom':fields[0],'start':max(0,tss-1000),'end':tss+1000,'tss':tss,
                         'gene_id':a['gene_id'].split('.')[0],'gene':a.get('gene_name',a['gene_id']),
                         'gene_type':a['gene_type']})
    promoters=pd.DataFrame(rows)
    promoters[['chrom','start','end','gene_id','gene','gene_type','tss']].to_csv(
        data/'Promoters_AllBiotypes.bed',sep='\t',index=False,header=False)
    bedtools=shutil.which('bedtools') or '/opt/anaconda3/envs/cutrun_env/bin/bedtools'
    with (data/'PeakPromoter_overlaps.tsv').open('w') as h:
        subprocess.run([bedtools,'intersect','-a',str(data/'OriginalPeaks.bed'),'-b',str(data/'Promoters_AllBiotypes.bed'),'-wo'],stdout=h,check=True)
    overlaps=pd.read_csv(data/'PeakPromoter_overlaps.tsv',sep='\t',header=None,
        names=['chrom','start','end','peak_id','promoter_chrom','promoter_start','promoter_end','gene_id','gene','gene_type','tss','promoter_overlap_bp'])
    passing=overlaps.loc[overlaps.promoter_overlap_bp.ge(250)].copy()
    qualifying=set(passing.peak_id)
    passing['assignment_method']='promoter_overlap'
    passing['nearest_tss_distance_bp']=abs((passing.start+passing.end)//2-passing.tss)
    cols=COORDS+['peak_id','gene_id','gene','gene_type','assignment_method','promoter_overlap_bp','nearest_tss_distance_bp']
    passing=passing[cols].drop_duplicates(['peak_id','gene_id'])
    index=gene_tss()
    positions={chrom:[r[0] for r in entries] for chrom,entries in index.items()}
    annotation=defaultdict(list)
    for r in promoters.itertuples(index=False):annotation[(r.chrom,r.gene)].append(r)
    fallback=[]
    for r in peaks.loc[~peaks.peak_id.isin(qualifying)].itertuples(index=False):
        center=(int(r.start)+int(r.end))//2
        entries=index.get(r.chrom,[])
        if not entries:
            fallback.append((r.chrom,r.start,r.end,r.peak_id,'','','unassigned','unassigned',0,-1))
            continue
        point=bisect_left(positions[r.chrom],center)
        tss,gene=min(entries[max(0,point-1):point+1],key=lambda e:(abs(e[0]-center),e[1]))
        matches=[a for a in annotation[(r.chrom,gene)] if abs(a.tss-center)==abs(tss-center)]
        if not matches:raise ValueError(f'Unresolved nearest-TSS locus: {r.peak_id}')
        # Retain every tied original annotation record; do not redirect to coding genes.
        for a in matches:
            fallback.append((r.chrom,r.start,r.end,r.peak_id,a.gene_id,gene,a.gene_type,'nearest_TSS_fallback',0,abs(tss-center)))
    audit=pd.concat([passing,pd.DataFrame(fallback,columns=cols)],ignore_index=True).drop_duplicates(['peak_id','gene_id'])
    audit['retained']=audit.gene_type.eq('protein_coding')
    audit.to_csv(data/'PeakGeneAssignments_all_biotypes_audit.tsv',sep='\t',index=False)
    assignments=audit.loc[audit.retained].copy()
    assignments.to_csv(data/'PeakGeneAssignments.tsv',sep='\t',index=False)
    assert set(audit.peak_id)==set(peaks.peak_id)
    assert not set(audit.loc[audit.assignment_method.eq('nearest_TSS_fallback'),'peak_id'])&qualifying
    peak_tables,gene_tables,summaries={},{},[]
    for factor in ASSAYS:
        links=assignments.loc[assignments.peak_id.isin(peak_members[factor])].copy()
        kept=peaks.loc[peaks.peak_id.isin(set(links.peak_id))].copy()
        genes=links.sort_values(['gene_id','nearest_tss_distance_bp','chrom','start','end']).drop_duplicates('gene_id')
        assert kept.peak_id.is_unique and genes.gene_id.is_unique
        kept.to_csv(data/f'Peaks_{factor}.tsv',sep='\t',index=False)
        kept[COORDS].to_csv(data/f'Peaks_{factor}.bed',sep='\t',index=False,header=False)
        genes.to_csv(data/f'PeakAssociatedGenes_{factor}.tsv',sep='\t',index=False)
        links.to_csv(data/f'PeakGeneAssignments_{factor}.tsv',sep='\t',index=False)
        peak_tables[factor],gene_tables[factor]=kept,genes
        summaries.append({'factor':factor,'original_peaks':len(sources[factor]),'retained_peaks':len(kept),
                          'associated_gene_ids':len(genes),'promoter_associated_gene_ids':links.loc[links.assignment_method.eq('promoter_overlap'),'gene_id'].nunique()})
    pd.DataFrame(summaries).to_csv(data/'AssignmentSummary.tsv',sep='\t',index=False)
    return gene_tables, peak_tables


def render_venn(venn, tables, column: str, name: str, output: Path, visual_dir: Path):
    a, b, c = (set(tables[f][column]) for f in FACTORS)
    regions = {"100": a-b-c, "010": b-a-c, "001": c-a-b,
               "110": (a & b)-c, "101": (a & c)-b, "011": (b & c)-a,
               "111": a & b & c}
    pd.DataFrame([{"region": key, "count": len(value)} for key, value in regions.items()]).to_csv(
        output / f"{name}_counts.tsv", sep="\t", index=False)
    pd.DataFrame([{"region": key, "member": item} for key, value in regions.items()
                  for item in sorted(value)], columns=["region", "member"]).to_csv(
        output / f"{name}_members.tsv", sep="\t", index=False)
    venn.plot_triple_venn((a, b, c), FACTORS, visual_dir / f"{name}.png",
                          show_numbers=not TEXT_FREE, show_totals=not TEXT_FREE)


def shared_peak_distribution(renderer, categories, peaks, output: Path) -> None:
    shared_ids = set.intersection(*(set(peaks[f].peak_id) for f in FACTORS))
    shared = peaks["MCM3"].loc[peaks["MCM3"].peak_id.isin(shared_ids)].copy()
    shared.to_csv(output / "SharedPeaks.tsv", sep="\t", index=False)
    counts = annotate_unique_genes(renderer, categories, shared, "SharedPeaks")
    if sum(counts.values()) != len(shared):
        raise ValueError("Shared peak distribution totals do not match the Venn")
    renderer.PEAK_UNIT_LABEL = "shared peaks"
    renderer.save_single(PEAK_VISUALS / "Pie_SharedPeaks_Distribution.png", "MCM3-NONO-PSPC1", counts, len(shared))
    pd.DataFrame([{"category": key, "n_peaks": value, "percent": 100*value/len(shared)}
                  for key, value in counts.items()]).to_csv(
        output / "Pie_SharedPeaks_Distribution_counts.tsv", sep="\t", index=False)


def main() -> None:
    spec = importlib.util.spec_from_file_location("peak_venn", PROJECT / "peak_venn.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load peak_venn.py")
    venn = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(venn)

    output = DATA / "protein_coding_peak_associations"
    output.mkdir(parents=True, exist_ok=True)
    if RENDER_GENES:
        GENE_VISUALS.mkdir(parents=True, exist_ok=True)
    if RENDER_PEAKS:
        PEAK_VISUALS.mkdir(parents=True, exist_ok=True)
    genes, peaks = protein_coding_assignments(output)
    locus_defined_gene_groups(output)
    retained_ids = set().union(*(set(peaks[f].peak_id) for f in FACTORS))
    loci = pd.read_csv(CUTRUN / "data/PeakLoci/Venn_Peaks_loci.tsv", sep="\t")
    loci.loc[loci.locus_id.isin(retained_ids)].to_csv(
        output / "Venn_Peaks_loci.tsv", sep="\t", index=False
    )
    renderer, categories = gene_annotation_renderer(output)
    for is_gene, tables in ((False, peaks), (True, genes)):
        if (is_gene and not RENDER_GENES) or (not is_gene and not RENDER_PEAKS):
            continue
        visual_dir = GENE_VISUALS if is_gene else PEAK_VISUALS
        prefix = "Pie_PeakAssociatedGenes" if is_gene else "Pie"
        combined = "Pie_PeakAssociatedGenes" if is_gene else "Pie_PeakDistribution"
        renderer.PEAK_UNIT_LABEL = "peak-associated genes" if is_gene else "peaks"
        counts, totals, rows = {}, {}, []
        for factor, table in tables.items():
            counts[factor] = annotate_unique_genes(renderer, categories, table, f"{prefix}_{factor}")
            totals[factor] = len(table)
            if sum(counts[factor].values()) != len(table):
                raise ValueError(f"Pie totals mismatch: {prefix}/{factor}")
            renderer.save_single(visual_dir / f"{prefix}_{factor}.png", factor, counts[factor], len(table))
            rows.extend({"factor": factor, "category": key, "count": value, "total": len(table)}
                        for key, value in counts[factor].items())
        renderer.save_three_panel(visual_dir / f"{combined}.png", counts, totals)
        pd.DataFrame(rows).to_csv(output / f"{combined}_counts.tsv", sep="\t", index=False)
    if RENDER_PEAKS:
        render_venn(venn, peaks, "peak_id", "Venn_Peaks", output, PEAK_VISUALS)
        shared_peak_distribution(renderer, categories, peaks, output)
    if RENDER_GENES and GENE_MEMBERSHIP == "shared_locus":
        gene_output = output / GENE_MEMBERSHIP
        gene_output.mkdir(parents=True, exist_ok=True)
        counts = pd.read_csv(output / "LocusDefinedGeneGroups_counts.tsv", sep="\t")
        counts.assign(region=counts.region.str.split("_").str[0],
                      counting_basis="distinct genes per locus class; non-disjoint categories").to_csv(
            gene_output / "Venn_PeakAssociatedGenes_counts.tsv", sep="\t", index=False)
        groups = pd.read_csv(output / "LocusDefinedGeneGroups.tsv", sep="\t")
        groups.assign(region=groups.region.str.split("_").str[0]).rename(columns={"gene_id": "member"}).to_csv(
            gene_output / "Venn_PeakAssociatedGenes_members.tsv", sep="\t", index=False)
        indexed = counts.set_index("region_mask")
        values = [int(indexed.loc[mask, "n_genes"]) if mask in indexed.index else 0 for mask in range(1, 8)]
        totals = tuple(genes[f].gene_id.nunique() for f in FACTORS)
        venn.plot_region_counts(values, FACTORS, totals, GENE_VISUALS / "Venn_PeakAssociatedGenes.png",
                                show_numbers=not TEXT_FREE, show_totals=not TEXT_FREE, locus_gene_groups=True)
        (gene_output / "Venn_PeakAssociatedGenes_parameters.json").write_text(json.dumps({
            "same_locus_required": True, "groups_are_disjoint": False,
            "membership": "Genes assigned after exact peak-locus intersection",
            "figure": "Schematic locus-defined gene categories, not a conventional gene Venn; counts are not additive",
        }, indent=2) + "\n")
    elif RENDER_GENES:
        gene_output = output / GENE_MEMBERSHIP
        gene_output.mkdir(parents=True, exist_ok=True)
        render_venn(venn, genes, "gene_id", "Venn_PeakAssociatedGenes", gene_output, GENE_VISUALS)
        (gene_output / "Venn_PeakAssociatedGenes_parameters.json").write_text(json.dumps({
            "membership": "Independent per-factor protein-coding peak-associated gene sets",
            "shared_genes": "Intersection of all three gene lists; peaks may occupy different genomic locations",
            "counting_unit": "Unique Ensembl gene ID; seven disjoint gene Venn regions",
            "same_locus_required": False,
            "peak_calling_and_assignment": "Unchanged; only gene-set intersection definition differs from locus-first groups",
        }, indent=2) + "\n")
    print("[DONE] Peak-associated gene figures:", GENE_VISUALS)
    print("[DONE] Whole-genome peak figures:", PEAK_VISUALS)


if __name__ == "__main__":
    main()
