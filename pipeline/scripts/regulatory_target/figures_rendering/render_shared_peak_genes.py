#!/usr/bin/env python3
"""Describe shared protein-coding-associated peaks and intersect their genes with DEGs."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mcm3-promotion-mpl')
import matplotlib
matplotlib.use('Agg')
import pandas as pd

FACTORS = ('MCM3', 'NONO', 'PSPC1')
COLORS = ('#A80F14', '#F39C12', '#1F4AA8')

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def render_directions(root, out, inputs, publication_figures):
    data,vis=out/'data',out
    data.mkdir(parents=True,exist_ok=True)
    vis.mkdir(parents=True,exist_ok=True)
    if publication_figures:
        (vis/'publication_figures').mkdir(parents=True,exist_ok=True)
    members=pd.read_csv(inputs/'Venn_Peaks_members.tsv',sep='\t',dtype=str)
    ids=set(members.loc[members.region.eq('111'),'member'])
    tables={f:pd.read_csv(inputs/f'Peaks_{f}.tsv',sep='\t') for f in FACTORS}
    shared=tables['MCM3'].loc[tables['MCM3'].peak_id.isin(ids)].copy()
    assert len(shared)==len(ids)
    assignments=pd.read_csv(inputs/'PeakGeneAssignments.tsv',sep='\t')
    links=assignments.loc[assignments.peak_id.isin(ids),['peak_id','gene_id','gene']].drop_duplicates()
    assert set(links.peak_id)==ids
    genes=set(links.gene_id)
    symbols=links.drop_duplicates('gene_id').set_index('gene_id').gene.to_dict()
    def gene_list(ids):
        return pd.DataFrame([{'gene_id':g,'gene':symbols[g]} for g in sorted(ids)],columns=['gene_id','gene'])
    links.to_csv(data/'SharedPeakGeneIDAssignments.tsv',sep='\t',index=False)
    shared.to_csv(data/'SharedPeaks.tsv',sep='\t',index=False)
    shared[['chrom','start','end']].to_csv(data/'SharedPeaks.bed',sep='\t',index=False,header=False)
    links.groupby(['gene_id','gene']).agg(n_shared_peaks=('peak_id','nunique')).reset_index().to_csv(
        data/'SharedPeakAssociatedGenes.tsv',sep='\t',index=False)

    helper=root/'scripts/CUT&RUN/figures_rendering'
    renderer=load(helper/'peak_annotation.py','annotation')
    distribution=load(helper/'render_peak_distribution.py','distribution')
    renderer.OUTDIR=data/'annotation'
    renderer.TMPDIR=renderer.OUTDIR/'intermediate'
    renderer.PAPER_TMPDIR=renderer.OUTDIR/'paperfig_intermediate'
    renderer.GTF_CANDIDATES=[root/'reference/gencode.vM25.annotation.gtf']
    renderer.BEDTOOLS=Path(shutil.which('bedtools') or '/opt/anaconda3/envs/cutrun_env/bin/bedtools')
    renderer.ensure_dirs()
    categories=distribution.category_beds(renderer,data/'SharedPeaks.bed')
    mid=renderer.TMPDIR/'SharedPeaks_mid.bed'
    renderer.bed_to_peak_midpoints(data/'SharedPeaks.bed',mid)
    counts=renderer.assign_categories(mid,categories)
    pd.DataFrame([{'category':k,'n_peaks':v,'percent':100*v/len(shared)} for k,v in counts.items()]).to_csv(
        data/'SharedPeaks_distribution.tsv',sep='\t',index=False)
    promoter=counts['Promoter (\u00b11kb)']

    direction=load(root/'scripts/regulatory_target/figures_rendering/direction_panel.py','direction')
    summary=[]
    duplicate_rows=[]
    matching_changes=[]
    mapping_audit=pd.read_csv(root/'regulatory_work/data/GeneSymbolMappingAudit.tsv',sep='\t')
    for i,f in enumerate(FACTORS):
        sets={}
        for deg_direction in ('UP','DOWN'):
            table=pd.read_csv(root/f'deg_work/data/VennDiagram_{deg_direction}_{f}_genes.tsv',sep='\t')
            table['gene_id']=table.gene_id.str.replace(r'\..*','',regex=True)
            assert table.gene_id.is_unique
            assert (table['adj.P.Val']<=0.05).all() and (table.logFC.abs()>=0.28).all()
            assert (table.logFC.gt(0) if deg_direction=='UP' else table.logFC.lt(0)).all()
            sets[deg_direction]=set(table.gene_id)
            audit=mapping_audit.loc[mapping_audit.factor.eq(f)&mapping_audit.direction.eq(deg_direction.title())]
            assert set(audit.gene_id)==sets[deg_direction]
            duplicates=audit.loc[audit.final_gene_symbol.notna()&audit.final_gene_symbol.duplicated(False)]
            duplicate_rows.extend(duplicates.to_dict('records'))
            table.loc[table.gene_id.isin(genes)].to_csv(data/f'SharedPeakGenes_{f}_{deg_direction}_DEGstatistics.tsv',sep='\t',index=False)
        up,down=sets['UP'],sets['DOWN']
        if up & down:
            raise ValueError(f'Overlapping Up and Down gene IDs for {f}')
        gu,gd=genes&up,genes&down
        for label,id_set in [('Up',up),('Down',down)]:
            audit=mapping_audit.loc[mapping_audit.factor.eq(f)&mapping_audit.direction.eq(label)]
            old_symbols=set(audit.final_gene_symbol.dropna())
            for gene_id in sorted(genes):
                symbol=symbols[gene_id]
                old_match=symbol in old_symbols
                current_match=gene_id in id_set
                if old_match!=current_match:
                    matching_changes.append({'factor':f,'direction':label,'gene_id':gene_id,'assigned_gene_symbol':symbol,
                        'previous_symbol_match':old_match,'current_gene_id_match':current_match,
                        'DEG_IDs_with_same_symbol':';'.join(audit.loc[audit.final_gene_symbol.eq(symbol),'gene_id']),
                        'DEG_symbol_for_same_ID':';'.join(audit.loc[audit.gene_id.eq(gene_id),'final_gene_symbol'].dropna())})
        for label,subset in [('Up',gu),('Down',gd),('NotDEG',genes-up-down)]:
            gene_list(subset).to_csv(data/f'SharedPeakGenes_{f}_{label}.tsv',sep='\t',index=False)
        row={'factor':f,'shared_peak_associated_genes':len(genes),'up_DEGs':len(up),'down_DEGs':len(down),
             'overlap_up':len(gu),'overlap_down':len(gd),'not_DEG':len(genes-up-down)}
        summary.append(row)
        pd.DataFrame({'gene_id':sorted(genes),'gene':[symbols[g] for g in sorted(genes)],'upregulated':[g in up for g in sorted(genes)],
                      'downregulated':[g in down for g in sorted(genes)]}).to_csv(
            data/f'SharedPeakGenes_vs_DEG_{f}.tsv',sep='\t',index=False)
        for hidden in ([False,True] if publication_figures else [False]):
            dest=vis/'publication_figures' if hidden else vis
            suffix='_noTexts' if hidden else ''
            from argparse import Namespace
            a=Namespace(out=str(dest/f'SharedPeakGenes_vs_DEG_direction_{f}{suffix}.png'),panel_letter='ABC'[i],
                factor=f,n_center=len(genes),n_up=len(up),n_down=len(down),n_only_center=len(genes-up-down),
                n_only_up=len(up-genes),n_only_down=len(down-genes),n_overlap_up=len(gu),n_overlap_down=len(gd),
                center_col=COLORS[i],col_up='#D87070',col_down='#6FA37A',col_overlap='#BFBFBF',stroke_col='#222222',
                stroke_lwd=2,r_side=.170,r_center=.205,cx=.50,cy=.53,lx=.305,ly=.53,rx=.695,ry=.53,
                width_px=3200,height_px=2400,dpi=300,fs_panel_letter=20,fs_title=24,fs_subtitle=20,
                fs_num_center=22,fs_num_overlap=17,fs_bottom_lab=20,fs_bottom_n=18,hide_numbers=False,no_text=hidden)
            direction._args=lambda:a
            # Preserve the requested geometry while accurately naming the central gene set.
            from matplotlib.axes import Axes
            original=Axes.text
            def text(self,x,y,s,*args,**kwargs):
                if s==f'{f} CUT&RUN': s='Shared-peak-associated genes'
                return original(self,x,y,s,*args,**kwargs)
            Axes.text=text
            try: direction.main()
            finally: Axes.text=original
    pd.DataFrame(summary).to_csv(data/'SharedPeakGenes_vs_DEG_summary.tsv',sep='\t',index=False)
    pd.DataFrame(duplicate_rows).to_csv(data/'DEG_DuplicateSymbolAudit.tsv',sep='\t',index=False)
    pd.DataFrame(matching_changes).to_csv(data/'IdentifierMatchingChanges.tsv',sep='\t',index=False)
    # Use the same R raster stacking and dimensions as promoter_vs_DEG_direction_all.png.
    rscript=shutil.which('Rscript') or '/opt/anaconda3/envs/mcm3-pipeline/bin/Rscript'
    for hidden in ([False,True] if publication_figures else [False]):
        dest=vis/'publication_figures' if hidden else vis
        suffix='_noTexts' if hidden else ''
        paths=[dest/f'SharedPeakGenes_vs_DEG_direction_{f}{suffix}.png' for f in FACTORS]
        code='''args <- commandArgs(TRUE)
library(png); library(grid); library(gridExtra)
grobs <- lapply(args[1:3], function(p) rasterGrob(readPNG(p), interpolate=TRUE))
g <- arrangeGrob(grobs[[1]],grobs[[2]],grobs[[3]],ncol=1)
png(args[4],width=2200,height=5200,res=300)
grid.newpage();grid.draw(g);dev.off()
'''
        subprocess.run([rscript,'-',*[str(p) for p in paths],str(dest/f'SharedPeakGenes_vs_DEG_direction_all{suffix}.png')],
            input=code,text=True,check=True,cwd=out)
    parameters={'n_shared_peaks':len(shared),'n_gene_assigned_peaks':len(shared),'n_unique_associated_genes':len(genes),
        'promoter_midpoint_peaks':promoter,'outside_promoter_midpoint_peaks':len(shared)-promoter,
        'gene_assignment':'All gene promoters at TSS +/-1kb overlapping by >=250bp; nearest-TSS fallback only without any qualifying promoter; protein-coding filter after assignment',
        'DEG_criteria':'BH-adjusted p <=0.05; abs(log2FC) >=0.28; existing final directional gene lists and symbol mapping',
        'counting_unit':'Distinct Ensembl gene IDs (version removed) throughout both shared-peak and DEG sets; symbols are labels only',
        'interpretation':'Exploratory peak-associated DEG candidates, not replacements for promoter-bound regulatory targets',
        'not_DEG':'Not in either significant DEG list; includes nonsignificant and untested genes'}
    (data/'AnalysisParameters.json').write_text(json.dumps(parameters,indent=2)+'\n')
    tracked=[inputs/'Venn_Peaks_members.tsv',inputs/'PeakGeneAssignments.tsv',*[inputs/f'Peaks_{f}.tsv' for f in FACTORS],
        root/'regulatory_work/data/GeneSymbolMappingAudit.tsv',
        *[root/f'deg_work/data/VennDiagram_{d}_{f}_genes.tsv' for f in FACTORS for d in ('UP','DOWN')]]
    (data/'DirectionInputChecksums.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tracked},indent=2)+'\n')
    (out/'README.md').write_text('# Shared-peak-associated genes: exploratory analysis\n\n'
        f'All {len(shared):,} shared peaks have existing protein-coding gene assignments by construction. '
        f'They represent {len(genes):,} distinct Ensembl gene IDs. Assignment is promoter-priority with nearest-TSS fallback, not proof of regulation.\n\n'
        'All sets and overlaps use Ensembl gene IDs with version suffixes removed. Symbols are labels only. '
        'DEG_DuplicateSymbolAudit.tsv explains differences from symbol-deduplicated counts. '
        'Counts are read from deg_work/data/VennDiagram_<UP|DOWN>_<FACTOR>_genes.tsv, without rerunning DE analysis.\n\n'
        'The distribution pie is generated only under cutrun_work/visuals, using canonical-transcript midpoint categories. '
        'Promoter localization means midpoint within TSS +/-1 kb, not the promoter-bound gene filter. '
        'A peak outside an annotated gene body can still be assigned to a nearby gene.\n\n'
        'The directional diagrams use the same central associated-gene set in every panel and the '
        'existing target-specific Up/Down DEG sets (BH-adjusted p <=0.05, abs(log2FC) >=0.28). '
        'The central-only count means not a significant DEG, including untested genes. '
        'These are descriptive intersections, not new significance or enrichment tests. '
        'No production regulatory-target membership is changed.\n\n'
        'The layout, colors and dimensions follow promoter_vs_DEG_direction_all.png; '
        'the title identifies shared-peak-associated genes instead of promoter-bound genes. '
        'Data include peak assignments, unique genes, directional overlaps and input checksums.\n')
    print(json.dumps(parameters,indent=2))
    print(pd.DataFrame(summary).to_string(index=False))
    print('[DONE] Shared-peak DEG direction figures:',out)


def render_targets(root, out, publication_figures):
    source=out
    data,vis=out/'data',out
    data.mkdir(parents=True,exist_ok=True)
    vis.mkdir(parents=True,exist_ok=True)
    if publication_figures:
        (vis/'publication_figures').mkdir(parents=True,exist_ok=True)
    gene_table=pd.read_csv(source/'data/SharedPeakAssociatedGenes.tsv',sep='\t')
    assert gene_table.gene_id.is_unique
    genes=set(gene_table.gene_id)
    names={}
    gtf=root/'reference/gencode.vM25.annotation.gtf'
    with gtf.open() as handle:
        for line in handle:
            fields=line.rstrip().split('\t')
            if len(fields)!=9 or fields[2]!='gene':continue
            a=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
            names[a['gene_id'].split('.')[0]]=a.get('gene_name',a['gene_id'])
    def write_list(path,ids):
        pd.DataFrame([{'gene_id':g,'gene':names.get(g,g)} for g in sorted(ids)],columns=['gene_id','gene']).to_csv(path,sep='\t',index=False)
    sets,summary={},[]
    known_ids=set(genes)
    inputs=[source/'data/SharedPeakAssociatedGenes.tsv',source/'data/SharedPeakGenes_vs_DEG_summary.tsv',gtf]
    for factor in FACTORS:
        directions={}
        for label in ['Up','Down']:
            path=root/f'deg_work/data/VennDiagram_{label.upper()}_{factor}_genes.tsv'
            inputs.append(path)
            table=pd.read_csv(path,sep='\t')
            table['gene_id']=table.gene_id.str.replace(r'\..*','',regex=True)
            assert table.gene_id.is_unique
            assert (table['adj.P.Val']<=.05).all() and (table.logFC.abs()>=.28).all()
            directions[label]=set(table.gene_id)
            write_list(data/f'SharedPeakGenes_{factor}_{label}.tsv',genes&directions[label])
            table.loc[table.gene_id.isin(genes)].to_csv(data/f'{factor}_{label}_DEGstatistics.tsv',sep='\t',index=False)
            known_ids.update(directions[label])
        assert not directions['Up']&directions['Down']
        all_degs=directions['Up']|directions['Down']
        sets[factor]=genes&all_degs
        write_list(data/f'Venn_target_{factor}_genes.tsv',sets[factor])
        write_list(data/f'SharedPeakGenes_{factor}_OtherDEGs.tsv',all_degs-genes)
        summary.append({'factor':factor,'up_DEGs':len(directions['Up']),'down_DEGs':len(directions['Down']),
                        'overlap_up':len(genes&directions['Up']),'overlap_down':len(genes&directions['Down']),
                        'associated_DEGs':len(sets[factor])})
    summary=pd.DataFrame(summary)
    previous=pd.read_csv(source/'data/SharedPeakGenes_vs_DEG_summary.tsv',sep='\t')
    columns=['factor','up_DEGs','down_DEGs','overlap_up','overlap_down']
    pd.testing.assert_frame_equal(summary[columns],previous[columns])
    summary.to_csv(data/'SharedPeakDEGs_summary.tsv',sep='\t',index=False)
    write_list(data/'GeneIdentifiers.tsv',known_ids)
    a,b,c=(sets[f] for f in FACTORS)
    regions={'100':a-b-c,'010':b-a-c,'001':c-a-b,'110':(a&b)-c,'101':(a&c)-b,'011':(b&c)-a,'111':a&b&c}
    labels={'100':'MCM3_only','010':'NONO_only','001':'PSPC1_only','110':'MCM3_NONO_only',
            '101':'MCM3_PSPC1_only','011':'NONO_PSPC1_only','111':'MCM3_NONO_PSPC1'}
    pd.DataFrame([{'region':labels[k],'gene_id':g,'gene':names.get(g,g)} for k,values in regions.items() for g in sorted(values)],
                 columns=['region','gene_id','gene']).to_csv(data/'Venn_target_genes.tsv',sep='\t',index=False)
    pd.DataFrame([{'region':k,'n':len(v)} for k,v in regions.items()]).to_csv(data/'Venn_target_counts.tsv',sep='\t',index=False)
    write_list(data/'Venn_target_shared_genes.tsv',regions['111'])
    (data/'AnalysisParameters.tsv').write_text(
        'analysis_role\texploratory shared-peak-associated DEGs, not promoter-bound direct targets\n'
        'gene_set\tprotein-coding genes assigned to canonical peak loci shared by MCM3, NONO and PSPC1\n'
        'gene_assignment\tall promoters overlapping >=250bp at TSS +/-1kb; nearest-TSS fallback only if no qualifying promoter; protein-coding filter afterward\n'
        'counting_unit\tEnsembl gene ID with version removed; symbols are labels\n'
        'DEG_threshold\tBH-adjusted p <=0.05; abs(log2FC) >=0.28\n'
        'groups\tA: all three KDs; B: MCM3+PSPC1 KDs only; C: NONO+PSPC1 KDs only; D: MCM3+NONO KDs only\n'
        'group_meaning\tDEG overlap, not factor-exclusive binding; all input genes have three-factor shared-peak association\n'
        'profiles\tGENCODE M25 gene bodies; 3kb upstream, 3kb scaled body, 3kb downstream; 25bp bins; existing mean bigWigs\n')
    env=os.environ.copy()
    env['MPLCONFIGDIR']=str(out/'.matplotlib')
    env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
    env['PYTHONDONTWRITEBYTECODE']='1'
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[name]='1'
    helper=root/'scripts/regulatory_target/figures_rendering/venn.py'
    for hidden in ([False,True] if publication_figures else [False]):
        dest=vis/'publication_figures/Venn_target_noTexts.png' if hidden else vis/'Venn_target.png'
        command=[sys.executable,str(helper),'--out',str(dest),'--title','Shared-peak-associated DEGs']
        for key,f in zip(['a','b','c'],FACTORS):command.extend([f'--{key}-name',f,f'--{key}-total',str(len(sets[f]))])
        for key,values in regions.items():command.extend([f'--n{key}',str(len(values))])
        if hidden:command.append('--hide-numbers')
        subprocess.run(command,check=True,env=env)
    rscript=shutil.which('Rscript',path=env['PATH'])
    if not rscript:
        raise FileNotFoundError('Rscript')
    script=Path(__file__).with_name('render_shared_peak_targets.R')
    command=[rscript,str(script),str(root),str(out)]
    if publication_figures:
        command.append(str(vis/'publication_figures'))
    subprocess.run(command,check=True,env=env,cwd=out)
    inputs.extend([helper,script,Path(__file__)])
    (data/'InputChecksums.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2)+'\n')
    (out/'README.md').write_text(
        '# Shared-peak-associated DEG figures\n\n'
        'All inputs are genes associated with three-factor shared CUT&RUN peaks. The target sets '
        'intersect those genes with each KD DEG list by Ensembl ID. These are candidate associations, '
        'not established promoter-bound direct regulatory targets. The separate promoter-bound analysis is unchanged.\n\n'
        'Group A is DEG in all three KDs; B in MCM3 and PSPC1 only; C in NONO and PSPC1 only; '
        'D in MCM3 and NONO only. Only refers to DEG status, not CUT&RUN binding. '
        'Up/Down/Mix compares the RNA-seq directions of the indicated KDs.\n\n'
        'Pie_TargetBinding filenames are retained for comparison with the production design, '
        'but labels and tables describe shared-peak gene association among DEGs. '
        'The metaprofiles use one GENCODE gene body per Ensembl ID, the original bigWigs, '
        'and the same 25bp, +/-3kb / 3kb-scaled-body settings and colors as production. '
        'They are gene-body profiles, not averages centered on the shared peak loci.\n\n'
        'Figures are in this folder; noTexts counterparts are in publication_figures. '
        'Source sets, group directions, counts, BEDs, matrices, and parameters are in data.\n')
    print(summary.to_string(index=False),flush=True)
    print('Venn regions:',{k:len(v) for k,v in regions.items()},flush=True)
    print('[DONE] Shared-peak-associated DEG figures:',out,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pipeline-root',type=Path,required=True)
    parser.add_argument('--publication-figures',action='store_true')
    args=parser.parse_args()
    root=args.pipeline_root.resolve()
    inputs=root/'cutrun_work/data/figure_inputs/protein_coding_peak_associations'
    if not (inputs/'PeakGeneAssignments.tsv').is_file():
        raise FileNotFoundError('Run CUT&RUN Step 06 to generate promoter-priority assignments first')
    out=root/'regulatory_work/visuals/shared_peaks_visuals'
    render_directions(root,out,inputs,args.publication_figures)
    render_targets(root,out,args.publication_figures)


if __name__=='__main__':
    main()
