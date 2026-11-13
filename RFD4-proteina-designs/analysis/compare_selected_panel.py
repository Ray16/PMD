#!/usr/bin/env python3
"""Audit the user-confirmed 20-mutant panel and compare global ipTM and geometry.

python compare_selected_panel.py --project-root /homes/rzhu/PMD
Reads existing full-complex Boltz-2 predictions only; performs no new cofolding.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from adjustText import adjust_text
from Bio.Data.PDBData import protein_letters_3to1
from Bio import SeqIO

from plot_attack_geometry import measure


EXPECTED = '''atherm_S190T atherm_T191N atherm_G189I atherm_I15G atherm_I103F
haur_I112D haur_I19A haur_I19V haur_I112L atherm_S138T atherm_I103D haur_F119P
haur_I19C atherm_G102V atherm_T191Q atherm_I15S atherm_I103M atherm_P98C
haur_I112M haur_A17G'''.split()


def read(path, delimiter='\t'):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_csv(path, rows):
    with path.open('w') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def audit_pose(path, expected_sequence):
    data=path.read_bytes(); lig={}; residues={}
    for line in data.decode().splitlines():
        if line.startswith('ATOM  ') and line[21]=='A':
            residues[(int(line[22:26]),line[26])]=line[17:20].strip()
        if not line.startswith('HETATM'):
            continue
        element=(line[76:78].strip() or line[12:16].strip()[0]).upper()
        if element=='H':
            continue
        key=(line[21],line[12:16].strip())
        assert key not in lig, f'Ambiguous ligand atom {path} {key}'
        lig[key]=np.array([float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]])
    sequence=''.join(protein_letters_3to1[residues[k]] for k in sorted(residues))
    assert sequence==expected_sequence, f'Sequence mismatch: {path}'
    expected_mvap={'C18','C19','C20','C21','C22','C23','O11','O12','O13','O14','O15','O16','O17','P24'}
    assert {n for c,n in lig if c=='L'}==expected_mvap
    assert len([n for c,n in lig if c=='T'])==31
    assert all(len([n for c,n in lig if c==chain])==1 for chain in ['M','N'])
    # Validate the hydroxyl on the tertiary MVAP carbon and the beta-gamma bridge.
    for a,b in [(('L','O16'),('L','C23')), (('T','PG'),('T','O3B')), (('T','PB'),('T','O3B'))]:
        assert 1.1 < np.linalg.norm(lig[a]-lig[b]) < 1.95, f'Unexpected ligand bond: {path}'
    assert {n for (c,n),p in lig.items() if c=='L' and n.startswith('C') and
            n!='C23' and np.linalg.norm(p-lig['L','C23'])<1.9}=={'C20','C21','C22'}
    assert {n for (c,n),p in lig.items() if c=='T' and n.startswith('P') and
            np.linalg.norm(p-lig['T','O3B'])<1.95}=={'PB','PG'}
    u=lig['L','O16']-lig['T','PG']; v=lig['T','O3B']-lig['T','PG']
    d=float(np.linalg.norm(u))
    angle=math.degrees(math.acos(float(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))
    return d,angle,sha(data)


def savefig(fig,out,name):
    fig.savefig(out/f'{name}.png',dpi=180,facecolor='white')
    path=out/f'{name}.svg';fig.savefig(path,facecolor='white')
    path.write_text('\n'.join(l.rstrip() for l in path.read_text().splitlines())+'\n')
    plt.close(fig)


def plots(out, summary, poses, rfd):
    fig,axes=plt.subplots(2,2,figsize=(18,12.5),gridspec_kw={'height_ratios':[1,1.1]})
    fig.subplots_adjust(top=.86,bottom=.15,hspace=.28,wspace=.38,left=.06,right=.93)
    for col,arm in enumerate(['atherm','haur']):
        rs=[r for r in summary if r['arm']==arm]
        ax,tab=axes[:,col]
        ax.add_patch(Rectangle((165,0),15,4,color='#edf7ef',zorder=0))
        ax.axvline(165,color='#789879',ls='--',lw=.8);ax.axhline(4,color='#789879',ls='--',lw=.8)
        labels=[]
        for r in rs:
            x,y=r['angle_median_deg'],r['distance_median_A']
            color='#697786' if r['mutation']=='WT' else '#0072b2'
            ax.errorbar(x,y,xerr=np.array([[x-r['angle_q25_deg']],[r['angle_q75_deg']-x]]),
                        yerr=np.array([[y-r['distance_q25_A']],[r['distance_q75_A']-y]]),
                        fmt='D' if r['mutation']=='WT' else 'o',ms=5,color=color,
                        ecolor=color,elinewidth=1,alpha=.65,capsize=2,zorder=2)
            labels.append(ax.text(x,y,r['mutation'],fontsize=8.5,color=color))
        ax.scatter([r['angle_deg'] for r in rfd],[r['distance_A'] for r in rfd],
                   s=100,marker='*',color='#d55e00',edgecolor='white',linewidth=.5,zorder=4)
        ymin=min([r['distance_q25_A'] for r in rs]+[r['distance_A'] for r in rfd])-.12
        ymax=max(r['distance_q75_A'] for r in rs)+.25
        xmin=min(r['angle_q25_deg'] for r in rs)-4
        ax.set(xlim=(xmin,181),ylim=(ymin,ymax),xlabel='Attack angle (°)',ylabel='Attack distance (Å)')
        ax.set_title(('At' if arm=='atherm' else 'Ha')+f': {len(rs)-1} selected mutants + WT',fontsize=13)
        ax.spines[['top','right']].set_visible(False)
        # Keep group labels in a dedicated gutter. Their points stay at the exact
        # measured coordinates; labels are ordered by height to reduce crossings.
        adjust_text(labels,ax=ax,expand=(1.5,1.6),force_text=(.35,.6),
                    force_static=(.3,.5),time_lim=4,
                    arrowprops=dict(arrowstyle='-',color='#808080',lw=.55))
        for index,r in enumerate(sorted(rfd,key=lambda r:(r['distance_A'],r['angle_deg']))):
            ax.annotate(r['group'],xy=(r['angle_deg'],r['distance_A']),
                        xycoords='data',xytext=(1.035,.055+.068*index),
                        textcoords='axes fraction',ha='left',va='center',fontsize=10,
                        weight='bold',color='#a34000',annotation_clip=False,
                        arrowprops=dict(arrowstyle='-',color='#ba7953',lw=.65,
                                        shrinkA=3,shrinkB=4,connectionstyle='arc3,rad=0'))
        ax.text(1.035,.625,'RFD4',transform=ax.transAxes,fontsize=10,
                weight='bold',color='#a34000',ha='left',va='bottom')
        tab.axis('off')
        cells=[[r['mutation'],f"{r['global_iptm_mean']:.5f}",f"{r['global_iptm_delta_vs_WT']:+.5f}",
                f"{r['distance_median_A']:.3f}",f"{r['angle_median_deg']:.2f}",f"{r['NAC_count']}/50"] for r in rs]
        table=tab.table(cellText=cells,colLabels=['Mutant','Mean global ipTM','Δ vs WT','Median d (Å)','Median angle','NAC'],
                        cellLoc='center',colWidths=[.15,.19,.17,.19,.2,.1],bbox=[0,0,1,1])
        table.auto_set_font_size(False);table.set_fontsize(9.5)
        for (i,j),cell in table.get_celld().items():
            cell.set_edgecolor('white')
            if i==0:
                cell.set_facecolor('#18334d');cell.set_text_props(color='white',weight='bold')
            else:cell.set_facecolor('#edf2f6' if i%2 else '#fafbfc')
    fig.suptitle('The selected 20 mutants vs RFD4: geometry and global ipTM',fontsize=20,weight='bold',y=.97)
    fig.text(.5,.92,'User-confirmed October 15 testing panel · all 50 saved full-complex Boltz-2 poses per mutant',ha='center',fontsize=12)
    fig.text(.06,.08,'Blue: selected-mutant medians; gray: WT; whiskers: middle 50% of poses. Orange stars: 8 unique RFD4 structures.\n'
             'All RFD4 groups G1–G8 are labelled in gutters; connector lines point to exact coordinates. No coordinate jitter.\n'
             'NAC: distance ≤4.0 Å AND angle ≥165°. Median coordinates need not represent a single pose.\n'
             'Global ipTM is native confidence JSON [iptm]. RFD4 ligand atoms carry coordinate-conditioning annotations; Boltz-2 global ipTM is unavailable.',
             fontsize=10,color='#44546a',linespacing=1.7)
    savefig(fig,out,'rfd4_vs_selected20')
    fig,axes=plt.subplots(1,2,figsize=(13,8))
    fig.subplots_adjust(left=.1,right=.98,top=.85,bottom=.15,wspace=.27)
    rng=np.random.default_rng(0)
    for ax,arm in zip(axes,['atherm','haur']):
        rs=[r for r in summary if r['arm']==arm and r['mutation']!='WT']
        rs.sort(key=lambda r:-r['global_iptm_mean'])
        wt=next(r for r in summary if r['arm']==arm and r['mutation']=='WT')
        ax.axvline(wt['global_iptm_mean'],ls='--',color='#777',label=f"WT mean {wt['global_iptm_mean']:.5f}")
        for y,r in enumerate(rs):
            xs=[p['global_iptm'] for p in poses if p['id']==r['id']]
            ax.scatter(xs,y+rng.uniform(-.15,.15,len(xs)),s=10,color='#0072b2',alpha=.25)
            ax.errorbar(r['global_iptm_mean'],y,xerr=r['global_iptm_SD'],fmt='o',ms=5,color='#0072b2',capsize=3)
        ax.set(yticks=range(len(rs)),yticklabels=[r['mutation'] for r in rs],xlabel='Full-complex native global ipTM')
        ax.invert_yaxis();ax.spines[['top','right']].set_visible(False)
        ax.set_title('At' if arm=='atherm' else 'Ha');ax.legend(fontsize=10,frameon=False)
    fig.suptitle('Boltz-2 global ipTM for the experimental testing panel',fontsize=19,weight='bold',y=.95)
    fig.text(.1,.055,'Small dots: 50 saved poses; large dot: mean; whiskers: ±1 sample SD.\n'
             'Spread is across saved diffusion samples, not biological replicates or a measure of seed robustness.',fontsize=10,color='#44546a',linespacing=1.6)
    savefig(fig,out,'selected20_boltz2_global_iptm')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--project-root',type=Path,required=True)
    ap.add_argument('--plots-only',action='store_true',help='Render the previously audited CSVs without rescoring unchanged structures')
    args=ap.parse_args();root=args.project_root.resolve();sat=root/'saturation_mutagenesis'
    package=sat/'proposed_candidate_final/titer_final_20261015';out=root/'RFD4-proteina-designs/analysis'
    if args.plots_only:
        summary=read(out/'selected20_boltz2_full_summary.csv',',')
        strings={'id','arm','mutation','selection_group','sequence_sha256'}
        for row in summary:
            for key,value in row.items():
                if key not in strings and value!='':
                    row[key]=int(value) if key in {'n_models','NAC_count'} else float(value)
        poses=read(out/'selected20_boltz2_full_per_model.csv',',')
        for row in poses:row['global_iptm']=float(row['global_iptm'])
        rfd=read(out/'rfd4_unique_geometry_reference.csv',',')
        for row in rfd:
            row['distance_A']=float(row['distance_A']);row['angle_deg']=float(row['angle_deg'])
        assert [r['id'] for r in summary]==EXPECTED+['atherm_WT','haur_WT']
        assert len(poses)==1100 and len(rfd)==8
        plots(out,summary,poses,rfd)
        print('Re-rendered audited plots with every RFD4 group labelled.')
        return
    panel=read(package/'results/screening_panel.tsv')
    assert [r['id'] for r in panel]==EXPECTED, 'Panel differs from user-specified IDs'
    seqs={r.id:str(r.seq) for r in SeqIO.parse(package/'results/screening_panel.fasta','fasta')}
    assert seqs=={r['id']:r['sequence'] for r in panel}
    for r in SeqIO.parse(package/'results/wt_controls.fasta','fasta'):seqs[r.id]=str(r.seq)
    selected=EXPECTED+['atherm_WT','haur_WT']
    old={(r['id'],int(r['model'])):r for r in read(sat/'results/tables/nac_per_model.tsv') if r['id'] in selected}
    scores={(r['id'],int(r['sample'])):r for r in read(sat/'results/tables/models_long.tsv')
            if r['mode']=='full' and r['id'] in selected}
    poses=[]
    for vid in selected:
        candidates=[package/'audit/reopened'/f'{vid}.json',
                    sat/'proposed_candidate_final/archive/pre_titer_selection_20261015/audit/ensembles'/f'{vid}.json']
        cache=json.loads(next(p for p in candidates if p.exists()).read_text())
        assert len(cache['source_hashes'])==50
        for s in sorted(cache['source_hashes'],key=lambda s:s['model']):
            model=int(s['model']);path=sat/s['path']
            d,a,digest=audit_pose(path,seqs[vid]);assert digest==s['sha256']
            previous=old[vid,model]
            assert abs(d-float(previous['d_att']))<=.00051 and abs(a-float(previous['ang_att']))<=.0051
            nac=int(d<=4 and a>=165);assert nac==int(previous['nac'])
            confidence=path.with_name('confidence_'+path.stem+'.json')
            data=confidence.read_bytes();c=json.loads(data);iptm=float(c['iptm'])
            assert 0<=iptm<=1 and abs(iptm-float(scores[vid,model]['iptm']))<1e-12
            poses.append(dict(id=vid,arm=vid.split('_')[0],model=model,global_iptm=iptm,
                distance_A=d,angle_deg=a,NAC=nac,atom_checks_pass=True,sequence_verified=True,
                structure_path=str(path.relative_to(root)),structure_sha256=digest,
                confidence_path=str(confidence.relative_to(root)),confidence_sha256=sha(data)))
        print(f'Verified {vid}: 50 raw structures and native global ipTM JSONs',flush=True)
    assert len(poses)==1100 and len({(p['id'],p['model']) for p in poses})==1100
    means={vid:float(np.mean([p['global_iptm'] for p in poses if p['id']==vid])) for vid in selected}
    summary=[];panelmap={r['id']:r for r in panel}
    for vid in selected:
        rs=[p for p in poses if p['id']==vid];arm,mutation=vid.split('_',1)
        d=np.array([p['distance_A'] for p in rs]);a=np.array([p['angle_deg'] for p in rs]);v=np.array([p['global_iptm'] for p in rs])
        row=dict(id=vid,arm=arm,mutation=mutation,n_models=len(rs),global_iptm_mean=means[vid],
            global_iptm_SD=float(np.std(v,ddof=1)),global_iptm_median=float(np.median(v)),
            global_iptm_min=float(v.min()),global_iptm_max=float(v.max()),
            global_iptm_WT_mean=means[arm+'_WT'],global_iptm_delta_vs_WT=means[vid]-means[arm+'_WT'],
            distance_mean_A=float(d.mean()),distance_median_A=float(np.median(d)),
            distance_q25_A=float(np.quantile(d,.25)),distance_q75_A=float(np.quantile(d,.75)),
            angle_mean_deg=float(a.mean()),angle_median_deg=float(np.median(a)),
            angle_q25_deg=float(np.quantile(a,.25)),angle_q75_deg=float(np.quantile(a,.75)),
            NAC_count=sum(p['NAC'] for p in rs),NAC_fraction=sum(p['NAC'] for p in rs)/len(rs),
            selection_group=panelmap.get(vid,{}).get('testing_priority','WT control'),
            ESM_MVAP_global_delta_vs_WT=panelmap.get(vid,{}).get('esm_delta',''),
            sequence_sha256=sha(seqs[vid].encode()))
        if vid in panelmap:assert abs(row['global_iptm_delta_vs_WT']-float(panelmap[vid]['boltz_full_delta']))<1e-10
        summary.append(row)
    new=read(out/'attack_geometry.csv',',');groups=defaultdict(list)
    for r in new:groups[r['geometry_group']].append(r)
    rfd=[]
    for group,members in sorted(groups.items()):
        row,audit=measure(root/'RFD4-proteina-designs'/members[0]['file'])
        assert row['sha256']==members[0]['sha256']
        rfd.append(dict(group=group,candidate_ids=';'.join(r['candidate'] for r in members),
                        distance_A=row['distance_A'],angle_deg=row['angle_deg'],
                        global_iptm='',global_iptm_status='No Boltz-2 full-complex result located in supplied RFD4 inputs'))
    write_csv(out/'selected20_boltz2_full_summary.csv',summary)
    write_csv(out/'selected20_boltz2_full_per_model.csv',poses)
    write_csv(out/'rfd4_unique_geometry_reference.csv',rfd)
    plots(out,summary,poses,rfd)
    provenance=dict(panel_path=str((package/'results/screening_panel.tsv').relative_to(root)),
        panel_sha256=sha((package/'results/screening_panel.tsv').read_bytes()),confirmed_ids=EXPECTED,
        verified_raw_structures=1100,verified_confidence_JSONs=1100,models_per_sequence=50,
        global_metric='Native top-level confidence JSON iptm; no ps_iptm, pair ipTM, pTM or confidence_score substitution',
        previous_atoms='MVAP chain L O16; ATP chain T PG; beta-gamma bridge chain T O3B; vertex PG',
        RFD4_atoms='MVAP IMV O16; ATP IAT P41; beta-gamma bridge IAT O27; vertex P41',
        selection_policy=json.loads((package/'selection_policy.json').read_text()),
        limitations=['RFD4 inputs are coordinate-conditioned generated structures, not Boltz-2 cofold predictions.',
                     'The panel primary metric was ESMFold2 MVAP-only global ipTM; Boltz full-complex global ipTM was secondary.',
                     'Pose spread is conditional on saved campaign sampling; no fresh independent-seed confirmation.'])
    (out/'selected20_comparison_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Saved panel comparison, global ipTM plot, 22-sequence summary, and 1100 audited per-model records.')


if __name__=='__main__':
    main()
