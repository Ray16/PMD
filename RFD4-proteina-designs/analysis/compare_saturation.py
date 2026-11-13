#!/usr/bin/env python3
"""Compare generated RFD4 geometry with the historical Boltz-2 saturation screen.

python compare_saturation.py --project-root /homes/rzhu/PMD
Uses previously audited RFD4 measurements and historical per-pose NAC data.
Rechecks source hashes, selected raw WT poses and ligand-conditioning annotations.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from adjustText import adjust_text
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def read(path, delimiter=','):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_wt(root, lookup):
    audits = []
    base = root / 'saturation_mutagenesis'
    for arm in ['atherm', 'haur']:
        cache = json.loads((base / 'proposed_candidate_final/archive/pre_titer_selection_20261015'
                            / 'audit/ensembles' / f'{arm}_WT.json').read_text())
        for model in [0, 49]:
            source = next(s for s in cache['source_hashes'] if s['model'] == model)
            path = base / source['path']
            assert sha(path) == source['sha256'], 'Historical raw structure hash changed'
            lig = {}
            for line in path.read_text().splitlines():
                if line.startswith('HETATM'):
                    key = (line[21], line[12:16].strip())
                    assert key not in lig, 'Ambiguous ligand atom'
                    lig[key] = np.array([float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]])
            # Recheck the chemical neighbourhood of the historical atom mapping.
            for a, b in [(('L', 'O16'), ('L', 'C23')), (('T', 'PG'), ('T', 'O3B')),
                         (('T', 'PB'), ('T', 'O3B')), (('T', 'PA'), ('T', 'O3A')),
                         (('T', 'PB'), ('T', 'O3A'))]:
                assert 1.1 < np.linalg.norm(lig[a] - lig[b]) < 1.95, f'Unexpected bond {a}, {b}'
            assert {name for (chain, name), xyz in lig.items() if chain == 'L' and
                    name.startswith('C') and np.linalg.norm(xyz-lig['L', 'O16']) < 1.9} == {'C23'}
            assert {name for (chain, name), xyz in lig.items() if chain == 'L' and
                    name.startswith('C') and name != 'C23' and
                    np.linalg.norm(xyz-lig['L', 'C23']) < 1.9} == {'C20', 'C21', 'C22'}
            assert {name for (chain, name), xyz in lig.items() if chain == 'T' and
                    name.startswith('P') and np.linalg.norm(xyz-lig['T', 'O3B']) < 1.95} == {'PB', 'PG'}
            o, p, b = lig['L', 'O16'], lig['T', 'PG'], lig['T', 'O3B']
            u, v = o-p, b-p
            d = float(np.linalg.norm(u))
            angle = math.degrees(math.acos(float(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v), -1, 1))))
            stored = lookup[f'{arm}_WT', model]
            assert abs(d-float(stored['d_att'])) <= .00051
            assert abs(angle-float(stored['ang_att'])) <= .0051
            assert int(d <= 4 and angle >= 165) == int(stored['nac'])
            audits.append(dict(file=str(path.relative_to(root)), sha256=sha(path),
                               distance_A=d, angle_deg=angle, historical_table_matches=True))
    return audits


def summary(label, rows):
    d = np.array([r['d'] for r in rows]); a = np.array([r['a'] for r in rows])
    return dict(cohort=label, n_structures=len(rows), n_sequences=len({r['id'] for r in rows}),
                distance_median_A=float(np.median(d)), distance_q10_A=float(np.quantile(d,.1)),
                distance_q90_A=float(np.quantile(d,.9)), distance_min_A=float(d.min()), distance_max_A=float(d.max()),
                angle_median_deg=float(np.median(a)), angle_q10_deg=float(np.quantile(a,.1)),
                angle_q90_deg=float(np.quantile(a,.9)), angle_min_deg=float(a.min()), angle_max_deg=float(a.max()),
                NAC_fraction=sum(r['nac'] for r in rows)/len(rows))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,required=True)
    args=parser.parse_args(); root=args.project_root.resolve()
    out=root/'RFD4-proteina-designs/analysis'
    source=root/'saturation_mutagenesis/results/tables/nac_per_model.tsv'
    raw=read(source,'\t'); lookup={(r['id'],int(r['model'])):r for r in raw}
    assert len(lookup)==len(raw), 'Duplicate historical variant/model'
    counts=Counter(r['id'] for r in raw)
    assert set(counts.values())=={50}, 'Historical ensembles are not all n=50'
    old=[dict(id=r['id'],arm=r['arm'],d=float(r['d_att']),a=float(r['ang_att']),
              nac=int(r['nac'])) for r in raw]
    new_raw=read(out/'attack_geometry.csv')
    new=[dict(id=r['candidate'],d=float(r['distance_A']),a=float(r['angle_deg']),
              nac=int(r['near_attack']=='True'),group=r['geometry_group']) for r in new_raw]
    flags=[]; coordinate_groups={}; sequence_hashes=set()
    for r in new_raw:
        path=root/'RFD4-proteina-designs'/r['file']
        assert sha(path)==r['sha256'], 'RFD4 input changed since atom audit'
        cif=MMCIF2Dict(str(path))
        sequence=''.join(''.join(cif['_entity_poly.pdbx_seq_one_letter_code_can']).split())
        sequence_hashes.add(hashlib.sha256(sequence.encode()).hexdigest())
        signature=list(zip(*(cif['_atom_site.'+key] for key in
            ['label_asym_id','label_seq_id','label_comp_id','label_atom_id','Cartn_x','Cartn_y','Cartn_z'])))
        coordinate_hash=hashlib.sha256(repr(signature).encode()).hexdigest()
        coordinate_groups.setdefault(coordinate_hash,[]).append(r['file'])
        for comp,name in [('IMV','O16'),('IAT','P41'),('IAT','O27')]:
            hits=[i for i,(c,a) in enumerate(zip(cif['_atom_site.label_comp_id'],cif['_atom_site.label_atom_id'])) if (c,a)==(comp,name)]
            assert len(hits)==1
            region=cif['_atom_site.annotation_src-coordinate-conditioned-region_1_atom'][hits[0]]
            assert region not in ['.', '?', ''], 'Missing coordinate-conditioned-region annotation'
        flags.append(r['file'])
    audits=audit_wt(root,lookup)
    unique=list({r['group']:r for r in new}.values())
    assert len(unique)==len(coordinate_groups)==len(sequence_hashes)==8
    cohorts=[summary('RFD4 unique generated structures',unique)]
    for arm in ['atherm','haur']:
        for wt in [True,False]:
            rs=[r for r in old if r['arm']==arm and r['id'].endswith('_WT')==wt]
            cohorts.append(summary(f'Boltz-2 {arm} '+('WT' if wt else 'all mutants'),rs))
    with (out/'saturation_comparison_summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(cohorts[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(cohorts)
    fig,axes=plt.subplots(2,2,figsize=(13,10.5))
    fig.subplots_adjust(top=.84,bottom=.17,hspace=.35,wspace=.2)
    for col,(arm,name) in enumerate([('atherm','At'),('haur','Ha')]):
        mutants=[r for r in old if r['arm']==arm and not r['id'].endswith('_WT')]
        wt=[r for r in old if r['id']==arm+'_WT']
        for j,ax in enumerate(axes[:,col]):
            ax.add_patch(Rectangle((165,0),15,4,color='#edf7ef',zorder=0))
            ax.axvline(165,color='#658b68',ls='--',lw=1)
            ax.axhline(4,color='#658b68',ls='--',lw=1)
            ax.scatter([r['a'] for r in mutants],[r['d'] for r in mutants],s=3,
                       color='#6e7d89',alpha=.12,rasterized=True,zorder=1)
            ax.scatter([r['a'] for r in wt],[r['d'] for r in wt],s=24,
                       facecolors='none',edgecolors='#0072b2',linewidths=.8,label='Previous WT (50 poses)',zorder=2)
            ax.scatter([r['a'] for r in unique],[r['d'] for r in unique],s=95,
                       marker='*',color='#d55e00',edgecolor='white',linewidths=.6,
                       label='RFD4 (8 distinct geometries)',zorder=4)
            ax.set(xlabel='Attack angle (°)',ylabel='Attack distance (Å)')
            ax.spines[['top','right']].set_visible(False)
            if j==0:
                ax.set(xlim=(0,181),ylim=(2.5,14))
                ax.set_title(f'{name}: {len(mutants)//50:,} mutants × 50 poses',fontsize=12)
                ax.text(.03,.96,f"WT NAC: {sum(r['nac'] for r in wt)}/50\n"
                        f"All mutant poses: {100*sum(r['nac'] for r in mutants)/len(mutants):.1f}% NAC",
                        transform=ax.transAxes,va='top',fontsize=11)
            else:
                ax.set(xlim=(150,181),ylim=(3.3,4.5))
                ax.set_title('Zoom around the near-attack region',fontsize=12)
                labels=[ax.text(r['a'],r['d'],r['group'],fontsize=8,color='#803800') for r in unique]
                adjust_text(labels,ax=ax,expand=(1.5,1.5),
                            arrowprops=dict(arrowstyle='-',color='#803800',lw=.5))
                ax.text(.03,.97,'NAC: angle ≥165° AND distance ≤4.0 Å',transform=ax.transAxes,va='top',fontsize=9)
        axes[0,col].legend(loc='upper right',fontsize=8,frameon=False)
    fig.suptitle('RFD4 geometry compared with saturation mutagenesis',fontsize=19,weight='bold',y=.96)
    fig.text(.5,.9,'Same chemical atom definitions · historical Boltz-2 full complexes · all raw poses shown in top panels',ha='center',fontsize=11)
    fig.text(.07,.08,'Gray: saturation mutants; blue: WT; orange: RFD4. Lower panels magnify part of the top panels.\n'
             'RFD4: 37 files contain 8 distinct complete structures (counted once each); historical variants: 50 cofold poses each.\n'
             'RFD4 ligand atoms are annotated as coordinate-conditioned regions; favorable generated geometry is not an independent cofold result.',
             fontsize=10,color='#44546a',linespacing=1.7)
    fig.savefig(out/'rfd4_vs_saturation.png',dpi=180,facecolor='white')
    fig.savefig(out/'rfd4_vs_saturation.svg',facecolor='white')
    svg=out/'rfd4_vs_saturation.svg';svg.write_text('\n'.join(l.rstrip() for l in svg.read_text().splitlines())+'\n')
    plt.close(fig)
    provenance=dict(historical_table=str(source.relative_to(root)),historical_table_sha256=sha(source),
        historical_total_poses=len(raw),historical_total_sequences=len(counts),models_per_sequence=50,
        raw_WT_spot_checks=audits,RFD4_coordinate_conditioned_files=flags,
        RFD4_unique_sequences=len(sequence_hashes),RFD4_unique_complete_coordinate_sets=len(coordinate_groups),
        RFD4_duplicate_coordinate_groups=coordinate_groups,RFD4_summary_weighting='Each distinct complete coordinate set counted once',
        atom_mapping={'historical':'chain L O16 / chain T PG / chain T O3B',
                      'RFD4':'IMV O16 / IAT P41 / IAT O27'},
        limitation='RFD4 generated coordinate-conditioned inputs and independent 50-pose cofold ensembles are not equivalent evidence.')
    (out/'saturation_comparison_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(cohorts,indent=2))


if __name__=='__main__':
    main()
