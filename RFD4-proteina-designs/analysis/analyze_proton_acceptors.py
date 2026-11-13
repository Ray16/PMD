#!/usr/bin/env python3
"""Measure both candidate proton-acceptor contacts and prepare active-site views.

All distances are heavy-atom proxies, not proton-transfer barriers or probabilities.
Run in the pmd environment; render_sites_pymol.py renders the prepared structures.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def read(path):
    return list(csv.DictReader(path.open()))


def label(a):
    return f"{a['chain']}:{a['resn']}{a['resi']}:{a['name']}"


def distance(a,b):
    return float(np.linalg.norm(np.array(a['xyz'])-np.array(b['xyz'])))


def best(atoms,targets):
    return min((distance(a,t),a,t) for a in atoms for t in targets)


def plot_all_acceptors(rows,out):
    series=[('Asp_O_to_MVAP_O16_A','Nearest Asp O','#7b3fb2'),
            ('Glu_O_to_MVAP_O16_A','Nearest Glu O','#168249'),
            ('His_N_to_MVAP_O16_A','Nearest His N','#2071b5'),
            ('gamma_O_to_MVAP_O16_A','Nearest ATP γ O','#d55e00')]
    fig=plt.figure(figsize=(17,11));grid=fig.add_gridspec(2,2,height_ratios=[1.25,1])
    full=fig.add_subplot(grid[0,0]);zoom=fig.add_subplot(grid[0,1]);tab=fig.add_subplot(grid[1,:])
    fig.subplots_adjust(left=.065,right=.98,bottom=.15,top=.84,hspace=.38,wspace=.2)
    x=np.arange(len(rows));offsets=[-.24,-.08,.08,.24]
    for (key,name,color),offset in zip(series,offsets):
        full.scatter(x+offset,[float(r[key]) for r in rows],s=48,color=color,label=name,zorder=3)
    full.set(xticks=x,xticklabels=[r['group'] for r in rows],ylabel='Distance from MVAP hydroxyl IMV:O16 (Å)',xlabel='RFD4 group',xlim=(-.65,7.65))
    full.set_ylim(2.9,max(float(r[k]) for r in rows for k,_,_ in series)+1.0)
    full.set_title('All four candidate acceptor classes',fontsize=13,pad=12)
    full.legend(frameon=False,ncol=2,fontsize=9,loc='upper left')
    for i,r in enumerate(rows):
        a=float(r['Asp_O_to_MVAP_O16_A']);b=float(r['gamma_O_to_MVAP_O16_A'])
        zoom.plot([i-.14,i+.14],[a,b],color='#abb5bf',lw=1,zorder=1)
        zoom.scatter(i-.14,a,s=50,color='#7b3fb2',zorder=3)
        zoom.scatter(i+.14,b,s=50,color='#d55e00',zorder=3)
        zoom.annotate(f'{a:.2f}',(i-.14,a),xytext=(-4,8),textcoords='offset points',ha='right',fontsize=9,color='#7b3fb2')
        zoom.annotate(f'{b:.2f}',(i+.14,b),xytext=(4,-15),textcoords='offset points',ha='left',fontsize=9,color='#d55e00')
    zoom.set(xticks=x,xticklabels=[r['group'] for r in rows],xlabel='RFD4 group',ylabel='O···O distance (Å)',xlim=(-.65,7.65))
    values=[float(r[k]) for r in rows for k in ['Asp_O_to_MVAP_O16_A','gamma_O_to_MVAP_O16_A']]
    zoom.set_ylim(min(values)-.14,max(values)+.15);zoom.set_title('Zoom: Asp versus ATP γ oxygen',fontsize=13,pad=12)
    for ax in [full,zoom]:ax.grid(axis='y',alpha=.18);ax.spines[['top','right']].set_visible(False)
    tab.axis('off');cells=[]
    for r in rows:
        short=lambda key:r[key].split(':',1)[1]
        cells.append([r['group'],short('Asp_atom'),f"{float(r['Asp_O_to_MVAP_O16_A']):.3f}",
            short('nearest_Glu_to_OH_atom'),f"{float(r['Glu_O_to_MVAP_O16_A']):.3f}",
            short('nearest_His_to_OH_atom'),f"{float(r['His_N_to_MVAP_O16_A']):.3f}",
            short('gamma_acceptor_atom'),f"{float(r['gamma_O_to_MVAP_O16_A']):.3f}"])
    table=tab.table(cellText=cells,colLabels=['Group','Asp O atom','Å','Glu O atom','Å','His N atom','Å','ATP γ O atom','Å'],
        cellLoc='center',colWidths=[.06,.14,.075,.14,.075,.14,.075,.15,.075],bbox=[0,0,1,1])
    table.auto_set_font_size(False);table.set_fontsize(10.5)
    for (i,j),cell in table.get_celld().items():
        cell.set_edgecolor('white')
        if i==0:cell.set_facecolor('#18334d');cell.set_text_props(color='white',weight='bold')
        else:cell.set_facecolor('#edf2f6' if i%2 else '#fafbfc')
    fig.suptitle('Candidate proton acceptors in all eight RFD4 structures',fontsize=21,weight='bold',y=.97)
    fig.text(.5,.915,'Asp carboxylate O · Glu carboxylate O · His ring N · ATP nonbridging γ-phosphate O',ha='center',fontsize=12)
    fig.text(.065,.075,'Same donor atom in every measurement: MVAP IMV:O16. Each class reports its nearest eligible atom.\n'
        'Histidine protonation is unresolved: its nearest ring N is not automatically an available proton acceptor.\n'
        'Distances assess geometric access, not mechanism probabilities or reaction barriers. Water-mediated routes are not resolved.',
        fontsize=10.5,color='#44546a',linespacing=1.7)
    for ext in ['png','svg']:fig.savefig(out/f'proton_acceptor_distances.{ext}',dpi=180,facecolor='white')
    svg=out/'proton_acceptor_distances.svg';svg.write_text('\n'.join(l.rstrip() for l in svg.read_text().splitlines())+'\n');plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--project-root',type=Path,required=True)
    ap.add_argument('--plots-only',action='store_true')
    args=ap.parse_args();root=args.project_root.resolve();out=root/'RFD4-proteina-designs/analysis'
    if args.plots_only:
        plot_all_acceptors(read(out/'proton_acceptor_distances.csv'),out)
        print('Plotted Asp, Glu, His and ATP gamma-oxygen distances.');return
    sites=out/'active_sites';(sites/'structures').mkdir(parents=True,exist_ok=True)
    refs=read(out/'rfd4_unique_geometry_reference.csv');rows=[];records=[];reference=None
    original={r['file']:r for r in read(out/'attack_geometry.csv')}
    for ref in refs:
        filename=ref['candidate_ids'].split(';')[0]+'_generated.cif';path=root/'RFD4-proteina-designs'/filename
        assert hashlib.sha256(path.read_bytes()).hexdigest()==original[filename]['sha256']
        c=MMCIF2Dict(str(path));atoms=[]
        for i,name in enumerate(c['_atom_site.label_atom_id']):
            atoms.append(dict(serial=i+1,name=name,resn=c['_atom_site.label_comp_id'][i],
                chain=c['_atom_site.auth_asym_id'][i],resi=int(c['_atom_site.auth_seq_id'][i]),
                element=c['_atom_site.type_symbol'][i].upper(),
                xyz=[float(c['_atom_site.Cartn_'+a][i]) for a in 'xyz'],
                conditioned_region=c['_atom_site.annotation_src-coordinate-conditioned-region_1_atom'][i]))
        def one(comp,name):
            hits=[a for a in atoms if a['resn']==comp and a['name']==name];assert len(hits)==1
            return hits[0]
        hydroxyl=one('IMV','O16');pg=one('IAT','P41');bridge=one('IAT','O27')
        gamma=[one('IAT',name) for name in ['O13','O16','O17']]
        mvap=[a for a in atoms if a['resn']=='IMV'];atp=[a for a in atoms if a['resn']=='IAT']
        protein=[a for a in atoms if a['chain']=='H']
        acids=[a for a in protein if a['resn'] in ['ASP','GLU'] and a['name'] in ['OD1','OD2','OE1','OE2']]
        asps=[a for a in acids if a['resn']=='ASP'];lys=[a for a in protein if a['resn']=='LYS' and a['name']=='NZ']
        args=[a for a in protein if a['resn']=='ARG' and a['name'] in ['NE','NH1','NH2']]
        his=[a for a in protein if a['resn']=='HIS' and a['name'] in ['ND1','NE2']]
        dbase,base,_=best(asps,[hydroxyl]);dgamma,gamma_acceptor,_=best(gamma,[hydroxyl])
        carbox=[one('IMV',name) for name in ['O11','O13']]
        phosphate=[one('IMV',name) for name in ['O12','O14','O15','O17']]
        dlys,lysine,lys_target=best(lys,gamma)
        dcat,cation,cat_target=best(lys+args,carbox)
        dhis,histidine,his_target=best(his,phosphate)
        dglu,glutamate,_=best([a for a in acids if a['resn']=='GLU'],[hydroxyl])
        dhisOH,hisOH,_=best(his,[hydroxyl])
        mg=[a for a in atoms if a['element']=='MG'];assert len(mg)==2
        metal=[]
        for m in mg:
            oxygen=[a for a in atoms if a['element']=='O']
            contacts=sorted([(distance(m,a),a['serial']) for a in oxygen])[:8]
            metal.append(dict(atom=m,nearest_O=[dict(distance_A=d,atom=atoms[j-1]) for d,j in contacts],
                              nearest_gamma_O_A=min(distance(m,a) for a in gamma),
                              oxygen_count_within_2p6A=sum(distance(m,a)<=2.6 for a in oxygen)))
        row=dict(group=ref['group'],representative=filename,candidate_ids=ref['candidate_ids'],
            Asp_atom=label(base),Asp_O_to_MVAP_O16_A=dbase,
            Asp_OD1_to_MVAP_O16_A=distance(next(a for a in protein if a['resi']==base['resi'] and a['name']=='OD1'),hydroxyl),
            Asp_OD2_to_MVAP_O16_A=distance(next(a for a in protein if a['resi']==base['resi'] and a['name']=='OD2'),hydroxyl),
            gamma_acceptor_atom=label(gamma_acceptor),gamma_O_to_MVAP_O16_A=dgamma,
            **{f'IAT_{a["name"]}_to_MVAP_O16_A':distance(a,hydroxyl) for a in gamma},
            Asp_minus_gamma_distance_A=dbase-dgamma,
            attack_O_to_Pgamma_A=distance(hydroxyl,pg),
            attack_angle_deg=float(original[filename]['angle_deg']),
            nearest_Lys_atom=label(lysine),Lys_to_gamma_O_A=dlys,
            nearest_carboxylate_cation=label(cation),cation_to_MVAP_carboxylate_O_A=dcat,
            nearest_His_atom=label(histidine),His_to_MVAP_phosphate_O_A=dhis,
            nearest_Glu_to_OH_atom=label(glutamate),Glu_O_to_MVAP_O16_A=dglu,
            nearest_His_to_OH_atom=label(hisOH),His_N_to_MVAP_O16_A=dhisOH,
            Mg_M_to_gamma_O_A=metal[0]['nearest_gamma_O_A'],Mg_N_to_gamma_O_A=metal[1]['nearest_gamma_O_A'],
            Asp_coordinate_annotation=base['conditioned_region'],Lys_coordinate_annotation=lysine['conditioned_region'])
        rows.append(row)
        # Align every structure on the same 45 ligand atoms. Distances are unchanged.
        ligand=mvap+atp;keys=sorted((a['resn'],a['name']) for a in ligand)
        bykey={(a['resn'],a['name']):a for a in ligand};x=np.array([bykey[k]['xyz'] for k in keys])
        if reference is None:reference=x.copy()
        xc=x.mean(0);yc=reference.mean(0);u,_,vt=np.linalg.svd((x-xc).T@(reference-yc))
        rotation=u@vt
        if np.linalg.det(rotation)<0:u[:,-1]*=-1;rotation=u@vt
        for a in atoms:a['xyz']=((np.array(a['xyz'])-xc)@rotation+yc).tolist()
        assert abs(distance(base,hydroxyl)-dbase)<1e-10
        display={base['resi'],lysine['resi']}
        if dcat<=4.5:display.add(cation['resi'])
        if dhis<=4.5:display.add(histidine['resi'])
        # Include additional polar residues only when they contact either putative route.
        for a in protein:
            if a['element'] in ['N','O'] and a['name'] not in ['N','O'] and min(distance(a,t) for t in [hydroxyl,*gamma])<3.5:
                display.add(a['resi'])
        pdb=[]
        for a in atoms:
            rec='ATOM  ' if a['chain']=='H' else 'HETATM'
            x,y,z=a['xyz'];name=f" {a['name']:<3}" if len(a['name'])<4 else a['name']
            pdb.append(f"{rec}{a['serial']:5d} {name:4s} {a['resn']:3s} {a['chain']:1s}{a['resi']:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{20.:6.2f}          {a['element']:>2s}")
        bonds=[]
        for comp,a,b,order in zip(*(c['_chem_comp_bond.'+k] for k in ['comp_id','atom_id_1','atom_id_2','value_order'])):
            if comp in ['IMV','IAT']:
                aa=one(comp,a);bb=one(comp,b);bonds.append([aa['serial'],bb['serial'],2 if order=='DOUB' else 1])
                pdb.append(f"CONECT{aa['serial']:5d}{bb['serial']:5d}")
        pdb.extend(['END','']);pdbtext='\n'.join(pdb)
        (sites/'structures'/f"{ref['group']}.pdb").write_text(pdbtext)
        record=dict(group=ref['group'],candidate_ids=ref['candidate_ids'],metrics=row,pdb=pdbtext,
            hydroxyl=hydroxyl,gamma_P=pg,bridge=bridge,Asp=base,gamma_acceptor=gamma_acceptor,
            Lys=lysine,carboxylate_cation=cation,histidine=histidine,gamma_O=gamma,metals=metal,
            display_residues=sorted(display),bonds=bonds)
        records.append(record)
        print(ref['group'],f"Asp {dbase:.3f} Å; gamma O {dgamma:.3f} Å ({gamma_acceptor['name']}); cation {dcat:.3f} Å",flush=True)
    with (out/'proton_acceptor_distances.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    (sites/'active_site_data.json').write_text(json.dumps(records,indent=2)+'\n')
    plot_all_acceptors(rows,out)


if __name__=='__main__':main()
