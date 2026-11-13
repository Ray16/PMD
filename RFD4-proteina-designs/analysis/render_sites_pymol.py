#!/usr/bin/env python3
"""Render all eight active sites with identical camera settings using CPU PyMOL."""
import argparse
import json
from pathlib import Path
import pymol2


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--site-dir',type=Path,required=True)
    folder=ap.parse_args().site_dir;records=json.loads((folder/'active_site_data.json').read_text())
    images=folder/'images';images.mkdir(exist_ok=True)
    common=None
    with pymol2.PyMOL() as pymol:
        c=pymol.cmd
        for record in records:
            g=record['group'];c.reinitialize();c.set('max_threads',1);c.set('ray_trace_mode',1)
            c.set('antialias',1);c.set('ray_shadows',0);c.set('specular',.15);c.bg_color('white')
            c.load(str(folder/'structures'/f'{g}.pdb'),'design')
            c.hide('everything','all');c.remove('hydro')
            c.set_color('mvap_green',[.13,.64,.40]);c.set_color('atp_gold',[.88,.62,.14])
            c.set_color('base_purple',[.48,.24,.70]);c.set_color('gamma_orange',[.84,.25,.02])
            c.set_color('metal_teal',[.1,.68,.67]);c.set_color('lys_blue',[.23,.43,.70])
            c.set_color('pocket_gray',[.70,.74,.78])
            # Enforce the embedded chemical graph for both nonstandard ligands.
            c.unbond('resn IMV or resn IAT','resn IMV or resn IAT')
            for first,second,order in record['bonds']:
                c.bond(f'design and id {first}',f'design and id {second}',order)
            c.select('ligands','resn IMV or resn IAT')
            c.select('pocket','byres (chain H within 6 of ligands)')
            c.show('cartoon','pocket');c.color('pocket_gray','pocket');c.set('cartoon_transparency',.72)
            c.show('sticks','ligands');c.set('stick_radius',.16,'ligands')
            c.color('mvap_green','resn IMV and elem C');c.color('atp_gold','resn IAT and elem C')
            c.color('oxygen','ligands and elem O');c.color('nitrogen','ligands and elem N');c.color('atp_gold','elem P')
            display='+'.join(str(i) for i in record['display_residues'])
            c.select('display_residues',f'chain H and resi {display}')
            c.show('sticks','display_residues');c.set('stick_radius',.13,'display_residues')
            c.color('pocket_gray','display_residues and elem C')
            c.color('oxygen','display_residues and elem O');c.color('nitrogen','display_residues and elem N')
            def atom(a):return f"design and chain {a['chain']} and resi {a['resi']} and name {a['name']}"
            base=record['Asp'];lys=record['Lys'];oh=record['hydroxyl'];gamma=record['gamma_acceptor']
            c.color('base_purple',f"chain H and resi {base['resi']} and elem C")
            c.color('lys_blue',f"chain H and resi {lys['resi']} and elem C")
            c.show('spheres','elem Mg');c.set('sphere_scale',.38,'elem Mg');c.color('metal_teal','elem Mg')
            c.show('spheres',f"({atom(oh)}) or ({atom(gamma)}) or ({atom(base)})")
            c.set('sphere_scale',.22,f"({atom(oh)}) or ({atom(gamma)}) or ({atom(base)})")
            c.distance('Asp_route',atom(oh),atom(base),cutoff=100);c.set('dash_color','base_purple','Asp_route')
            c.distance('ATP_route',atom(oh),atom(gamma),cutoff=100);c.set('dash_color','gamma_orange','ATP_route')
            c.hide('labels','all');c.set('dash_radius',.055);c.set('dash_gap',.18)
            c.set('label_color','black');c.set('label_size',15);c.set('label_outline_color','white')
            c.label(f"chain H and resi {base['resi']} and name CA",f'"Asp{base["resi"]}"')
            c.label(f"chain H and resi {lys['resi']} and name CA",f'"Lys{lys["resi"]}"')
            for m in record['metals']:
                c.label(atom(m['atom']),f'"Mg {m["atom"]["chain"]}"')
            # Same view for every group after a rigid-body ligand alignment.
            if common is None:
                c.orient('ligands');c.rotate('x',25);c.rotate('z',-15)
                c.zoom('ligands or (chain H and resi '+str(base['resi'])+')',buffer=3.0)
                common=c.get_view()
            else:c.set_view(common)
            c.png(str(images/f'{g}_active_site.png'),width=1400,height=1050,dpi=160,ray=1,quiet=1)
            c.save(str(images/f'{g}_active_site.pse'))
            # A second view emphasizes the two acceptor contacts, not the entire ATP.
            c.hide('cartoon');c.hide('everything','resn IAT and not (name P41+P43+O13+O16+O17+O27+O15+O19+O28)')
            c.zoom(f"({atom(oh)}) or ({atom(base)}) or ({atom(record['gamma_P'])}) or ({atom(gamma)})",buffer=4.0)
            c.png(str(images/f'{g}_proton_routes.png'),width=1200,height=1000,dpi=160,ray=1,quiet=1)
            c.save(str(images/f'{g}_proton_routes.pse'))
            print(f'{g}: active-site and proton-route images rendered',flush=True)


if __name__=='__main__':main()
