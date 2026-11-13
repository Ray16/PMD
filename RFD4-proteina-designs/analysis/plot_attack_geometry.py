#!/usr/bin/env python3
"""Audit ligand atom identities and plot RFD4 attack geometry.

Run with Python containing biopython, numpy, matplotlib and adjustText:
  python plot_attack_geometry.py --input-dir RFD4-proteina-designs
The atom mapping is verified from each CIF's embedded chemical bond graph.
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
import numpy as np
from adjustText import adjust_text
from Bio.PDB import MMCIFParser, calc_angle
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def check(condition, message):
    if not condition:
        raise ValueError(message)


def measure(path):
    cif = MMCIF2Dict(str(path))
    for comp, expected in [('IMV', {'C': 6, 'O': 7, 'P': 1}),
                           ('IAT', {'C': 10, 'N': 5, 'O': 13, 'P': 3})]:
        elements = Counter(e.upper() for c, e in zip(cif['_atom_site.label_comp_id'],
                           cif['_atom_site.type_symbol']) if c == comp and e.upper() != 'H')
        check(elements == expected, f'{path.name}: unexpected {comp} heavy-atom formula')
    graphs = defaultdict(lambda: defaultdict(set))
    for comp, a, b in zip(*(cif['_chem_comp_bond.' + key] for key in
                            ['comp_id', 'atom_id_1', 'atom_id_2'])):
        graphs[comp][a].add(b)
        graphs[comp][b].add(a)
    mv, atp = graphs['IMV'], graphs['IAT']
    # MVAP C3 is bonded to three carbons and one terminal oxygen.
    check(mv['O16'] == {'C23'}, 'IMV O16 is not a terminal hydroxyl oxygen')
    check(mv['C23'] == {'C21', 'C22', 'C20', 'O16'}, 'Wrong MVAP tertiary carbon')
    check({'C18', 'C23'} <= mv['C21'] and {'C19', 'C23'} <= mv['C22'],
          'MVAP C3 backbone does not reach carboxyl and phosphorylated ends')
    check({'O11', 'O13'} <= mv['C18'] and mv['C20'] == {'C23'} and
          mv['O17'] == {'C19', 'P24'}, 'Wrong MVAP backbone or phosphate')
    # Trace from ribose to alpha, beta, gamma phosphate, independently of names.
    check(atp['O26'] == {'C36', 'P42'}, 'ATP alpha phosphate is not ribose-linked')
    check(atp['O28'] == {'P42', 'P43'}, 'ATP alpha-beta bridge is wrong')
    check(atp['O27'] == {'P43', 'P41'}, 'ATP beta-gamma bridge is wrong')
    check(atp['P41'] == {'O27', 'O13', 'O16', 'O17'}, 'Wrong terminal gamma phosphate')
    check(all(atp[o] == {'P41'} for o in ['O13', 'O16', 'O17']),
          'Gamma phosphate has another bridge')
    selections = [('IMV', 'O16', 'O'), ('IAT', 'P41', 'P'), ('IAT', 'O27', 'O')]
    coords, identifiers, selected_indices = [], [], []
    for comp, atom, element in selections:
        hits = [i for i, (c, a) in enumerate(zip(cif['_atom_site.label_comp_id'],
                                                cif['_atom_site.label_atom_id']))
                if c == comp and a == atom]
        check(len(hits) == 1, f'{path.name}: ambiguous or missing {comp}:{atom}')
        i = hits[0]
        selected_indices.append(i)
        check(cif['_atom_site.type_symbol'][i].upper() == element, 'Wrong element')
        check(float(cif['_atom_site.occupancy'][i]) == 1.0, 'Partial occupancy')
        check(cif['_atom_site.label_alt_id'][i] in ['.', '?'], 'Alternate conformer')
        coords.append(np.array([float(cif['_atom_site.Cartn_' + a][i]) for a in 'xyz']))
        identifiers.append(f"{cif['_atom_site.auth_asym_id'][i]}:"
                           f"{cif['_atom_site.auth_seq_id'][i]}:{comp}:{atom}")
    o, pg, bridge = coords
    u, v = o - pg, bridge - pg
    distance = float(np.linalg.norm(u))
    angle = math.degrees(math.acos(float(np.clip(np.dot(u, v) /
                      np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))
    check(1.2 < np.linalg.norm(v) < 2.0, 'Implausible P-gamma/bridge bond length')
    # Independent structure parser plus Bio.PDB's angle implementation.
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    check(len(structure) == 1, 'Expected a single coordinate model')
    independent = []
    for comp, name, _ in selections:
        hits = [a for a in structure.get_atoms()
                if a.name == name and a.parent.resname == comp]
        check(len(hits) == 1, 'Independent parser found ambiguous atoms')
        independent.append(hits[0])
    other_d = float(independent[0] - independent[1])
    other_a = math.degrees(float(calc_angle(*(a.get_vector() for a in independent))))
    check(abs(distance - other_d) < 1e-4 and abs(angle - other_a) < 1e-3,
          'Independent measurements disagree')
    row = dict(candidate=path.stem.replace('_generated', ''), file=path.name,
               distance_A=distance, angle_deg=angle,
               near_attack=distance <= 4.0 and angle >= 165.0,
               nucleophile=identifiers[0], gamma_P=identifiers[1],
               leaving_bridge=identifiers[2], topology_verified=True,
               independent_distance_A=other_d, independent_angle_deg=other_a,
               sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    audit = dict(file=path.name, atom_mapping=dict(zip(['MVAP_O16', 'ATP_PG', 'ATP_O3B'], identifiers)),
                 coordinates_A=[p.tolist() for p in coords],
                 terminal_P_bridge_bond_A=float(np.linalg.norm(v)),
                 bonds={comp: {atom: sorted(neighbors) for atom, neighbors in graphs[comp].items()}
                        for comp in ['IMV', 'IAT']})
    return row, audit


def make_plot(rows, output):
    groups = defaultdict(list)
    for row in rows:
        groups[(round(row['distance_A'], 9), round(row['angle_deg'], 9))].append(row)
    ordered = sorted(groups.items(), key=lambda x: (-x[0][1], x[0][0]))
    fig, (ax, table_ax) = plt.subplots(1, 2, figsize=(14.4, 7.8),
                                      gridspec_kw={'width_ratios': [1, 1.15]})
    fig.subplots_adjust(left=.07, right=.98, bottom=.23, top=.83, wspace=.14)
    colors = plt.get_cmap('tab10').colors
    texts, table_rows = [], []
    for i, ((distance, angle), members) in enumerate(ordered):
        group = f'G{i+1}'
        pareto = not any(r['distance_A'] <= members[0]['distance_A'] and
                        r['angle_deg'] >= members[0]['angle_deg'] and
                        (r['distance_A'] < members[0]['distance_A'] or
                         r['angle_deg'] > members[0]['angle_deg']) for r in rows)
        ax.scatter(angle, distance, s=70 + 20*len(members), color=colors[i],
                   edgecolor='#172b4d' if pareto else 'white',
                   linewidth=2 if pareto else 1, zorder=3)
        texts.append(ax.text(angle, distance, group, fontsize=12, weight='bold'))
        ids = [r['candidate'].split('_')[1] for r in members]
        label = '\n'.join(', '.join(ids[j:j+3]) for j in range(0, len(ids), 3))
        table_rows.append([group, label, f'{distance:.3f}', f'{angle:.2f}'])
        for r in members:
            r['geometry_group'] = group
            r['pareto_geometry'] = pareto
    ax.set(xlim=(169.5, 179.3), ylim=(3.53, 3.744),
           xlabel='Attack angle: MVAP O16–ATP PG–O3B (°)',
           ylabel='Attack distance: MVAP O16–ATP PG (Å)')
    ax.set_title('One point per distinct geometry', fontsize=13, pad=16)
    ax.grid(alpha=.18)
    ax.spines[['top', 'right']].set_visible(False)
    adjust_text(texts, ax=ax, expand=(1.8, 1.8),
                arrowprops=dict(arrowstyle='-', color='#777777', lw=.7))
    ax.text(.025, .97, 'Shorter distance ↓\nStraighter attack →', transform=ax.transAxes,
            va='top', fontsize=10, color='#44546a')
    table_ax.axis('off')
    table_ax.set_title('Candidate IDs sharing each point', fontsize=13, pad=16)
    table = table_ax.table(cellText=table_rows,
        colLabels=['Group', 'Candidate suffixes', 'Distance (Å)', 'Angle (°)'],
        colWidths=[.12, .49, .21, .18], cellLoc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    for (i, j), cell in table.get_celld().items():
        cell.set_edgecolor('white')
        if i == 0:
            cell.set_facecolor('#18334d'); cell.set_text_props(color='white', weight='bold')
        else:
            cell.set_facecolor('#f0f4f7' if i % 2 else '#fafbfc')
            if j == 0:
                cell.set_text_props(color=colors[i-1], weight='bold')
    fig.suptitle('RFD4 designs: attack distance vs angle', fontsize=21, weight='bold', y=.965)
    fig.text(.5, .89, f'{len(rows)} candidates · {len(groups)} distinct geometries · '
             f"{sum(r['near_attack'] for r in rows)}/{len(rows)} pass distance ≤4.0 Å and angle ≥165°",
             ha='center', fontsize=12, color='#34506a')
    fig.text(.07, .09, 'Point area reflects the number of overlapping candidates; no coordinate jitter.\n'
             'Dark outlines mark geometry trade-offs: no other group is both shorter and straighter.',
             fontsize=10, color='#44546a', linespacing=1.6)
    fig.text(.07, .025, 'RFD4 atom mapping, verified from bonds in every file: '
             'IMV:O16 → IAT:P41 (PG); angle reference IAT:O27 (O3B).\n'
             'Measurements describe the supplied generated structures; they do not independently validate catalytic activity.',
             fontsize=9, color='#55616c', linespacing=1.5)
    for ext in ['png', 'svg']:
        fig.savefig(output / f'attack_distance_vs_angle.{ext}', dpi=200, facecolor='white')
    svg = output / 'attack_distance_vs_angle.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    output = args.output_dir or args.input_dir / 'analysis'
    paths = sorted(args.input_dir.glob('candidate_*_generated.cif'))
    check(bool(paths), f'No input CIFs in {args.input_dir}')
    output.mkdir(parents=True, exist_ok=True)
    rows, audits = [], []
    for path in paths:
        row, audit = measure(path)
        rows.append(row); audits.append(audit)
    make_plot(rows, output)
    with (output / 'attack_geometry.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    (output / 'atom_identity_audit.json').write_text(json.dumps(audits, indent=2) + '\n')
    print(f'Verified and plotted {len(rows)} candidates; {sum(r["near_attack"] for r in rows)} pass NAC.')
    for row in rows:
        if row['pareto_geometry']:
            print(f"{row['candidate']}: {row['distance_A']:.3f} Å, {row['angle_deg']:.2f}°")


if __name__ == '__main__':
    main()
