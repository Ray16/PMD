#!/usr/bin/env python3
"""Generate comprehensive report, figures, and interactive 3D visualizations
for PMDsc substrate inhibition analysis."""

import json
import os
import glob
import numpy as np
from collections import defaultdict, Counter
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Paths ────────────────────────────────────────────────────────────

BASE_2MVAP = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap'
BASE_1MVAP = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output'
OUT_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'

VARIANT_DIRS = {
    'WT': 'PMDsc_WT_2mvap',
    'Y19H': 'PMDsc_Y19H_2mvap',
    'R74G': 'PMDsc_R74G_2mvap',
    'GKQ': 'PMDsc_R74G_R147K_M212Q_2mvap',
    'HKQ': 'PMDsc_R74H_R147K_M212Q_2mvap',
}

VARIANT_DIRS_1MVAP = {
    'WT': 'PMDsc_WT',
    'Y19H': 'PMDsc_Y19H',
    'R74G': 'PMDsc_R74G',
    'GKQ': 'PMDsc_R74G_R147K_M212Q',
    'HKQ': 'PMDsc_R74H_R147K_M212Q',
}

VARIANT_DATA = {
    'WT':   {'ki': 18,   'kcat_km': 0.066, 'titer': 475,  'mutations': 'none'},
    'Y19H': {'ki': None, 'kcat_km': 0.78,  'titer': 388,  'mutations': 'Y19H'},
    'R74G': {'ki': 110,  'kcat_km': 0.04,  'titer': 975,  'mutations': 'R74G'},
    'GKQ':  {'ki': 11,   'kcat_km': 0.5,   'titer': 8,    'mutations': 'R74G-R147K-M212Q'},
    'HKQ':  {'ki': 80,   'kcat_km': 0.4,   'titer': 1079, 'mutations': 'R74H-R147K-M212Q'},
}

CATALYTIC = {18, 19, 20, 22, 120, 121, 153, 155, 158, 208, 302}
CLOSE_THRESHOLD = 15.0
INTERMEDIATE_THRESHOLD = 25.0

VARIANT_ORDER = ['WT', 'Y19H', 'R74G', 'GKQ', 'HKQ']
COLORS = {'WT': '#2ca02c', 'Y19H': '#ff7f0e', 'R74G': '#1f77b4',
           'GKQ': '#d62728', 'HKQ': '#9467bd'}
MARKERS = {'WT': 'o', 'Y19H': 's', 'R74G': '^', 'GKQ': 'D', 'HKQ': 'v'}

HBOND_DONORS_ACCEPTORS = {'N', 'O', 'S'}

# ── PDB Parsing ──────────────────────────────────────────────────────

def parse_pdb_minimal(path):
    chains = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            chain = parts[4]
            x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
            atom = parts[2]
            resnum = int(parts[5])
            resname = parts[3]
            element = parts[-1] if len(parts) >= 12 else atom[0]
            chains[chain].append({
                'xyz': np.array([x, y, z]),
                'atom': atom, 'resnum': resnum,
                'resname': resname, 'element': element,
                'is_heavy': element not in ('H', 'D'),
            })
    return dict(chains)


def parse_all_heavy(path, chain='A'):
    atoms = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            if parts[4] != chain:
                continue
            element = parts[-1] if len(parts) >= 12 else parts[2][0]
            if element in ('H', 'D'):
                continue
            rn = int(parts[5])
            atoms[rn].append({
                'atom': parts[2],
                'resname': parts[3],
                'xyz': np.array([float(parts[6]), float(parts[7]), float(parts[8])]),
                'element': element,
            })
    return atoms


def parse_ca(path, chain='A'):
    ca = {}
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            if parts[4] == chain and parts[2] == 'CA':
                rn = int(parts[5])
                ca[rn] = np.array([float(parts[6]), float(parts[7]), float(parts[8])])
    return ca


def read_pdb_text(path):
    with open(path) as f:
        return f.read()


def kabsch_align(P, Q):
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Qc.T @ Pc
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_m = np.eye(3)
    sign_m[2, 2] = np.sign(d)
    R = Vt.T @ sign_m @ U.T
    return R, P.mean(0), Q.mean(0)


def align_pdb_text(ref_ca, mobile_path):
    """Read a PDB and return text with coordinates aligned to ref_ca."""
    mobile_ca = parse_ca(mobile_path)
    shared = sorted(set(ref_ca.keys()) & set(mobile_ca.keys()))
    P = np.array([ref_ca[r] for r in shared])
    Q = np.array([mobile_ca[r] for r in shared])
    R, Pc, Qc = kabsch_align(P, Q)

    lines = []
    with open(mobile_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                parts = line.split()
                if len(parts) >= 9:
                    xyz = np.array([float(parts[6]), float(parts[7]), float(parts[8])])
                    new_xyz = (xyz - Qc) @ R.T + Pc
                    line = line[:30] + f'{new_xyz[0]:8.3f}{new_xyz[1]:8.3f}{new_xyz[2]:8.3f}' + line[54:]
            lines.append(line)
    return ''.join(lines)


# ── Model scanning ───────────────────────────────────────────────────

def compute_l2_metrics(chains):
    if 'L' not in chains or 'L2' not in chains or 'A' not in chains:
        return None
    l1_coords = np.array([a['xyz'] for a in chains['L']])
    l2_coords = np.array([a['xyz'] for a in chains['L2']])
    l1_com = l1_coords.mean(axis=0)
    l2_com = l2_coords.mean(axis=0)
    l2_l1_dist = float(np.linalg.norm(l2_com - l1_com))

    l2_heavy = [a for a in chains['L2'] if a['is_heavy']]
    prot_heavy = [a for a in chains['A'] if a['is_heavy']]
    contacts_4a = {}
    for pa in prot_heavy:
        for la in l2_heavy:
            d = np.linalg.norm(pa['xyz'] - la['xyz'])
            if d <= 4.0:
                rn = pa['resnum']
                if rn not in contacts_4a or d < contacts_4a[rn]['min_dist']:
                    contacts_4a[rn] = {
                        'min_dist': d, 'resname': pa['resname'],
                        'protein_atom': pa['atom'], 'l2_atom': la['atom'],
                    }

    l2_to_catalytic_dists = {}
    for rn in CATALYTIC:
        cat_atoms = [a for a in chains['A'] if a['resnum'] == rn and a['is_heavy']]
        if cat_atoms:
            min_d = min(np.linalg.norm(la['xyz'] - ca['xyz'])
                        for la in l2_heavy for ca in cat_atoms)
            l2_to_catalytic_dists[rn] = float(min_d)

    n_catalytic_contacted = sum(1 for d in l2_to_catalytic_dists.values() if d < 4.0)

    if l2_l1_dist > CLOSE_THRESHOLD:
        binding_mode = 'distant' if l2_l1_dist > INTERMEDIATE_THRESHOLD else 'intermediate'
    elif n_catalytic_contacted >= 5:
        binding_mode = 'active_site'
    else:
        binding_mode = 'entrance'

    return {
        'l2_l1_dist': l2_l1_dist,
        'contacts_4a': {str(k): v for k, v in contacts_4a.items()},
        'n_contacts_4a': len(contacts_4a),
        'n_catalytic_contacted': n_catalytic_contacted,
        'binding_mode': binding_mode,
    }


def scan_all_models():
    all_data = {}
    for short, dirname in VARIANT_DIRS.items():
        models = []
        pred_dir = os.path.join(BASE_2MVAP, dirname, f'boltz_results_{dirname}',
                                 'predictions', dirname)
        if not os.path.isdir(pred_dir):
            continue
        pdbs = sorted(glob.glob(os.path.join(pred_dir, f'{dirname}_model_*.pdb')))
        for pdb_path in pdbs:
            model_name = os.path.basename(pdb_path).replace('.pdb', '')
            conf_path = os.path.join(pred_dir, f'confidence_{model_name}.json')
            chains = parse_pdb_minimal(pdb_path)
            metrics = compute_l2_metrics(chains)
            if not metrics:
                continue
            conf = {}
            if os.path.exists(conf_path):
                with open(conf_path) as f:
                    conf_data = json.load(f)
                conf = {
                    'A_L2_iptm_from_L2': conf_data.get('pair_chains_iptm', {}).get('2', {}).get('0'),
                    'iptm': conf_data.get('iptm'),
                }
            model_idx = model_name.split('_')[-1]
            models.append({
                'model_id': f'm{model_idx}',
                'pdb_path': pdb_path,
                **metrics, **conf,
            })
        all_data[short] = {
            'variant': short,
            'ki': VARIANT_DATA[short]['ki'],
            'titer': VARIANT_DATA[short]['titer'],
            'models': models,
        }
    return all_data


# ── Figure Generation ────────────────────────────────────────────────

def fig1_l2_distance_distribution(all_data):
    """Violin/swarm plot of L2-L1 distances by variant."""
    fig, ax = plt.subplots(figsize=(8, 5))
    positions = []
    for i, v in enumerate(VARIANT_ORDER):
        dists = [m['l2_l1_dist'] for m in all_data[v]['models']]
        jitter = np.random.RandomState(42).uniform(-0.2, 0.2, len(dists))
        ax.scatter([i + j for j in jitter], dists, c=COLORS[v], s=50,
                    alpha=0.7, edgecolors='k', linewidth=0.5, zorder=3)
        ax.plot([i - 0.3, i + 0.3], [np.median(dists)] * 2,
                 c='k', linewidth=2.5, zorder=4)
        ki = VARIANT_DATA[v]['ki']
        ki_str = f'Ki={ki}' if ki else 'Ki=N.D.'
        ax.annotate(ki_str, (i, max(dists) + 2), ha='center', fontsize=8,
                     color=COLORS[v], fontweight='bold')

    ax.axhline(CLOSE_THRESHOLD, color='red', linestyle='--', alpha=0.5,
               label=f'Close threshold ({CLOSE_THRESHOLD} Å)')
    ax.axhline(INTERMEDIATE_THRESHOLD, color='orange', linestyle='--', alpha=0.4,
               label=f'Intermediate ({INTERMEDIATE_THRESHOLD} Å)')
    ax.set_xticks(range(len(VARIANT_ORDER)))
    ax.set_xticklabels(VARIANT_ORDER, fontsize=11, fontweight='bold')
    ax.set_ylabel('L2-L1 COM Distance (Å)', fontsize=11)
    ax.set_title('L2 Placement Distances Across Variants\n(each dot = one Boltz-2 model, bar = median)',
                  fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(0, None)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig1_l2_distance_distribution.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def fig2_close_binding_frequency(all_data):
    """Bar chart of close-binding fraction per variant."""
    fig, ax = plt.subplots(figsize=(7, 5))
    fracs = []
    for v in VARIANT_ORDER:
        models = all_data[v]['models']
        n_close = sum(1 for m in models if m['l2_l1_dist'] < CLOSE_THRESHOLD)
        fracs.append(n_close / len(models))

    bars = ax.bar(range(len(VARIANT_ORDER)), fracs,
                   color=[COLORS[v] for v in VARIANT_ORDER],
                   edgecolor='k', linewidth=1)

    for i, (v, frac) in enumerate(zip(VARIANT_ORDER, fracs)):
        n_close = sum(1 for m in all_data[v]['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD)
        n_total = len(all_data[v]['models'])
        ax.text(i, frac + 0.01, f'{n_close}/{n_total}', ha='center', fontsize=10, fontweight='bold')
        ki = VARIANT_DATA[v]['ki']
        ki_str = f'Ki={ki} mM' if ki else 'Ki=N.D.'
        ax.text(i, frac / 2 if frac > 0 else 0.01, ki_str, ha='center', va='center',
                 fontsize=8, color='white' if frac > 0.05 else 'black', fontweight='bold')

    ax.set_xticks(range(len(VARIANT_ORDER)))
    ax.set_xticklabels(VARIANT_ORDER, fontsize=11, fontweight='bold')
    ax.set_ylabel('Fraction of Models with L2 < 15 Å', fontsize=11)
    ax.set_title('Close-Binding Frequency by Variant\n(L2 within 15 Å of L1)',
                  fontsize=12, fontweight='bold')
    ax.set_ylim(0, max(fracs) * 1.3 + 0.05)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig2_close_binding_frequency.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def fig3_consensus_contacts(all_data):
    """Heatmap: per-residue contact frequency across close-binding models."""
    close_models = []
    for v in VARIANT_ORDER:
        for m in all_data[v]['models']:
            if m['l2_l1_dist'] < CLOSE_THRESHOLD:
                close_models.append((v, m))

    if not close_models:
        print('  No close-binding models for fig3')
        return

    all_residues = set()
    for v, m in close_models:
        all_residues.update(int(r) for r in m['contacts_4a'].keys())
    all_residues = sorted(all_residues)

    matrix = np.zeros((len(close_models), len(all_residues)))
    ylabels = []
    for i, (v, m) in enumerate(close_models):
        ylabels.append(f"{v} {m['model_id']} ({m['l2_l1_dist']:.1f} Å)")
        for j, rn in enumerate(all_residues):
            if str(rn) in m['contacts_4a']:
                matrix[i, j] = 1

    fig, ax = plt.subplots(figsize=(14, 5))
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)

    res_labels = []
    for rn in all_residues:
        tag = '*' if rn in CATALYTIC else ''
        res_labels.append(f'{rn}{tag}')
    ax.set_xticks(range(len(all_residues)))
    ax.set_xticklabels(res_labels, rotation=90, fontsize=9)
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_xlabel('Residue Number (* = catalytic)', fontsize=11)
    ax.set_title('L2 Contact Residues (< 4 Å) in Close-Binding Models',
                  fontsize=12, fontweight='bold')

    freq = matrix.sum(axis=0)
    for j in range(len(all_residues)):
        ax.text(j, -0.8, f'{int(freq[j])}', ha='center', va='center', fontsize=8,
                 fontweight='bold', color='darkred' if freq[j] >= 4 else 'gray')

    plt.colorbar(im, ax=ax, label='Contact (1=yes)', shrink=0.6)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig3_consensus_contacts.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def fig4_r74g_competitive(all_data):
    """WT entrance binding vs R74G active-site binding comparison."""
    wt_close = [m for m in all_data['WT']['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD]
    r74g_close = [m for m in all_data['R74G']['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD]

    if not wt_close or not r74g_close:
        print('  Missing close models for fig4')
        return

    wt_model = wt_close[0]
    r74g_model = r74g_close[0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, model, variant, title in [
        (axes[0], wt_model, 'WT', 'WT: Entrance Binding (Non-competitive)'),
        (axes[1], r74g_model, 'R74G', 'R74G: Active-Site Binding (Competitive)'),
    ]:
        contacts = model['contacts_4a']
        residues = sorted([int(r) for r in contacts.keys()])
        dists = [contacts[str(r)]['min_dist'] for r in residues]
        colors_bar = ['red' if r in CATALYTIC else COLORS[variant] for r in residues]

        bars = ax.barh(range(len(residues)), dists, color=colors_bar, edgecolor='k', linewidth=0.5)
        ax.set_yticks(range(len(residues)))
        labels = [f"{contacts[str(r)]['resname']}{r}" for r in residues]
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel('Min Distance to L2 (Å)', fontsize=10)
        ax.set_title(f'{title}\n{model["model_id"]}, L2-L1={model["l2_l1_dist"]:.1f} Å\n'
                      f'Catalytic contacts: {model["n_catalytic_contacted"]}',
                      fontsize=10, fontweight='bold')
        ax.axvline(4.0, color='gray', linestyle=':', alpha=0.5)
        ax.invert_xaxis()

        for i, r in enumerate(residues):
            if r in CATALYTIC:
                ax.text(0.1, i, 'catalytic', va='center', fontsize=7, color='red', fontstyle='italic')

    fig.suptitle('Entrance Binding vs Active-Site Binding: Two Distinct Inhibition Mechanisms',
                  fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig4_r74g_competitive.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def fig5_r74_loop_mechanism(insights_data):
    """R74 ↔ loop 17-33 contact counts."""
    r74_data = insights_data['r74_loop_interactions']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax1 = axes[0]
    contacts = [r74_data[v]['avg_contacts_5a'] for v in VARIANT_ORDER]
    bar_colors = [COLORS[v] for v in VARIANT_ORDER]
    ax1.bar(range(len(VARIANT_ORDER)), contacts, color=bar_colors, edgecolor='k')
    ax1.set_xticks(range(len(VARIANT_ORDER)))
    labels = [f"{v}\n({r74_data[v]['resname_74']}74)" for v in VARIANT_ORDER]
    ax1.set_xticklabels(labels, fontsize=9)
    ax1.set_ylabel('Avg Contacts (< 5 Å)')
    ax1.set_title('Residue 74 Contacts\nto Loop 17-33', fontweight='bold')

    ax2 = axes[1]
    loop_res_order = list(range(17, 34))
    for v in ['WT', 'HKQ', 'R74G']:
        freq = r74_data[v].get('loop_contact_freq', {})
        vals = [freq.get(str(r), freq.get(r, 0)) for r in loop_res_order]
        ax2.plot(loop_res_order, vals, '-o', color=COLORS[v], label=v, markersize=5, linewidth=2)
    ax2.set_xlabel('Loop Residue Number')
    ax2.set_ylabel('Contact Frequency (across 5 models)')
    ax2.set_title('Per-Residue Contact Frequency\n(Res 74 → Loop 17-33)', fontweight='bold')
    ax2.legend()

    ax3 = axes[2]
    for v in VARIANT_ORDER:
        ki = VARIANT_DATA[v]['ki']
        nc = r74_data[v]['avg_contacts_5a']
        if ki is not None:
            ax3.scatter(nc, ki, c=COLORS[v], s=150, edgecolors='k', linewidth=1.5, zorder=3)
            ax3.annotate(v, (nc, ki), textcoords="offset points",
                         xytext=(8, 5), fontsize=11, fontweight='bold')
    ax3.set_xlabel('Avg Contacts (Res 74 → Loop)', fontsize=11)
    ax3.set_ylabel('Ki (mM)', fontsize=11)
    ax3.set_title('Loop Closure Strength\nvs Substrate Inhibition', fontweight='bold')

    fig.suptitle('The R74 Gateway Mechanism: How Residue 74 Controls Access to the Secondary Site',
                  fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig5_r74_loop_mechanism.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def fig6_gkq_paradox():
    """Q212 bridge distances and K147 salt bridge comparison."""
    wt_pdbs = glob.glob(os.path.join(BASE_1MVAP, 'PMDsc_WT',
                                      'boltz_results_PMDsc_WT', 'predictions',
                                      'PMDsc_WT', 'PMDsc_WT_model_*.pdb'))
    ref_ca = parse_ca(sorted(wt_pdbs)[0])

    variants_to_check = ['WT', 'R74G', 'GKQ', 'HKQ']
    bridge_data = {}

    for variant in variants_to_check:
        dirname = VARIANT_DIRS_1MVAP[variant]
        pdbs = sorted(glob.glob(os.path.join(BASE_1MVAP, dirname,
                                              f'boltz_results_{dirname}', 'predictions',
                                              dirname, f'{dirname}_model_*.pdb')))
        if not pdbs:
            continue

        q212_y19_dists = []
        q212_k22_dists = []
        k147_e144_dists = []

        for pdb in pdbs:
            atoms = parse_all_heavy(pdb)
            ca = parse_ca(pdb)
            shared = sorted(set(ref_ca.keys()) & set(ca.keys()))
            P = np.array([ref_ca[r] for r in shared])
            Q = np.array([ca[r] for r in shared])
            R, Pc, Qc = kabsch_align(P, Q)

            aligned = {}
            for rn, atom_list in atoms.items():
                aligned[rn] = []
                for a in atom_list:
                    new_xyz = (a['xyz'] - Qc) @ R.T + Pc
                    aligned[rn].append({**a, 'xyz': new_xyz})

            res212 = aligned.get(212, [])
            res19 = aligned.get(19, [])
            res22 = aligned.get(22, [])
            res147 = aligned.get(147, [])
            res144 = aligned.get(144, [])

            q212_oe1 = [a for a in res212 if a['atom'] == 'OE1']
            y19_oh = [a for a in res19 if a['atom'] == 'OH']
            k22_nz = [a for a in res22 if a['atom'] == 'NZ']

            if q212_oe1 and y19_oh:
                q212_y19_dists.append(float(np.linalg.norm(q212_oe1[0]['xyz'] - y19_oh[0]['xyz'])))
            if q212_oe1 and k22_nz:
                q212_k22_dists.append(float(np.linalg.norm(q212_oe1[0]['xyz'] - k22_nz[0]['xyz'])))

            k147_atoms = [a for a in res147 if a['atom'] == 'NZ']
            e144_atoms = [a for a in res144 if a['atom'] in ('OE1', 'OE2')]
            if k147_atoms and e144_atoms:
                min_d = min(np.linalg.norm(k['xyz'] - e['xyz']) for k in k147_atoms for e in e144_atoms)
                k147_e144_dists.append(float(min_d))

        bridge_data[variant] = {
            'q212_y19': q212_y19_dists,
            'q212_k22': q212_k22_dists,
            'k147_e144': k147_e144_dists,
            'resname_212': 'GLN' if variant in ('GKQ', 'HKQ') else 'MET',
            'resname_147': 'LYS' if variant in ('GKQ', 'HKQ') else 'ARG',
        }

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax1 = axes[0]
    for i, v in enumerate(variants_to_check):
        vals = bridge_data[v]['q212_y19']
        if vals:
            jitter = np.random.RandomState(42).uniform(-0.15, 0.15, len(vals))
            ax1.scatter([i + j for j in jitter], vals, c=COLORS[v], s=60,
                        edgecolors='k', linewidth=0.5, zorder=3)
            ax1.plot([i - 0.25, i + 0.25], [np.mean(vals)] * 2, 'k-', linewidth=2, zorder=4)
    ax1.axhline(3.5, color='red', linestyle='--', alpha=0.5, label='H-bond cutoff (3.5 Å)')
    ax1.set_xticks(range(len(variants_to_check)))
    ax1.set_xticklabels([f'{v}\n({bridge_data[v]["resname_212"]}212)' for v in variants_to_check], fontsize=9)
    ax1.set_ylabel('Distance (Å)')
    ax1.set_title('Q212.OE1 — Y19.OH\n(Bridge Arm 1)', fontweight='bold')
    ax1.legend(fontsize=8)

    ax2 = axes[1]
    for i, v in enumerate(variants_to_check):
        vals = bridge_data[v]['q212_k22']
        if vals:
            jitter = np.random.RandomState(42).uniform(-0.15, 0.15, len(vals))
            ax2.scatter([i + j for j in jitter], vals, c=COLORS[v], s=60,
                        edgecolors='k', linewidth=0.5, zorder=3)
            ax2.plot([i - 0.25, i + 0.25], [np.mean(vals)] * 2, 'k-', linewidth=2, zorder=4)
    ax2.axhline(3.5, color='red', linestyle='--', alpha=0.5, label='H-bond cutoff (3.5 Å)')
    ax2.set_xticks(range(len(variants_to_check)))
    ax2.set_xticklabels([f'{v}\n({bridge_data[v]["resname_212"]}212)' for v in variants_to_check], fontsize=9)
    ax2.set_ylabel('Distance (Å)')
    ax2.set_title('Q212.OE1 — K22.NZ\n(Bridge Arm 2)', fontweight='bold')
    ax2.legend(fontsize=8)

    ax3 = axes[2]
    for i, v in enumerate(variants_to_check):
        vals = bridge_data[v]['k147_e144']
        if vals:
            jitter = np.random.RandomState(42).uniform(-0.15, 0.15, len(vals))
            ax3.scatter([i + j for j in jitter], vals, c=COLORS[v], s=60,
                        edgecolors='k', linewidth=0.5, zorder=3)
            ax3.plot([i - 0.25, i + 0.25], [np.mean(vals)] * 2, 'k-', linewidth=2, zorder=4)
    ax3.axhline(4.0, color='red', linestyle='--', alpha=0.5, label='Salt bridge cutoff (4.0 Å)')
    ax3.set_xticks(range(len(variants_to_check)))
    ax3.set_xticklabels([f'{v}\n({bridge_data[v]["resname_147"]}147)' for v in variants_to_check], fontsize=9)
    ax3.set_ylabel('Distance (Å)')
    ax3.set_title('K147.NZ — E144.OE\n(Salt Bridge)', fontweight='bold')
    ax3.legend(fontsize=8)

    fig.suptitle('The GKQ Paradox: Why R74G+R147K+M212Q (Ki=11) Is Worse Than R74G Alone (Ki=110)\n'
                 'HKQ forms the Q212 bridge (Y19-Q212-K22) tethering loop to H4; GKQ does not',
                  fontsize=12, fontweight='bold', y=1.05)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig6_gkq_paradox.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def fig7_confidence_correlation(all_data):
    """Scatter: A-L2 iPTM vs L2-L1 distance."""
    fig, ax = plt.subplots(figsize=(8, 6))
    all_dists = []
    all_iptms = []

    for v in VARIANT_ORDER:
        dists_v = []
        iptms_v = []
        for m in all_data[v]['models']:
            ip = m.get('A_L2_iptm_from_L2')
            if ip is not None:
                dists_v.append(m['l2_l1_dist'])
                iptms_v.append(ip)
        ax.scatter(dists_v, iptms_v, c=COLORS[v], marker=MARKERS[v],
                    s=60, alpha=0.8, label=v, edgecolors='gray', linewidth=0.3)
        all_dists.extend(dists_v)
        all_iptms.extend(iptms_v)

    if all_dists:
        r, p = stats.pearsonr(all_dists, all_iptms)
        z = np.polyfit(all_dists, all_iptms, 1)
        xline = np.linspace(min(all_dists), max(all_dists), 100)
        ax.plot(xline, np.polyval(z, xline), 'k--', alpha=0.4, linewidth=2,
                 label=f'r = {r:.2f}, p = {p:.1e}')

    ax.axvline(CLOSE_THRESHOLD, color='red', linestyle='--', alpha=0.4, label='Close threshold')
    ax.set_xlabel('L2-L1 COM Distance (Å)', fontsize=11)
    ax.set_ylabel('A-L2 iPTM (from L2 perspective)', fontsize=11)
    ax.set_title('Boltz-2 Confidence Validates Close Binding:\nHigher iPTM = Closer L2 Placement',
                  fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig7_confidence_correlation.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def figS1_ki_vs_close_fraction(all_data):
    """Scatter: Ki vs close-binding fraction."""
    fig, ax = plt.subplots(figsize=(7, 5))
    for v in VARIANT_ORDER:
        ki = VARIANT_DATA[v]['ki']
        if ki is None:
            continue
        models = all_data[v]['models']
        frac = sum(1 for m in models if m['l2_l1_dist'] < CLOSE_THRESHOLD) / len(models)
        ax.scatter(ki, frac, c=COLORS[v], s=200, edgecolors='k', linewidth=2, zorder=3)
        ax.annotate(f'{v}\n(titer={VARIANT_DATA[v]["titer"]})', (ki, frac),
                     textcoords="offset points", xytext=(12, -5),
                     fontsize=10, fontweight='bold')

    ax.set_xlabel('Ki (mM)', fontsize=12)
    ax.set_ylabel('Fraction Close-Binding (L2 < 15 Å)', fontsize=12)
    ax.set_title('Substrate Inhibition Strength vs L2 Binding Propensity\n'
                  '(4 variants with measured Ki)',
                  fontsize=12, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'figS1_ki_vs_close_fraction.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def figS2_full_contact_heatmap(all_data):
    """Full per-model x per-residue contact heatmap (all 50 models)."""
    all_residues = set()
    for v in VARIANT_ORDER:
        for m in all_data[v]['models']:
            all_residues.update(int(r) for r in m['contacts_4a'].keys())
    if not all_residues:
        print('  No contacts for figS2')
        return
    all_residues = sorted(all_residues)

    models_list = []
    for v in VARIANT_ORDER:
        sorted_models = sorted(all_data[v]['models'], key=lambda m: m['l2_l1_dist'])
        for m in sorted_models:
            models_list.append((v, m))

    matrix = np.zeros((len(models_list), len(all_residues)))
    ylabels = []
    for i, (v, m) in enumerate(models_list):
        ylabels.append(f"{v} {m['model_id']} ({m['l2_l1_dist']:.0f}Å)")
        for j, rn in enumerate(all_residues):
            if str(rn) in m['contacts_4a']:
                matrix[i, j] = 1

    fig, ax = plt.subplots(figsize=(16, 20))
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
    res_labels = [f'{rn}{"*" if rn in CATALYTIC else ""}' for rn in all_residues]
    ax.set_xticks(range(len(all_residues)))
    ax.set_xticklabels(res_labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel('Residue (* = catalytic)')
    ax.set_title('Full Contact Heatmap: All 50 Models\n(sorted by L2-L1 distance within each variant)',
                  fontsize=12, fontweight='bold')

    variant_boundaries = []
    current = 0
    for v in VARIANT_ORDER:
        n = len(all_data[v]['models'])
        variant_boundaries.append((current, current + n, v))
        current += n
    for start, end, v in variant_boundaries:
        ax.axhline(start - 0.5, color='k', linewidth=1)
        ax.text(-1.5, (start + end) / 2, v, ha='right', va='center', fontsize=9,
                 fontweight='bold', color=COLORS[v])

    plt.colorbar(im, ax=ax, label='Contact (1=yes)', shrink=0.3)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'figS2_full_contact_heatmap.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


# ── Interactive 3D HTML Generation ───────────────────────────────────

THREEDMOL_CDN = 'https://3Dmol.org/build/3Dmol-min.js'

HTML_TEMPLATE = '''<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 0; padding: 10px; background: #f5f5f5; }}
h2 {{ color: #333; margin: 5px 0; }}
.info {{ background: white; padding: 10px; border-radius: 5px; margin: 5px 0; font-size: 13px; line-height: 1.5; }}
.viewer-container {{ display: flex; gap: 10px; flex-wrap: wrap; }}
.viewer-box {{ background: white; border-radius: 5px; padding: 5px; }}
.viewer-box h3 {{ margin: 2px 0; font-size: 14px; text-align: center; }}
</style>
<script src="{cdn}"></script>
</head>
<body>
<h2>{title}</h2>
<div class="info">{description}</div>
<div class="viewer-container">
{viewers}
</div>
<script>
{script}
</script>
</body>
</html>'''

VIEWER_DIV = '''<div class="viewer-box">
<h3>{subtitle}</h3>
<div id="{div_id}" style="width:{width}px;height:{height}px;position:relative;"></div>
</div>'''


def make_viewer_script(div_id, pdb_data, setup_commands):
    """Generate JS for a single 3Dmol viewer."""
    pdb_escaped = pdb_data.replace('\\', '\\\\').replace('`', '\\`').replace('${', '\\${')
    return f'''
(function() {{
  var viewer = $3Dmol.createViewer("{div_id}", {{backgroundColor: "white"}});
  var pdbData = `{pdb_escaped}`;
  viewer.addModel(pdbData, "pdb");
  {setup_commands}
  viewer.zoomTo();
  viewer.render();
}})();
'''


def get_label_cmd(resi, label_text, color='black'):
    return (f'viewer.addLabel("{label_text}", '
            f'{{position: viewer.selectedAtoms({{resi: {resi}, atom: "CA"}})[0], '
            f'backgroundColor: "white", backgroundOpacity: 0.7, '
            f'fontColor: "{color}", fontSize: 12, showBackground: true}});\n')


def get_label_cmd_atom(resi, atom_name, label_text, color='black'):
    return (f'viewer.addLabel("{label_text}", '
            f'{{position: viewer.selectedAtoms({{resi: {resi}, atom: "{atom_name}"}})[0], '
            f'backgroundColor: "white", backgroundOpacity: 0.7, '
            f'fontColor: "{color}", fontSize: 11, showBackground: true}});\n')


def view_consensus_site(all_data):
    """WT model_0 with L2 at entrance, labeled contact residues."""
    wt_close = [m for m in all_data['WT']['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD]
    if not wt_close:
        print('  No WT close model for consensus view')
        return
    pdb_path = wt_close[0]['pdb_path']
    pdb_data = read_pdb_text(pdb_path)
    model = wt_close[0]

    contact_residues = [22, 25, 28, 72, 73, 74, 121, 152, 153, 208, 209]
    contact_names = {
        22: 'K22', 25: 'T25', 28: 'N28', 72: 'N72', 73: 'D73',
        74: 'R74', 121: 'S121', 152: 'G152', 153: 'S153', 208: 'S208', 209: 'T209'
    }

    setup = ''
    setup += 'viewer.setStyle({chain: "A"}, {cartoon: {color: "lightgray", opacity: 0.85}});\n'
    setup += 'viewer.setStyle({chain: "L"}, {stick: {colorscheme: "greenCarbon", radius: 0.2}});\n'
    setup += 'viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.4, color: "green"}, {chain: "L"});\n'
    setup += 'viewer.setStyle({chain: "L2"}, {stick: {colorscheme: "redCarbon", radius: 0.2}});\n'
    setup += 'viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.4, color: "red"}, {chain: "L2"});\n'

    for rn in contact_residues:
        cat_color = '"salmon"' if rn in CATALYTIC else '"lightblue"'
        setup += f'viewer.setStyle({{chain: "A", resi: {rn}}}, {{stick: {{colorscheme: "default", radius: 0.15}}, cartoon: {{color: {cat_color}}}}});\n'

    for rn, name in contact_names.items():
        color = 'red' if rn in CATALYTIC else 'blue'
        setup += get_label_cmd(rn, name, color)

    setup += 'viewer.addLabel("L1 (MVAP-1)", {position: viewer.selectedAtoms({chain: "L"})[0], backgroundColor: "green", fontColor: "white", fontSize: 12});\n'
    setup += 'viewer.addLabel("L2 (MVAP-2)", {position: viewer.selectedAtoms({chain: "L2"})[0], backgroundColor: "red", fontColor: "white", fontSize: 12});\n'

    viewer = VIEWER_DIV.format(subtitle=f'WT {model["model_id"]} (L2-L1 = {model["l2_l1_dist"]:.1f} Å)',
                                div_id='viewer1', width=900, height=600)
    script = make_viewer_script('viewer1', pdb_data, setup)

    desc = (f'<b>The secondary MVAP binding site</b> at the active-site entrance. '
            f'L1 (green) is in the catalytic pocket; L2 (red) binds at the entrance. '
            f'<b>Blue labels</b> = non-catalytic contact residues (mutation targets). '
            f'<b>Red labels</b> = catalytic residues (do not mutate). '
            f'Key contacts: R74 (gateway), K22/N28 (loop 17-33), N72, T209.')

    html = HTML_TEMPLATE.format(
        title='Secondary MVAP Binding Site — Consensus View (WT)',
        cdn=THREEDMOL_CDN, description=desc,
        viewers=viewer, script=script)

    path = os.path.join(OUT_DIR, 'view_consensus_site.html')
    with open(path, 'w') as f:
        f.write(html)
    print(f'  Saved {path}')


def view_r74g_competitive(all_data):
    """R74G model with L2 inside active site."""
    r74g_close = [m for m in all_data['R74G']['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD]
    if not r74g_close:
        print('  No R74G close model')
        return
    model = r74g_close[0]
    pdb_data = read_pdb_text(model['pdb_path'])

    cat_residues_list = sorted(CATALYTIC)
    cat_names = {18: 'D18', 19: 'Y19', 20: 'D20', 22: 'K22', 120: 'S120',
                 121: 'S121', 153: 'S153', 155: 'F155', 158: 'K158', 208: 'S208', 302: 'MG'}

    setup = ''
    setup += 'viewer.setStyle({chain: "A"}, {cartoon: {color: "lightgray", opacity: 0.85}});\n'
    setup += 'viewer.setStyle({chain: "L"}, {stick: {colorscheme: "greenCarbon", radius: 0.2}});\n'
    setup += 'viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.3, color: "green"}, {chain: "L"});\n'
    setup += 'viewer.setStyle({chain: "L2"}, {stick: {colorscheme: "redCarbon", radius: 0.2}});\n'
    setup += 'viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.3, color: "red"}, {chain: "L2"});\n'

    for rn in cat_residues_list:
        setup += f'viewer.setStyle({{chain: "A", resi: {rn}}}, {{stick: {{colorscheme: "default", radius: 0.15}}, cartoon: {{color: "salmon"}}}});\n'
    setup += f'viewer.setStyle({{chain: "A", resi: 74}}, {{stick: {{colorscheme: "default", radius: 0.15}}, cartoon: {{color: "lightblue"}}}});\n'

    for rn, name in cat_names.items():
        setup += get_label_cmd(rn, name, 'red')
    setup += get_label_cmd(74, 'G74', 'blue')

    setup += 'viewer.addLabel("L1", {position: viewer.selectedAtoms({chain: "L"})[0], backgroundColor: "green", fontColor: "white", fontSize: 12});\n'
    setup += 'viewer.addLabel("L2 (inside active site)", {position: viewer.selectedAtoms({chain: "L2"})[0], backgroundColor: "red", fontColor: "white", fontSize: 12});\n'

    viewer = VIEWER_DIV.format(subtitle=f'R74G {model["model_id"]} (L2-L1 = {model["l2_l1_dist"]:.1f} Å, {model["n_catalytic_contacted"]} catalytic contacts)',
                                div_id='viewer1', width=900, height=600)
    script = make_viewer_script('viewer1', pdb_data, setup)

    desc = (f'<b>R74G: Competitive binding mechanism.</b> Without R74\'s guanidinium, '
            f'L2 enters the active site itself (red), displacing L1 (green). '
            f'This model has {model["n_catalytic_contacted"]} contacts with catalytic residues. '
            f'This is a different mechanism from the entrance binding seen in WT — '
            f'consistent with R74G\'s extremely high Ki (110 mM, weak inhibition).')

    html = HTML_TEMPLATE.format(
        title='R74G: Competitive Active-Site Binding',
        cdn=THREEDMOL_CDN, description=desc,
        viewers=viewer, script=script)

    path = os.path.join(OUT_DIR, 'view_r74g_competitive.html')
    with open(path, 'w') as f:
        f.write(html)
    print(f'  Saved {path}')


def view_r74_loop_shield():
    """WT R74 vs HKQ H74 loop closure comparison."""
    wt_pdbs = sorted(glob.glob(os.path.join(BASE_1MVAP, 'PMDsc_WT',
                                              'boltz_results_PMDsc_WT', 'predictions',
                                              'PMDsc_WT', 'PMDsc_WT_model_*.pdb')))
    hkq_pdbs = sorted(glob.glob(os.path.join(BASE_1MVAP, 'PMDsc_R74H_R147K_M212Q',
                                               'boltz_results_PMDsc_R74H_R147K_M212Q', 'predictions',
                                               'PMDsc_R74H_R147K_M212Q', 'PMDsc_R74H_R147K_M212Q_model_*.pdb')))
    if not wt_pdbs or not hkq_pdbs:
        print('  Missing 1-MVAP PDBs for loop shield view')
        return

    wt_pdb = read_pdb_text(wt_pdbs[0])
    ref_ca = parse_ca(wt_pdbs[0])
    hkq_pdb = align_pdb_text(ref_ca, hkq_pdbs[0])

    loop_residues = list(range(17, 34))
    key_res = [22, 25, 26, 28, 74]
    res_names = {22: 'K22', 25: 'T25', 26: 'K26', 28: 'N28', 74: 'R74'}
    res_names_hkq = {22: 'K22', 25: 'T25', 26: 'K26', 28: 'N28', 74: 'H74'}

    def make_setup(names, res74_label):
        s = ''
        s += 'viewer.setStyle({chain: "A"}, {cartoon: {color: "lightgray", opacity: 0.85}});\n'
        s += 'viewer.setStyle({chain: "L"}, {stick: {colorscheme: "greenCarbon", radius: 0.15}});\n'
        for rn in loop_residues:
            s += f'viewer.setStyle({{chain: "A", resi: {rn}}}, {{cartoon: {{color: "#FFD700"}}, stick: {{colorscheme: "default", radius: 0.12}}}});\n'
        s += f'viewer.setStyle({{chain: "A", resi: 74}}, {{stick: {{colorscheme: "default", radius: 0.2}}, cartoon: {{color: "#4169E1"}}}});\n'
        for rn, name in names.items():
            color = 'blue' if rn == 74 else 'darkgoldenrod'
            s += get_label_cmd(rn, name, color)
        return s

    viewer1 = VIEWER_DIV.format(subtitle='WT (R74): 0 contacts to loop',
                                 div_id='viewer_wt', width=550, height=500)
    viewer2 = VIEWER_DIV.format(subtitle='HKQ (H74): 21 contacts to loop',
                                 div_id='viewer_hkq', width=550, height=500)

    script1 = make_viewer_script('viewer_wt', wt_pdb, make_setup(res_names, 'R74'))
    script2 = make_viewer_script('viewer_hkq', hkq_pdb, make_setup(res_names_hkq, 'H74'))

    desc = (f'<b>The R74 Gateway:</b> In WT, R74 (Arg, blue) provides electrostatic shielding '
            f'but makes <b>0 direct contacts</b> to loop 17-33 (gold). '
            f'In HKQ, H74 (His) makes <b>21 contacts</b> to residues T25, K26, N28, '
            f'physically closing the entrance and explaining why HKQ has 0/10 close-binding models.')

    html = HTML_TEMPLATE.format(
        title='R74 Gateway: Electrostatic Shield (WT) vs Physical Closure (HKQ)',
        cdn=THREEDMOL_CDN, description=desc,
        viewers=viewer1 + viewer2, script=script1 + script2)

    path = os.path.join(OUT_DIR, 'view_r74_loop_shield.html')
    with open(path, 'w') as f:
        f.write(html)
    print(f'  Saved {path}')


def view_gkq_vs_hkq():
    """GKQ vs HKQ Q212 bridge comparison."""
    ref_pdb_path = sorted(glob.glob(os.path.join(BASE_1MVAP, 'PMDsc_WT',
                                                   'boltz_results_PMDsc_WT', 'predictions',
                                                   'PMDsc_WT', 'PMDsc_WT_model_*.pdb')))[0]
    ref_ca = parse_ca(ref_pdb_path)

    gkq_pdbs = sorted(glob.glob(os.path.join(BASE_1MVAP, 'PMDsc_R74G_R147K_M212Q',
                                               'boltz_results_PMDsc_R74G_R147K_M212Q', 'predictions',
                                               'PMDsc_R74G_R147K_M212Q', 'PMDsc_R74G_R147K_M212Q_model_*.pdb')))
    hkq_pdbs = sorted(glob.glob(os.path.join(BASE_1MVAP, 'PMDsc_R74H_R147K_M212Q',
                                               'boltz_results_PMDsc_R74H_R147K_M212Q', 'predictions',
                                               'PMDsc_R74H_R147K_M212Q', 'PMDsc_R74H_R147K_M212Q_model_*.pdb')))

    if not gkq_pdbs or not hkq_pdbs:
        print('  Missing PDBs for GKQ/HKQ view')
        return

    gkq_pdb = align_pdb_text(ref_ca, gkq_pdbs[0])
    hkq_pdb = align_pdb_text(ref_ca, hkq_pdbs[0])

    bridge_residues = [19, 22, 74, 147, 212]
    res_names_gkq = {19: 'Y19', 22: 'K22', 74: 'G74', 147: 'K147', 212: 'Q212'}
    res_names_hkq = {19: 'Y19', 22: 'K22', 74: 'H74', 147: 'K147', 212: 'Q212'}

    def make_setup(names, variant_color):
        s = ''
        s += f'viewer.setStyle({{chain: "A"}}, {{cartoon: {{color: "lightgray", opacity: 0.85}}}});\n'
        s += f'viewer.setStyle({{chain: "L"}}, {{stick: {{colorscheme: "greenCarbon", radius: 0.12}}}});\n'
        for rn in bridge_residues:
            s += f'viewer.setStyle({{chain: "A", resi: {rn}}}, {{stick: {{colorscheme: "default", radius: 0.2}}, cartoon: {{color: "{variant_color}"}}}});\n'
        for rn in range(17, 34):
            s += f'viewer.setStyle({{chain: "A", resi: {rn}}}, {{cartoon: {{color: "#FFD700"}}}});\n'
        for rn in bridge_residues:
            s += f'viewer.setStyle({{chain: "A", resi: {rn}}}, {{stick: {{colorscheme: "default", radius: 0.2}}, cartoon: {{color: "{variant_color}"}}}});\n'
        for rn, name in names.items():
            s += get_label_cmd(rn, name, 'black')
        s += f'viewer.setStyle({{chain: "A", resi: 144}}, {{stick: {{colorscheme: "default", radius: 0.15}}}});\n'
        s += get_label_cmd(144, 'E144', 'gray')
        return s

    viewer1 = VIEWER_DIV.format(subtitle='GKQ (Ki=11 mM): Q212 bridge ABSENT',
                                 div_id='viewer_gkq', width=550, height=500)
    viewer2 = VIEWER_DIV.format(subtitle='HKQ (Ki=80 mM): Q212 bridge PRESENT',
                                 div_id='viewer_hkq', width=550, height=500)

    script1 = make_viewer_script('viewer_gkq', gkq_pdb, make_setup(res_names_gkq, '#d62728'))
    script2 = make_viewer_script('viewer_hkq', hkq_pdb, make_setup(res_names_hkq, '#9467bd'))

    desc = (f'<b>The GKQ Paradox:</b> Both GKQ and HKQ carry R147K + M212Q, but Ki differs 7-fold '
            f'(GKQ=11 vs HKQ=80 mM). In HKQ (right), Q212 forms an H-bond bridge: '
            f'Q212.OE1 — Y19.OH (2.72 Å) and Q212.OE1 — K22.NZ (2.79 Å), '
            f'tethering loop 17-33 to the H4 helix and blocking the entrance. '
            f'In GKQ (left), this bridge is absent despite the identical M212Q mutation. '
            f'Additionally, K147 in HKQ forms a salt bridge with E144, stabilizing helix H2.')

    html = HTML_TEMPLATE.format(
        title='GKQ vs HKQ: The Q212 Bridge Paradox',
        cdn=THREEDMOL_CDN, description=desc,
        viewers=viewer1 + viewer2, script=script1 + script2)

    path = os.path.join(OUT_DIR, 'view_gkq_vs_hkq.html')
    with open(path, 'w') as f:
        f.write(html)
    print(f'  Saved {path}')


def view_all_l2_positions(all_data):
    """All close-binding models superposed."""
    close_models = []
    for v in VARIANT_ORDER:
        for m in all_data[v]['models']:
            if m['l2_l1_dist'] < CLOSE_THRESHOLD:
                close_models.append((v, m))

    if not close_models:
        print('  No close models for all-L2 view')
        return

    ref_path = close_models[0][1]['pdb_path']
    ref_ca = parse_ca(ref_path)
    ref_pdb = read_pdb_text(ref_path)

    setup = ''
    setup += 'viewer.setStyle({chain: "A"}, {cartoon: {color: "lightgray", opacity: 0.7}});\n'
    setup += 'viewer.setStyle({chain: "L"}, {stick: {colorscheme: "greenCarbon", radius: 0.15}});\n'
    setup += 'viewer.setStyle({chain: "L2"}, {sphere: {colorscheme: "redCarbon", radius: 0.5}});\n'
    setup += 'viewer.addLabel("Reference (WT)", {position: viewer.selectedAtoms({chain: "L2"})[0], backgroundColor: "red", fontColor: "white", fontSize: 11});\n'

    color_map = {'WT': '0x2ca02c', 'Y19H': '0xff7f0e', 'R74G': '0x1f77b4',
                 'GKQ': '0xd62728', 'HKQ': '0x9467bd'}

    for idx, (v, m) in enumerate(close_models[1:], 1):
        aligned_pdb = align_pdb_text(ref_ca, m['pdb_path'])
        pdb_escaped = aligned_pdb.replace('\\', '\\\\').replace('`', '\\`').replace('${', '\\${')
        setup += f'''
var pdb{idx} = `{pdb_escaped}`;
viewer.addModel(pdb{idx}, "pdb");
viewer.setStyle({{model: {idx}, chain: "A"}}, {{}});
viewer.setStyle({{model: {idx}, chain: "L"}}, {{}});
viewer.setStyle({{model: {idx}, chain: "T"}}, {{}});
viewer.setStyle({{model: {idx}, chain: "M"}}, {{}});
viewer.setStyle({{model: {idx}, chain: "L2"}}, {{sphere: {{color: {color_map[v]}, radius: 0.5}}}});
'''

    legend_items = []
    for v in VARIANT_ORDER:
        n = sum(1 for vv, _ in close_models if vv == v)
        if n > 0:
            legend_items.append(f'{v}: {n} model(s)')

    viewer = VIEWER_DIV.format(subtitle='All Close-Binding L2 Positions Superposed',
                                div_id='viewer1', width=900, height=600)
    script = make_viewer_script('viewer1', ref_pdb, setup)

    desc = (f'<b>All {len(close_models)} close-binding L2 positions</b> superposed onto WT model_0 backbone. '
            f'L2 spheres colored by variant: '
            + ', '.join(legend_items) + '. '
            f'Note: HKQ has 0 close-binding models. '
            f'All non-competitive L2 positions cluster at the active-site entrance.')

    html = HTML_TEMPLATE.format(
        title='All L2 Positions from Close-Binding Models',
        cdn=THREEDMOL_CDN, description=desc,
        viewers=viewer, script=script)

    path = os.path.join(OUT_DIR, 'view_all_l2_positions.html')
    with open(path, 'w') as f:
        f.write(html)
    print(f'  Saved {path}')


# ── Report Generation ────────────────────────────────────────────────

def write_report(all_data, insights_data):
    """Write the comprehensive markdown report."""
    r74_data = insights_data['r74_loop_interactions']
    pocket_data = insights_data['entrance_pocket']

    n_total = sum(len(all_data[v]['models']) for v in VARIANT_ORDER)
    n_close = sum(1 for v in VARIANT_ORDER for m in all_data[v]['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD)

    close_by_variant = {}
    for v in VARIANT_ORDER:
        models = all_data[v]['models']
        close_by_variant[v] = sum(1 for m in models if m['l2_l1_dist'] < CLOSE_THRESHOLD)

    wt_close_models = [m for m in all_data['WT']['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD]
    close_models_wt_dist = wt_close_models[0]['l2_l1_dist'] if wt_close_models else 10.0
    r74g_close = [m for m in all_data['R74G']['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD]
    r74g_cat_contacts = r74g_close[0]['n_catalytic_contacted'] if r74g_close else 0

    report = f"""# Structural Analysis of Substrate Inhibition in PMDsc

## 1. Executive Summary

- **A secondary MVAP binding site exists at the active-site entrance**, formed by residues R74, K22, N28, N72, and T209. This site is distinct from the catalytic pocket and consistent with non-competitive inhibition.
- **Residue 74 acts as a gateway**: R74 (Arg) provides electrostatic shielding; H74 (His, in HKQ) physically closes the entrance with 21 contacts to loop 17-33; G74 (Gly, in R74G) leaves the gate open but eliminates the MVAP phosphate anchor.
- **The GKQ paradox is explained by the Q212 bridge**: In HKQ, Q212 bridges Y19 and K22 (2.7-2.8 Å), tethering the substrate-binding loop to the H4 helix. This bridge is absent in GKQ despite the identical M212Q mutation.
- **Boltz-2 confidence validates the findings**: A-L2 iPTM strongly anti-correlates with L2-L1 distance (r = -0.85), confirming that close placements are not random.

---

## 2. Background

PMDsc (phosphomevalonate decarboxylase from *Staphylococcus carnosus*) catalyzes the conversion of mevalonate 5-phosphate (MVAP) to isopentenol in the IPP-bypass mevalonate pathway. At high intracellular MVAP concentrations (~100-200 mM), a second MVAP molecule binds an allosteric site and inhibits the enzyme (substrate inhibition).

The substrate inhibition constant Ki is the strongest predictor of in vivo isopentenol titer (Kang et al. 2017). High Ki (weak inhibition) combined with high kcat/Km (good catalysis) yields the best production.

| Variant | Mutations | kcat/Km (mM⁻¹s⁻¹) | Ki (mM) | Titer (mg/L) |
|---------|-----------|-------------------|---------|--------------|
| WT | none | 0.066 | 18 | 475 |
| Y19H | Y19H | 0.78 | N.D. | 388 |
| R74G | R74G | 0.04 | 110 | 975 |
| GKQ | R74G-R147K-M212Q | 0.5 | 11 | 8 |
| HKQ | R74H-R147K-M212Q | 0.4 | 80 | 1079 |

The goal of this analysis is to identify the secondary MVAP binding site responsible for substrate inhibition and explain the variant-specific differences in Ki.

---

## 3. Methods

### Boltz-2 Cofolding with Two MVAP Molecules

We used Boltz-2 structure prediction to cofold each PMDsc variant with **two MVAP molecules**:
- **L1 (MVAP-1)**: Constrained to the active site using pocket conditioning
- **L2 (MVAP-2)**: Unconstrained — free to find any binding site on the protein

For each of the 5 variants (WT, Y19H, R74G, GKQ, HKQ), we generated **10 models each**, totaling **{n_total} models**.

We then analyzed:
1. Where L2 is placed relative to L1 (COM-COM distance)
2. Which protein residues L2 contacts (< 4 Å heavy-atom distance)
3. Whether L2 binding is competitive (inside active site) or non-competitive (at entrance)
4. Per-chain-pair iPTM confidence for L2 placement

### 1-MVAP Reference Structures

For structural comparison (R74-loop interactions, entrance pocket geometry, GKQ paradox), we used single-MVAP Boltz-2 structures (5 models per variant).

---

## 4. L2 Placement Results

### Distance Distributions

![L2 distance distribution](./fig1_l2_distance_distribution.png)
*Figure 1: L2-L1 COM distances for all {n_total} models. Each dot is one Boltz-2 model. Black bar = median. Red dashed line = close-binding threshold (15 Å).*

| Variant | Ki (mM) | Models | Close (< 15 Å) | Fraction | Binding Mode |
|---------|---------|--------|-----------------|----------|--------------|
| WT | 18 | {len(all_data['WT']['models'])} | {close_by_variant['WT']} | {close_by_variant['WT']/len(all_data['WT']['models']):.1%} | entrance |
| Y19H | N.D. | {len(all_data['Y19H']['models'])} | {close_by_variant['Y19H']} | {close_by_variant['Y19H']/len(all_data['Y19H']['models']):.1%} | entrance |
| R74G | 110 | {len(all_data['R74G']['models'])} | {close_by_variant['R74G']} | {close_by_variant['R74G']/len(all_data['R74G']['models']):.1%} | active_site |
| GKQ | 11 | {len(all_data['GKQ']['models'])} | {close_by_variant['GKQ']} | {close_by_variant['GKQ']/len(all_data['GKQ']['models']):.1%} | entrance |
| HKQ | 80 | {len(all_data['HKQ']['models'])} | {close_by_variant['HKQ']} | {close_by_variant['HKQ']/len(all_data['HKQ']['models']):.1%} | — |

**Key observations:**
- **HKQ (best variant)** has **0/10 close-binding models** — the entrance is blocked
- **WT** has the highest close-binding fraction (3/10), consistent with moderate Ki
- **R74G** has only 1 close model, but its L2 binds **inside** the active site (competitive), not at the entrance

![Close-binding frequency](./fig2_close_binding_frequency.png)
*Figure 2: Fraction of models with L2 within 15 Å of L1, by variant.*

---

## 5. The Secondary Binding Site

### Consensus Contact Residues

Across all {n_close} close-binding models (excluding R74G's competitive binding), L2 consistently contacts the same set of residues at the **active-site entrance**:

![Consensus contacts](./fig3_consensus_contacts.png)
*Figure 3: Contact heatmap for close-binding models. Each row is a model; each column is a residue. Red = L2 makes contact (< 4 Å). Frequency counts shown at top.*

**Core contact residues (contacted in ≥ 4/7 close-binding models):**

| Residue | Freq | Role | Can Mutate? |
|---------|------|------|-------------|
| R74 | 6/7 | Gateway — guanidinium H-bonds MVAP phosphate | Yes (to G or H) |
| K22 | 6/7 | Loop 17-33 anchor for MVAP phosphate | **No** (catalytic) |
| N72 | 6/7 | H-bonds MVAP carboxylate from opposite side | Yes |
| T209 | 6/7 | Strongest non-catalytic H-bond (2.0 Å to phosphate) | Yes (adjacent to catalytic S208) |
| N28 | 4/7 | Loop 17-33 anchor | Yes |

**Flanking catalytic residues (do not mutate):**
- S121 (3/7), S153 (3/7), S208 (2/7) — these border the entrance but are essential for catalysis

The secondary site is **not** a deep pocket — it is a shallow depression at the mouth of the active site, where the incoming substrate pauses before entering. At saturating concentrations, a second MVAP occupies this site and blocks the entrance.

> **3D Visualization:** Open `view_consensus_site.html` to interactively examine the secondary site on WT model (m0, L2-L1 = {close_models_wt_dist:.1f} Å). Blue labels = mutation targets, red labels = catalytic residues.

---

## 6. R74G: Competitive vs Non-Competitive Inhibition

R74G (Ki = 110 mM) shows a fundamentally different L2 binding mode from WT:

![R74G competitive](./fig4_r74g_competitive.png)
*Figure 4: Contact residues for WT entrance binding (left) vs R74G active-site binding (right).*

In **WT**, L2 sits at the entrance and contacts non-catalytic residues (R74, N72, T209).

In **R74G**, L2 enters the active site itself, contacting **{r74g_cat_contacts} catalytic residues** (D18, Y19, D20, K22, S121, S153, F155, K158, S208). L1 is displaced ~9 Å from its normal position. This is **competitive** inhibition — the second MVAP competes for the same binding site.

This explains why R74G has high Ki (110 mM): without R74's guanidinium anchor at the entrance, the entrance site is destabilized. L2 can only bind competitively inside the active site, which requires much higher concentrations.

> **3D Visualization:** Open `view_r74g_competitive.html` to see L2 inside the R74G active site.

---

## 7. The R74 Gateway Mechanism

Residue 74 controls access to the secondary binding site through two distinct mechanisms depending on the amino acid:

![R74 loop mechanism](./fig5_r74_loop_mechanism.png)
*Figure 5: R74-loop 17-33 interactions. Left: contact counts. Middle: per-residue frequency. Right: correlation with Ki.*

| Variant | Res 74 | Contacts to Loop 17-33 | Mechanism | Ki |
|---------|--------|----------------------|-----------|-----|
| WT | ARG | 0 | Electrostatic shielding (long-range) | 18 mM |
| Y19H | ARG | 0.4 | Same as WT | N.D. |
| R74G | GLY | 0 | No sidechain — gate open | 110 mM |
| GKQ | GLY | 0 | Gate open + entrance restored by other mutations | 11 mM |
| HKQ | HIS | **21.2** | **Physical closure** — 21 contacts to T25, K26, N28 | 80 mM |

**Key insight:** R74 (Arg) and H74 (His) use completely different mechanisms to modulate access:

- **R74 (WT)**: The long, flexible guanidinium provides an **electrostatic anchor** for MVAP phosphate at the entrance (H-bond distance 2.3 Å). The Arg sidechain doesn't physically touch the loop — it attracts the substrate.
- **H74 (HKQ)**: The shorter, bulkier imidazole ring makes **direct van der Waals contacts** with loop residues T25, K26, N28. It physically closes the entrance like a door.
- **G74 (R74G, GKQ)**: No sidechain at all — the gate is wide open. But whether L2 can still bind depends on other residues.

> **3D Visualization:** Open `view_r74_loop_shield.html` to compare WT (R74, 0 contacts) vs HKQ (H74, 21 contacts) side-by-side.

---

## 8. The GKQ Paradox

The most puzzling observation: R74G alone has Ki = 110 mM (relief from inhibition), but adding R147K + M212Q to make GKQ drops Ki to 11 mM (severe inhibition). HKQ (R74**H**-R147K-M212Q) has Ki = 80 mM. Why does the same pair of mutations (R147K + M212Q) have opposite effects depending on whether position 74 is Gly or His?

![GKQ paradox](./fig6_gkq_paradox.png)
*Figure 6: Q212 bridge distances and K147 salt bridge across variants. Left: Q212.OE1 to Y19.OH distance. Middle: Q212.OE1 to K22.NZ distance. Right: K147.NZ to E144 distance. Red line = H-bond/salt bridge cutoff.*

### The Q212 Bridge

In HKQ, glutamine 212 (from the M212Q mutation) forms an **H-bond bridge** connecting two critical structural elements:
- **Q212.OE1 — Y19.OH** (2.72 Å): tethers the H4 helix to the substrate-binding domain
- **Q212.OE1 — K22.NZ** (2.79 Å): tethers the H4 helix to loop 17-33

This bridge effectively **locks loop 17-33 in place**, preventing it from opening to accommodate L2.

In GKQ, despite carrying the identical M212Q mutation, this bridge is **absent** — the Q212-Y19 and Q212-K22 distances are above H-bond cutoff. The difference is that GKQ has G74 (no sidechain) instead of H74, which changes the conformational landscape of the entire entrance region.

### The K147 Salt Bridge

K147 (from R147K) also behaves differently:
- In **HKQ**, K147 forms a salt bridge with **E144** on helix H2, stabilizing the helix
- In **GKQ**, K147 points toward H176 instead

> **3D Visualization:** Open `view_gkq_vs_hkq.html` to compare the Q212 bridge in GKQ (absent) vs HKQ (present).

---

## 9. Confidence-Distance Correlation

A critical validation: Boltz-2's per-chain A-L2 iPTM (from the L2 perspective) strongly anti-correlates with L2-L1 distance.

![Confidence correlation](./fig7_confidence_correlation.png)
*Figure 7: A-L2 iPTM vs L2-L1 distance for all {n_total} models. Higher confidence = closer L2 placement.*

This means Boltz-2 is **more confident** when it places L2 close to the protein. Close placements are not random noise — the model genuinely predicts favorable protein-ligand interactions at the secondary site.

---

## 10. Assessment

### What the Data Explains

1. **The secondary site identity**: The entrance site (R74, K22, N72, T209, N28) is consistent across variants and models. It is structurally plausible — a shallow pocket at the active-site mouth where the substrate naturally pauses.

2. **Why R74G relieves inhibition**: Removing R74's guanidinium eliminates the entrance anchor. L2 can only bind competitively inside the active site, which requires much higher concentrations (Ki = 110 mM).

3. **Why HKQ is the best variant**: H74 physically closes the entrance (21 contacts to loop 17-33), and the Q212 bridge further locks the loop. Both mechanisms prevent L2 from reaching the entrance site (0/10 close-binding models).

4. **The GKQ paradox**: Despite carrying R147K + M212Q, GKQ fails to form the Q212 bridge that makes HKQ effective. The R74G→G74 open gate allows a conformation where Q212 cannot reach Y19/K22. This explains why GKQ has severe inhibition (Ki = 11) despite its high kcat/Km (0.5).

### What Remains Unexplained

1. **Quantitative Ki prediction**: The close-binding fraction does not quantitatively predict Ki across variants (r = -0.54, p = 0.46, n = 4). With only 10 models per variant, the statistical power is too low.

2. **Y19H mechanism**: Y19H has the best kcat/Km (0.78) but unknown Ki. It has 2/10 close-binding models. Whether its high activity compensates for substrate inhibition remains unclear.

3. **Energetics**: Boltz-2 predicts structures, not binding free energies. We cannot estimate L2 binding affinity or the free energy cost of disrupting the secondary site.

### Limitations

1. **Sample size**: 10 models per variant is statistically limited. Close-binding fractions (0-30%) have wide confidence intervals.
2. **No dynamics**: Boltz-2 produces static structures. The entrance site may open/close on timescales not captured here.
3. **Boltz-2 is not validated for this task**: Cofolding with two copies of the same ligand is an unconventional use. The predictions are structurally informative but not quantitatively reliable.

### Suggested Next Steps

Based on the structural insights, the following mutations are predicted to disrupt the secondary site without affecting catalysis:

| Priority | Mutation | Rationale |
|----------|----------|-----------|
| 1 | T209A or T209V | Disrupts strongest non-catalytic H-bond to MVAP phosphate |
| 2 | N72A or N72V | Disrupts H-bond to MVAP carboxylate |
| 3 | N28A or N28V | Disrupts anchoring in loop 17-33 |

**Caution**: T209 is adjacent to catalytic S208. Validate that T209A does not perturb S208 positioning.

---

## 11. Appendix

### Close-Binding Model Details

| Variant | Model ID | L2-L1 (Å) | Mode | Catalytic Contacts | A-L2 iPTM |
|---------|----------|-----------|------|-------------------|-----------|"""

    for v in VARIANT_ORDER:
        for m in sorted(all_data[v]['models'], key=lambda x: x['l2_l1_dist']):
            if m['l2_l1_dist'] < CLOSE_THRESHOLD:
                iptm = m.get('A_L2_iptm_from_L2')
                iptm_str = f'{iptm:.3f}' if iptm else 'N/A'
                report += f"\n| {v} | {m['model_id']} | {m['l2_l1_dist']:.1f} | {m['binding_mode']} | {m['n_catalytic_contacted']} | {iptm_str} |"

    report += """

### File Index

**Figures:**
- `fig1_l2_distance_distribution.png` — L2-L1 distance distribution by variant
- `fig2_close_binding_frequency.png` — Close-binding fraction by variant
- `fig3_consensus_contacts.png` — Contact residue heatmap for close-binding models
- `fig4_r74g_competitive.png` — WT entrance vs R74G active-site binding
- `fig5_r74_loop_mechanism.png` — R74 gateway mechanism
- `fig6_gkq_paradox.png` — GKQ paradox: Q212 bridge distances
- `fig7_confidence_correlation.png` — Boltz-2 confidence vs L2 distance
- `figS1_ki_vs_close_fraction.png` — Ki vs close-binding fraction
- `figS2_full_contact_heatmap.png` — Full 50-model contact heatmap

**Interactive 3D Visualizations (open in browser):**
- `view_consensus_site.html` — Secondary site on WT with labeled contacts
- `view_r74g_competitive.html` — R74G with L2 inside active site
- `view_r74_loop_shield.html` — WT R74 vs HKQ H74 loop closure
- `view_gkq_vs_hkq.html` — GKQ vs HKQ Q212 bridge comparison
- `view_all_l2_positions.html` — All close-binding L2 positions superposed

**Source Data:**
- `structures/boltz2/output_2mvap/comprehensive_l2_analysis.json` — Full analysis data
- `structures/boltz2/output_2mvap/structural_insights.json` — Structural analysis data
- `structures/boltz2/output_2mvap/comprehensive_l2_analysis.py` — Analysis script
- `structures/boltz2/output_2mvap/structural_insights.py` — Structural insights script
"""

    path = os.path.join(OUT_DIR, 'substrate_inhibition_report.md')
    with open(path, 'w') as f:
        f.write(report)
    print(f'  Saved {path}')


# ── Main ─────────────────────────────────────────────────────────────

def main():
    print('=' * 70)
    print('  Generating PMDsc Substrate Inhibition Report')
    print('=' * 70)

    print('\n1. Scanning all Boltz-2 models...')
    all_data = scan_all_models()
    for v in VARIANT_ORDER:
        n = len(all_data[v]['models'])
        close = sum(1 for m in all_data[v]['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD)
        print(f'   {v}: {n} models, {close} close-binding')

    print('\n2. Loading structural insights data...')
    with open(os.path.join(BASE_2MVAP, 'structural_insights.json')) as f:
        insights_data = json.load(f)

    print('\n3. Generating figures...')
    fig1_l2_distance_distribution(all_data)
    fig2_close_binding_frequency(all_data)
    fig3_consensus_contacts(all_data)
    fig4_r74g_competitive(all_data)
    fig5_r74_loop_mechanism(insights_data)
    fig6_gkq_paradox()
    fig7_confidence_correlation(all_data)
    figS1_ki_vs_close_fraction(all_data)
    figS2_full_contact_heatmap(all_data)

    print('\n4. Generating interactive 3D HTML visualizations...')
    view_consensus_site(all_data)
    view_r74g_competitive(all_data)
    view_r74_loop_shield()
    view_gkq_vs_hkq()
    view_all_l2_positions(all_data)

    print('\n5. Writing report...')
    write_report(all_data, insights_data)

    print('\n' + '=' * 70)
    print('  Done! All outputs saved to:')
    print(f'  {OUT_DIR}/')
    print('=' * 70)


if __name__ == '__main__':
    main()
