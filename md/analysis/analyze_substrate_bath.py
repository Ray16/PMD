#!/usr/bin/env python3
"""Substrate bath MD analysis for PMDsc variants.

Analyzes MD trajectories of PMDsc variants solvated with ~38 free MVAP
molecules (~100 mM) to identify secondary substrate binding sites
relevant to non-competitive substrate inhibition.

Kang et al. 2017 showed that Ki (not kcat/Km) predicts in vivo titer,
and that R74 shields the active site by interacting with loop 17-33.
"""

import argparse
import itertools
import json
import os
import sys
import time
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

import mdtraj as md

# ── Constants ────────────────────────────────────────────────────────

VARIANTS = ['HKQ', 'GKQ', 'R74G', 'Y19H']
VARIANT_INFO = {
    'HKQ':  {'mutations': 'R74H-R147K-M212Q', 'Ki': 80,   'titer': 1079},
    'GKQ':  {'mutations': 'R74G-R147K-M212Q', 'Ki': 11,   'titer': 8},
    'R74G': {'mutations': 'R74G',             'Ki': 110,  'titer': 975},
    'Y19H': {'mutations': 'Y19H',             'Ki': None, 'titer': 388},
}
VARIANT_COLORS = {
    'HKQ': '#2ca02c', 'GKQ': '#d62728', 'R74G': '#1f77b4', 'Y19H': '#ff7f0e',
}

# Biological numbering (1-indexed)
CATALYTIC_BIO = [18, 19, 20, 22, 120, 121, 153, 155, 158, 208, 302]
MUTATION_SITES_BIO = [74, 147, 212]
LOOP_RANGE_BIO = (17, 33)

# 0-indexed for mdtraj
CATALYTIC_IDX = [x - 1 for x in CATALYTIC_BIO]
MUTATION_SITES_IDX = [x - 1 for x in MUTATION_SITES_BIO]

N_PROTEIN_RES = 396
CUTOFFS_NM = [0.4, 0.5, 0.6, 0.8, 1.0]

BASE_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/md'


# ── Data Loading ─────────────────────────────────────────────────────

def load_variant(variant, stride=1):
    """Load trajectory subsetted to protein + CHA residues only."""
    prmtop = f'{BASE_DIR}/md_prep/{variant}/system.prmtop'
    dcd = f'{BASE_DIR}/md_output/{variant}/trajectory.dcd'

    if not os.path.exists(dcd):
        print(f'  WARNING: {dcd} not found, skipping {variant}')
        return None

    print(f'  Loading topology...', flush=True)
    topo = md.load_prmtop(prmtop)

    protein_atoms = topo.select('protein')
    cha_atoms = topo.select('resname CHA')
    atom_subset = np.concatenate([protein_atoms, cha_atoms])

    print(f'  Loading trajectory (stride={stride})...', flush=True)
    traj = md.load(dcd, top=prmtop, atom_indices=atom_subset, stride=stride)
    sub = traj.topology

    print(f'  {traj.n_frames} frames, {sub.n_atoms} atoms (subsetted from {topo.n_atoms})',
          flush=True)

    # Identify protein heavy atoms
    protein_heavy = np.array(sub.select('protein and not element H'))

    # Map protein heavy atoms to residue indices
    atom_to_residue = {}
    residue_to_atoms = {i: [] for i in range(N_PROTEIN_RES)}
    for aidx in protein_heavy:
        ridx = sub.atom(aidx).residue.index
        if ridx < N_PROTEIN_RES:
            atom_to_residue[aidx] = ridx
            residue_to_atoms[ridx].append(aidx)

    # Identify bath MVAP: CHA with 27 atoms, excluding the first (active-site)
    cha_residues = [r for r in sub.residues if r.name == 'CHA']
    cha_27 = [r for r in cha_residues if r.n_atoms == 27]
    active_site_cha = cha_27[0]  # first 27-atom CHA is active-site
    bath_cha = cha_27[1:]        # rest are bath

    bath_heavy = []
    atom_to_mvap = {}
    mvap_to_atoms = {}
    mvap_to_residx = {}
    for mi, res in enumerate(bath_cha):
        atoms = [a.index for a in res.atoms if a.element.symbol != 'H']
        bath_heavy.extend(atoms)
        mvap_to_atoms[mi] = atoms
        mvap_to_residx[mi] = res.index
        for aidx in atoms:
            atom_to_mvap[aidx] = mi
    bath_heavy = np.array(bath_heavy)

    n_bath = len(bath_cha)
    print(f'  Bath MVAP: {n_bath} molecules, {len(bath_heavy)} heavy atoms', flush=True)

    # Image molecules: nearest periodic image of each CHA to protein
    protein_set = set(a for a in sub.atoms if a.residue.is_protein)
    cha_sets = []
    for r in cha_residues:
        cha_sets.append(set(a for a in r.atoms))

    print(f'  Imaging molecules...', flush=True)
    try:
        traj.image_molecules(inplace=True,
                             anchor_molecules=[protein_set],
                             other_molecules=cha_sets)
    except Exception as e:
        print(f'  WARNING: imaging failed ({e}), proceeding without', flush=True)

    # Align on protein CA for density analysis
    ca_atoms = np.array(sub.select('protein and name CA'))
    traj.superpose(traj, frame=0, atom_indices=ca_atoms)

    total_ns = traj.n_frames * traj.timestep / 1000.0 if traj.timestep else traj.n_frames * 10.0 / 1000.0

    return {
        'variant': variant,
        'traj': traj,
        'protein_heavy': protein_heavy,
        'bath_heavy': bath_heavy,
        'atom_to_residue': atom_to_residue,
        'residue_to_atoms': residue_to_atoms,
        'atom_to_mvap': atom_to_mvap,
        'mvap_to_atoms': mvap_to_atoms,
        'mvap_to_residx': mvap_to_residx,
        'n_bath': n_bath,
        'n_frames': traj.n_frames,
        'total_ns': total_ns,
        'active_site_res': active_site_cha.index,
        'ca_atoms': ca_atoms,
    }


# ── Analysis 1: Per-Residue Contact Frequency ───────────────────────

def compute_contact_frequency(data, cutoff_nm):
    """Fraction of frames each protein residue contacts any bath MVAP."""
    traj = data['traj']
    contact = np.zeros((traj.n_frames, N_PROTEIN_RES), dtype=bool)

    neighbors = md.compute_neighbors(
        traj, cutoff_nm, data['bath_heavy'],
        haystack_indices=data['protein_heavy'], periodic=True,
    )
    for fi, nbrs in enumerate(neighbors):
        for aidx in nbrs:
            aidx = int(aidx)
            if aidx in data['atom_to_residue']:
                contact[fi, data['atom_to_residue'][aidx]] = True

    return contact.mean(axis=0)


def run_multi_cutoff_contacts(data):
    """Run contact frequency at multiple cutoffs."""
    results = {}
    for c in CUTOFFS_NM:
        freq = compute_contact_frequency(data, c)
        n_contacting = (freq > 0).sum()
        print(f'    cutoff={c:.1f} nm: {n_contacting} residues with contacts '
              f'(max freq={freq.max():.3f})', flush=True)
        results[c] = freq
    return results


# ── Analysis 2: Per-Residue Minimum Distance ────────────────────────

def compute_min_distance_per_residue(data):
    """Average of per-frame minimum distance from each residue to nearest bath MVAP."""
    traj = data['traj']
    n_frames = traj.n_frames

    # Build residue pairs: each protein residue vs each bath CHA residue
    bath_res_indices = [data['mvap_to_residx'][mi] for mi in range(data['n_bath'])]
    prot_res = list(range(N_PROTEIN_RES))
    pairs = list(itertools.product(prot_res, bath_res_indices))

    print(f'    Computing {len(pairs)} residue-pair distances across {n_frames} frames...',
          flush=True)

    distances, _ = md.compute_contacts(
        traj, contacts=pairs, scheme='closest-heavy',
        ignore_nonprotein=False, periodic=True,
    )
    # Shape: (n_frames, n_pairs)
    dist_3d = distances.reshape(n_frames, N_PROTEIN_RES, data['n_bath'])
    # Min across bath MVAP, then mean across frames
    min_per_frame = dist_3d.min(axis=2)  # (n_frames, 396)
    avg_min = min_per_frame.mean(axis=0)  # (396,)
    std_min = min_per_frame.std(axis=0)

    return avg_min, std_min, min_per_frame


# ── Analysis 3: MVAP Residence / Approach ────────────────────────────

def compute_mvap_approach(data, cutoff_nm=0.8):
    """Per-MVAP minimum distance to protein over time + contact stats."""
    traj = data['traj']
    n_frames = traj.n_frames
    n_bath = data['n_bath']

    # For each MVAP, compute center-of-mass distance to protein COM
    # and minimum heavy-atom distance to protein
    mvap_min_dist = np.full((n_frames, n_bath), np.inf)
    mvap_contact = np.zeros((n_frames, n_bath), dtype=bool)

    # Use compute_neighbors per MVAP for contact detection
    for mi in range(n_bath):
        mvap_atoms = np.array(data['mvap_to_atoms'][mi])
        nbrs = md.compute_neighbors(
            traj, cutoff_nm, mvap_atoms,
            haystack_indices=data['protein_heavy'], periodic=True,
        )
        for fi, n in enumerate(nbrs):
            if len(n) > 0:
                mvap_contact[fi, mi] = True

    # Compute min-distance timeseries via pairwise distances
    # For efficiency, compute COM of each MVAP vs COM of protein per frame
    prot_xyz = traj.xyz[:, data['protein_heavy']]
    prot_com = prot_xyz.mean(axis=1)  # (n_frames, 3)

    for mi in range(n_bath):
        mvap_atoms = np.array(data['mvap_to_atoms'][mi])
        mvap_xyz = traj.xyz[:, mvap_atoms]
        mvap_com = mvap_xyz.mean(axis=1)  # (n_frames, 3)
        # COM-COM distance (approximate, not PBC-corrected, but good for imaging)
        mvap_min_dist[:, mi] = np.linalg.norm(mvap_com - prot_com, axis=1)

    # Residence stats
    dt_ps = 10.0  # 10 ps per frame
    stats = []
    for mi in range(n_bath):
        ct = mvap_contact[:, mi]
        total_contact_ps = ct.sum() * dt_ps
        diffs = np.diff(np.concatenate([[0], ct.astype(int), [0]]))
        starts = np.where(diffs == 1)[0]
        ends = np.where(diffs == -1)[0]
        durations = (ends - starts) * dt_ps if len(starts) > 0 else np.array([])
        stats.append({
            'mvap_idx': mi,
            'total_contact_ps': float(total_contact_ps),
            'n_events': len(starts),
            'max_duration_ps': float(durations.max()) if len(durations) > 0 else 0,
            'min_com_dist_nm': float(mvap_min_dist[:, mi].min()),
        })

    return {
        'min_dist': mvap_min_dist,
        'contact': mvap_contact,
        'stats': stats,
    }


# ── Analysis 4: Binding Site Clustering ──────────────────────────────

def identify_binding_sites(freq_array, traj, min_residues=3, cluster_dist_nm=1.0):
    """Cluster contacting residues into spatial binding sites."""
    # Use top residues by contact frequency
    nonzero = np.where(freq_array > 0)[0]
    if len(nonzero) < min_residues:
        # Fall back: use top 20 residues by frequency (even if zero)
        top_idx = np.argsort(freq_array)[-20:]
        if freq_array[top_idx[-1]] == 0:
            return []
        nonzero = top_idx[freq_array[top_idx] > 0]
        if len(nonzero) < min_residues:
            return []

    # Adaptive threshold: median of nonzero values
    threshold = np.median(freq_array[nonzero])
    hot = np.where(freq_array >= threshold)[0]
    if len(hot) < min_residues:
        hot = nonzero

    # Get CA coordinates
    ca_coords = []
    valid_hot = []
    for ridx in hot:
        res = traj.topology.residue(ridx)
        ca = [a for a in res.atoms if a.name == 'CA']
        if ca:
            ca_coords.append(traj.xyz[0, ca[0].index])
            valid_hot.append(ridx)
    if len(valid_hot) < 2:
        return [{'residues_bio': [int(r + 1) for r in valid_hot],
                 'avg_freq': float(freq_array[valid_hot].mean())}]

    ca_coords = np.array(ca_coords)
    valid_hot = np.array(valid_hot)

    dists = pdist(ca_coords)
    Z = linkage(dists, method='complete')
    labels = fcluster(Z, t=cluster_dist_nm, criterion='distance')

    sites = []
    for lid in np.unique(labels):
        mask = labels == lid
        members = valid_hot[mask]
        bio = [int(r + 1) for r in members]

        # Check overlap with key regions
        overlaps = []
        for r in members:
            if r in CATALYTIC_IDX:
                overlaps.append('catalytic')
            if r in MUTATION_SITES_IDX:
                overlaps.append('mutation_site')
            if LOOP_RANGE_BIO[0] - 1 <= r <= LOOP_RANGE_BIO[1] - 1:
                overlaps.append('loop_17-33')
        overlaps = list(set(overlaps))

        sites.append({
            'residues_bio': bio,
            'avg_freq': float(freq_array[members].mean()),
            'max_freq': float(freq_array[members].max()),
            'overlaps': overlaps,
        })

    return sorted(sites, key=lambda s: -s['avg_freq'])


# ── Analysis 5: 3D Density ───────────────────────────────────────────

def compute_density(data, spacing_nm=0.1, padding_nm=2.0):
    """3D density of bath MVAP heavy atoms around the protein."""
    traj = data['traj']
    prot_xyz = traj.xyz[:, data['protein_heavy']]
    prot_min = prot_xyz.min(axis=(0, 1)) - padding_nm
    prot_max = prot_xyz.max(axis=(0, 1)) + padding_nm

    n_bins = np.ceil((prot_max - prot_min) / spacing_nm).astype(int)

    bath_xyz = traj.xyz[:, data['bath_heavy']]
    bath_flat = bath_xyz.reshape(-1, 3)

    density, edges = np.histogramdd(
        bath_flat, bins=list(n_bins),
        range=list(zip(prot_min, prot_max)),
    )
    voxel_vol = spacing_nm ** 3
    density /= (traj.n_frames * voxel_vol)

    return {
        'density': density,
        'origin': prot_min,
        'spacing': spacing_nm,
        'n_bins': n_bins,
    }


def write_dx(density_data, filepath):
    """Write density in OpenDX format (coordinates in Angstroms)."""
    d = density_data
    nx, ny, nz = d['n_bins']
    ox, oy, oz = d['origin'] * 10  # nm -> A
    sp = d['spacing'] * 10

    with open(filepath, 'w') as f:
        f.write(f'object 1 class gridpositions counts {nx} {ny} {nz}\n')
        f.write(f'origin {ox:.4f} {oy:.4f} {oz:.4f}\n')
        f.write(f'delta {sp:.4f} 0.0000 0.0000\n')
        f.write(f'delta 0.0000 {sp:.4f} 0.0000\n')
        f.write(f'delta 0.0000 0.0000 {sp:.4f}\n')
        f.write(f'object 2 class gridconnections counts {nx} {ny} {nz}\n')
        n_total = int(nx * ny * nz)
        f.write(f'object 3 class array type double rank 0 items {n_total} data follows\n')
        flat = d['density'].flatten()
        for i in range(0, len(flat), 3):
            chunk = flat[i:i + 3]
            f.write(' '.join(f'{v:.6e}' for v in chunk) + '\n')
        f.write('attribute "dep" string "positions"\n')
        f.write('object "density" class field\n')
        f.write('component "positions" value 1\n')
        f.write('component "connections" value 2\n')
        f.write('component "data" value 3\n')


# ── Plotting ─────────────────────────────────────────────────────────

def _annotate_residue_axis(ax, ymin=None, ymax=None):
    """Add catalytic residue, mutation site, and loop annotations."""
    if ymin is None:
        ymin, ymax = ax.get_ylim()
    # Loop 17-33 background
    ax.axvspan(LOOP_RANGE_BIO[0], LOOP_RANGE_BIO[1],
               alpha=0.1, color='green', zorder=0)
    # Catalytic residues
    for r in CATALYTIC_BIO:
        ax.axvline(r, color='red', alpha=0.3, lw=0.5, zorder=1)
    # Mutation sites
    for r in MUTATION_SITES_BIO:
        ax.axvline(r, color='blue', alpha=0.5, lw=1.0, ls='--', zorder=1)


def plot_min_distance_heatmap(all_min_dist, output_dir):
    """Heatmap of per-residue average minimum distance to bath MVAP."""
    variants = list(all_min_dist.keys())
    n_var = len(variants)
    residues_bio = np.arange(1, N_PROTEIN_RES + 1)

    fig, axes = plt.subplots(n_var + 1, 1, figsize=(20, 3 * n_var + 4),
                             gridspec_kw={'height_ratios': [1] * n_var + [2]})

    vmin = min(d.min() for d in all_min_dist.values())
    vmax = np.percentile(np.concatenate(list(all_min_dist.values())), 95)

    for i, var in enumerate(variants):
        ax = axes[i]
        dist = all_min_dist[var]
        im = ax.imshow(dist.reshape(1, -1), aspect='auto', cmap='RdYlBu',
                       extent=[0.5, N_PROTEIN_RES + 0.5, 0, 1],
                       vmin=vmin, vmax=vmax)
        ki = VARIANT_INFO[var]['Ki']
        ki_str = f'Ki={ki}' if ki else 'Ki=N.D.'
        ax.set_ylabel(f'{var}\n({ki_str})', fontsize=10, rotation=0,
                      labelpad=60, va='center')
        ax.set_yticks([])
        if i < n_var - 1:
            ax.set_xticks([])
        _annotate_residue_axis(ax, 0, 1)

    # Colorbar
    cbar = fig.colorbar(im, ax=axes[:n_var], orientation='vertical',
                        fraction=0.01, pad=0.02)
    cbar.set_label('Avg min distance to bath MVAP (nm)')

    # Bottom panel: overlay line plot
    ax_bot = axes[-1]
    for var in variants:
        dist = all_min_dist[var]
        ax_bot.plot(residues_bio, dist, label=f'{var}', color=VARIANT_COLORS[var],
                    alpha=0.8, lw=0.8)
    _annotate_residue_axis(ax_bot)
    ax_bot.set_xlabel('Residue (biological numbering)')
    ax_bot.set_ylabel('Avg min distance (nm)')
    ax_bot.set_xlim(0.5, N_PROTEIN_RES + 0.5)
    ax_bot.legend(fontsize=8, loc='upper right')

    # Legend for annotations
    legend_elements = [
        Patch(facecolor='green', alpha=0.15, label='Loop 17-33'),
        plt.Line2D([0], [0], color='red', alpha=0.4, lw=1, label='Catalytic residues'),
        plt.Line2D([0], [0], color='blue', ls='--', lw=1, label='Mutation sites (74,147,212)'),
    ]
    ax_bot.legend(handles=legend_elements + [
        plt.Line2D([0], [0], color=VARIANT_COLORS[v], lw=1.5, label=v)
        for v in variants
    ], fontsize=7, loc='upper right', ncol=2)

    fig.suptitle('Per-Residue Average Minimum Distance to Bath MVAP', fontsize=14, y=0.98)
    plt.tight_layout(rect=[0, 0, 0.97, 0.96])
    path = os.path.join(output_dir, 'min_distance_heatmap.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def plot_contact_heatmap(all_contacts, output_dir):
    """Multi-cutoff contact frequency comparison."""
    variants = list(all_contacts.keys())
    cutoffs = sorted(list(all_contacts[variants[0]].keys()))
    residues_bio = np.arange(1, N_PROTEIN_RES + 1)

    # One figure per cutoff with significant contacts
    for cutoff in cutoffs:
        max_freq = max(all_contacts[v][cutoff].max() for v in variants)
        if max_freq < 1e-6:
            continue

        n_var = len(variants)
        fig, axes = plt.subplots(n_var, 1, figsize=(20, 2.5 * n_var), sharex=True)
        if n_var == 1:
            axes = [axes]

        for i, var in enumerate(variants):
            freq = all_contacts[var][cutoff]
            axes[i].bar(residues_bio, freq, width=1.0, color=VARIANT_COLORS[var], alpha=0.7)
            ki = VARIANT_INFO[var]['Ki']
            ki_str = f'Ki={ki}' if ki else 'Ki=N.D.'
            axes[i].set_ylabel(f'{var}\n({ki_str})', fontsize=9, rotation=0,
                               labelpad=55, va='center')
            _annotate_residue_axis(axes[i])
            axes[i].set_xlim(0.5, N_PROTEIN_RES + 0.5)

        axes[-1].set_xlabel('Residue (biological numbering)')
        fig.suptitle(f'Per-Residue Contact Frequency (cutoff = {cutoff:.1f} nm)',
                     fontsize=13, y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        path = os.path.join(output_dir, f'contact_freq_{cutoff:.1f}nm.png')
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved {path}')


def plot_gkq_vs_r74g(all_min_dist, all_contacts, output_dir):
    """Focused comparison of GKQ vs R74G."""
    if 'GKQ' not in all_min_dist or 'R74G' not in all_min_dist:
        print('  Skipping GKQ vs R74G comparison (missing data)')
        return

    residues_bio = np.arange(1, N_PROTEIN_RES + 1)
    fig, axes = plt.subplots(3, 1, figsize=(20, 10), sharex=True)

    # Panel 1: Min distance overlay
    ax = axes[0]
    ax.plot(residues_bio, all_min_dist['R74G'], color=VARIANT_COLORS['R74G'],
            label='R74G (Ki=110, titer=975)', lw=1)
    ax.plot(residues_bio, all_min_dist['GKQ'], color=VARIANT_COLORS['GKQ'],
            label='GKQ (Ki=11, titer=8)', lw=1)
    _annotate_residue_axis(ax)
    ax.set_ylabel('Avg min distance (nm)')
    ax.legend(fontsize=9)
    ax.set_title('GKQ vs R74G: Where does MVAP bind differently?', fontsize=12)

    # Panel 2: Difference (GKQ - R74G), negative = MVAP closer in GKQ
    ax = axes[1]
    diff = all_min_dist['GKQ'] - all_min_dist['R74G']
    colors = ['#d62728' if d < 0 else '#1f77b4' for d in diff]
    ax.bar(residues_bio, diff, width=1.0, color=colors, alpha=0.7)
    _annotate_residue_axis(ax)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_ylabel('Distance diff (nm)\nGKQ - R74G')
    ax.set_title('Negative (red) = MVAP closer in GKQ (low Ki, bad)', fontsize=10)

    # Panel 3: Zoom on mutation neighborhoods (residues 60-90, 135-160, 200-225)
    ax = axes[2]
    regions = [(60, 90, 'R74 neighborhood'), (135, 160, 'R147 neighborhood'),
               (200, 225, 'M212 neighborhood')]
    x_all = []
    d_all = []
    labels_all = []
    for start, end, label in regions:
        mask = (residues_bio >= start) & (residues_bio <= end)
        x_all.extend(residues_bio[mask])
        d_all.extend(diff[mask])
        labels_all.extend([label] * mask.sum())

    region_colors = {'R74 neighborhood': '#9467bd',
                     'R147 neighborhood': '#8c564b',
                     'M212 neighborhood': '#e377c2'}
    for start, end, label in regions:
        mask = (residues_bio >= start) & (residues_bio <= end)
        ax.bar(residues_bio[mask], diff[mask], width=1.0,
               color=region_colors[label], alpha=0.7, label=label)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_ylabel('Distance diff (nm)')
    ax.set_xlabel('Residue (biological numbering)')
    ax.legend(fontsize=9)
    ax.set_xlim(55, 230)
    _annotate_residue_axis(ax)

    plt.tight_layout()
    path = os.path.join(output_dir, 'gkq_vs_r74g.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def plot_mvap_approach(all_approach, output_dir):
    """Per-MVAP approach distance distributions."""
    variants = list(all_approach.keys())
    fig, axes = plt.subplots(1, len(variants), figsize=(5 * len(variants), 4),
                             sharey=True)
    if len(variants) == 1:
        axes = [axes]

    for i, var in enumerate(variants):
        stats = all_approach[var]['stats']
        min_dists = [s['min_com_dist_nm'] for s in stats]
        axes[i].hist(min_dists, bins=20, color=VARIANT_COLORS[var], alpha=0.7,
                     edgecolor='black', lw=0.5)
        ki = VARIANT_INFO[var]['Ki']
        ki_str = f'Ki={ki}' if ki else 'Ki=N.D.'
        axes[i].set_title(f'{var} ({ki_str})', fontsize=11)
        axes[i].set_xlabel('Min COM-COM distance (nm)')
        axes[i].axvline(np.median(min_dists), color='black', ls='--', lw=1,
                        label=f'median={np.median(min_dists):.2f}')
        axes[i].legend(fontsize=8)

    axes[0].set_ylabel('Count (MVAP molecules)')
    fig.suptitle('Closest Approach of Each Bath MVAP to Protein', fontsize=13, y=1.02)
    plt.tight_layout()
    path = os.path.join(output_dir, 'mvap_approach.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


def plot_closest_mvap_residues(all_min_dist_per_frame, all_data, output_dir):
    """Which residues are the closest 'landing spots' for bath MVAP?"""
    variants = list(all_min_dist_per_frame.keys())
    residues_bio = np.arange(1, N_PROTEIN_RES + 1)

    fig, axes = plt.subplots(len(variants), 1, figsize=(20, 3 * len(variants)),
                             sharex=True)
    if len(variants) == 1:
        axes = [axes]

    for i, var in enumerate(variants):
        # For each frame, which residue has the smallest min-distance?
        min_per_frame = all_min_dist_per_frame[var]  # (n_frames, 396)
        # Count how often each residue is the closest
        closest_residue = min_per_frame.argmin(axis=1)  # (n_frames,)
        counts = np.bincount(closest_residue, minlength=N_PROTEIN_RES)
        freq = counts / len(closest_residue)

        axes[i].bar(residues_bio, freq, width=1.0, color=VARIANT_COLORS[var], alpha=0.7)
        ki = VARIANT_INFO[var]['Ki']
        ki_str = f'Ki={ki}' if ki else 'Ki=N.D.'
        axes[i].set_ylabel(f'{var}\n({ki_str})', fontsize=9, rotation=0,
                           labelpad=55, va='center')
        _annotate_residue_axis(axes[i])

        # Label top 5
        top5 = np.argsort(freq)[-5:][::-1]
        for ridx in top5:
            if freq[ridx] > 0.01:
                res = all_data[var]['traj'].topology.residue(ridx)
                axes[i].annotate(f'{res.name}{ridx + 1}',
                                 xy=(ridx + 1, freq[ridx]),
                                 fontsize=6, ha='center', va='bottom')

    axes[-1].set_xlabel('Residue (biological numbering)')
    axes[-1].set_xlim(0.5, N_PROTEIN_RES + 0.5)
    fig.suptitle('Fraction of Frames Each Residue is the Closest to Any Bath MVAP',
                 fontsize=13, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    path = os.path.join(output_dir, 'closest_residue_freq.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {path}')


# ── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variants', nargs='+', default=VARIANTS)
    parser.add_argument('--stride', type=int, default=1)
    parser.add_argument('--output-dir', default=f'{BASE_DIR}/analysis')
    args = parser.parse_args()

    out = args.output_dir
    os.makedirs(out, exist_ok=True)

    # ── Load ─────────────────────────────────────────────────────────
    print('=' * 60)
    print('  Loading trajectories')
    print('=' * 60)
    all_data = {}
    for var in args.variants:
        print(f'\n[{var}]')
        data = load_variant(var, stride=args.stride)
        if data:
            all_data[var] = data

    if not all_data:
        print('No trajectories loaded. Exiting.')
        sys.exit(1)

    # ── 1. Multi-cutoff contact frequency ────────────────────────────
    print('\n' + '=' * 60)
    print('  Analysis 1: Per-residue contact frequency (multi-cutoff)')
    print('=' * 60)
    all_contacts = {}
    for var, data in all_data.items():
        print(f'\n[{var}] ({data["n_frames"]} frames, {data["total_ns"]:.2f} ns)')
        all_contacts[var] = run_multi_cutoff_contacts(data)

    # ── 2. Per-residue minimum distance ──────────────────────────────
    print('\n' + '=' * 60)
    print('  Analysis 2: Per-residue minimum distance')
    print('=' * 60)
    all_min_dist = {}
    all_min_dist_std = {}
    all_min_dist_per_frame = {}
    for var, data in all_data.items():
        print(f'\n[{var}]')
        avg, std, per_frame = compute_min_distance_per_residue(data)
        all_min_dist[var] = avg
        all_min_dist_std[var] = std
        all_min_dist_per_frame[var] = per_frame
        closest_res = np.argmin(avg)
        res_name = data['traj'].topology.residue(closest_res).name
        print(f'    Closest residue: {res_name}{closest_res + 1} '
              f'(avg min dist = {avg[closest_res]:.3f} nm)')

    # ── 3. MVAP approach analysis ────────────────────────────────────
    print('\n' + '=' * 60)
    print('  Analysis 3: MVAP approach / residence')
    print('=' * 60)
    all_approach = {}
    for var, data in all_data.items():
        print(f'\n[{var}]')
        all_approach[var] = compute_mvap_approach(data, cutoff_nm=0.8)
        stats = all_approach[var]['stats']
        n_contacting = sum(1 for s in stats if s['total_contact_ps'] > 0)
        closest = min(s['min_com_dist_nm'] for s in stats)
        print(f'    {n_contacting}/{len(stats)} MVAPs contacted protein (0.8 nm cutoff)')
        print(f'    Closest approach: {closest:.3f} nm')

    # ── 4. Binding site clustering ───────────────────────────────────
    print('\n' + '=' * 60)
    print('  Analysis 4: Binding site clustering')
    print('=' * 60)
    all_sites = {}
    # Use the best cutoff that shows contacts
    for var, data in all_data.items():
        best_cutoff = None
        for c in sorted(CUTOFFS_NM, reverse=True):
            if all_contacts[var][c].max() > 0:
                best_cutoff = c
                break
        if best_cutoff:
            freq = all_contacts[var][best_cutoff]
            sites = identify_binding_sites(freq, data['traj'])
            all_sites[var] = {'cutoff': best_cutoff, 'sites': sites}
            print(f'\n[{var}] (cutoff={best_cutoff} nm): {len(sites)} binding site clusters')
            for j, s in enumerate(sites[:5]):
                print(f'    Site {j + 1}: residues {s["residues_bio"]}, '
                      f'avg_freq={s["avg_freq"]:.3f}, overlaps={s["overlaps"]}')
        else:
            all_sites[var] = {'cutoff': None, 'sites': []}
            print(f'\n[{var}] No contacts at any cutoff')

    # ── 5. 3D density ────────────────────────────────────────────────
    print('\n' + '=' * 60)
    print('  Analysis 5: 3D density grids')
    print('=' * 60)
    for var, data in all_data.items():
        print(f'\n[{var}] Computing density grid...')
        dens = compute_density(data)
        dx_path = os.path.join(out, f'density_{var}.dx')
        write_dx(dens, dx_path)
        npz_path = os.path.join(out, f'density_{var}.npz')
        np.savez_compressed(npz_path, density=dens['density'],
                            origin=dens['origin'], spacing=dens['spacing'])
        max_dens = dens['density'].max()
        print(f'    Grid: {dens["n_bins"]}, max density: {max_dens:.2f} atoms/nm^3')
        print(f'    Saved {dx_path}')

    # ── Figures ──────────────────────────────────────────────────────
    print('\n' + '=' * 60)
    print('  Generating figures')
    print('=' * 60)

    print('\n  Min distance heatmap...')
    plot_min_distance_heatmap(all_min_dist, out)

    print('\n  Contact frequency heatmaps...')
    plot_contact_heatmap(all_contacts, out)

    print('\n  GKQ vs R74G comparison...')
    plot_gkq_vs_r74g(all_min_dist, all_contacts, out)

    print('\n  MVAP approach distributions...')
    plot_mvap_approach(all_approach, out)

    print('\n  Closest residue frequency...')
    plot_closest_mvap_residues(all_min_dist_per_frame, all_data, out)

    # ── Save results ─────────────────────────────────────────────────
    print('\n' + '=' * 60)
    print('  Saving results')
    print('=' * 60)

    results = {
        'variants': list(all_data.keys()),
        'trajectory_info': {
            var: {'n_frames': d['n_frames'], 'total_ns': round(d['total_ns'], 3)}
            for var, d in all_data.items()
        },
        'min_distance_per_residue': {
            var: {
                'top10_closest': [
                    {
                        'residue_bio': int(idx + 1),
                        'residue_name': all_data[var]['traj'].topology.residue(idx).name,
                        'avg_min_dist_nm': round(float(all_min_dist[var][idx]), 4),
                    }
                    for idx in np.argsort(all_min_dist[var])[:10]
                ]
            }
            for var in all_data
        },
        'contact_frequency': {
            var: {
                str(c): {
                    'n_contacting_residues': int((all_contacts[var][c] > 0).sum()),
                    'max_freq': round(float(all_contacts[var][c].max()), 4),
                }
                for c in CUTOFFS_NM
            }
            for var in all_data
        },
        'binding_sites': {
            var: all_sites.get(var, {}).get('sites', [])
            for var in all_data
        },
        'mvap_approach': {
            var: {
                'n_contacting_0.8nm': sum(1 for s in all_approach[var]['stats']
                                          if s['total_contact_ps'] > 0),
                'closest_approach_nm': round(
                    min(s['min_com_dist_nm'] for s in all_approach[var]['stats']), 4),
            }
            for var in all_data
        },
    }

    json_path = os.path.join(out, 'results.json')
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'  Saved {json_path}')

    # ── Summary ──────────────────────────────────────────────────────
    print('\n' + '=' * 60)
    print('  SUMMARY')
    print('=' * 60)
    print(f'\n  {"Variant":<8} {"Frames":>6} {"Time":>7} {"Ki":>5} {"Titer":>6} '
          f'{"Closest res":>20} {"Min dist":>9} {"Contacts@1nm":>13}')
    print('  ' + '-' * 85)
    for var in all_data:
        d = all_data[var]
        ki = VARIANT_INFO[var]['Ki']
        ki_str = str(ki) if ki else 'N.D.'
        closest_idx = np.argmin(all_min_dist[var])
        closest_name = d['traj'].topology.residue(closest_idx).name
        closest_bio = closest_idx + 1
        n_contacts_1nm = int((all_contacts[var][1.0] > 0).sum()) if 1.0 in all_contacts[var] else 0
        print(f'  {var:<8} {d["n_frames"]:>6} {d["total_ns"]:>6.2f}ns {ki_str:>5} '
              f'{VARIANT_INFO[var]["titer"]:>6} '
              f'{closest_name}{closest_bio:>4} ({all_min_dist[var][closest_idx]:.3f}nm) '
              f'{n_contacts_1nm:>13}')

    # Top-5 closest residues comparison
    print(f'\n  Top-5 closest residues per variant:')
    for var in all_data:
        top5 = np.argsort(all_min_dist[var])[:5]
        labels = []
        for idx in top5:
            name = all_data[var]['traj'].topology.residue(idx).name
            labels.append(f'{name}{idx + 1}({all_min_dist[var][idx]:.2f})')
        print(f'    {var:<6}: {", ".join(labels)}')

    # GKQ vs R74G key differences
    if 'GKQ' in all_min_dist and 'R74G' in all_min_dist:
        diff = all_min_dist['GKQ'] - all_min_dist['R74G']
        closer_in_gkq = np.argsort(diff)[:10]
        print(f'\n  Residues where MVAP is CLOSER in GKQ vs R74G (potential inhibitory site):')
        for idx in closer_in_gkq:
            if diff[idx] < 0:
                name = all_data['GKQ']['traj'].topology.residue(idx).name
                print(f'    {name}{idx + 1}: GKQ={all_min_dist["GKQ"][idx]:.3f} nm, '
                      f'R74G={all_min_dist["R74G"][idx]:.3f} nm, '
                      f'diff={diff[idx]:.3f} nm')

    print('\n  Done.')


if __name__ == '__main__':
    main()
