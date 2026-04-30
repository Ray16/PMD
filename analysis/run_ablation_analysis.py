#!/usr/bin/env python3
"""Analysis of secondary MVAP (S) placement from Boltz-2 cofolding ablation runs.

Handles both no_cofactor and no_constraint conditions.
Scans both main variant dirs and seed dirs for PDB models.
"""
import json
import os
import sys
import glob
import numpy as np
from collections import defaultdict, Counter
from scipy import stats

VARIANT_NAMES = {
    'WT': 'PMDsc_WT_2mvap',
    'Y19H': 'PMDsc_Y19H_2mvap',
    'K22M': 'PMDsc_K22M_2mvap',
    'K22Y': 'PMDsc_K22Y_2mvap',
    'R74G': 'PMDsc_R74G_2mvap',
    'R74H': 'PMDsc_R74H_2mvap',
    'I145A': 'PMDsc_I145A_2mvap',
    'I145F': 'PMDsc_I145F_2mvap',
    'R147K': 'PMDsc_R147K_2mvap',
    'S186C': 'PMDsc_S186C_2mvap',
    'S208E': 'PMDsc_S208E_2mvap',
    'T209D': 'PMDsc_T209D_2mvap',
    'M212Q': 'PMDsc_M212Q_2mvap',
    'I226V': 'PMDsc_I226V_2mvap',
    'V230E': 'PMDsc_V230E_2mvap',
    'R74G-R147K': 'PMDsc_R74G_R147K_2mvap',
    'HKQ': 'PMDsc_R74H_R147K_M212Q_2mvap',
    'GKQ': 'PMDsc_R74G_R147K_M212Q_2mvap',
    'R74G-R147K-Q140L': 'PMDsc_R74G_R147K_Q140L_2mvap',
}

VARIANT_DATA = {
    'WT':   {'ki': 18,   'kcat': 0.15, 'km': 2.3,  'kcat_km': 0.066, 'titer': 475},
    'Y19H': {'ki': None, 'kcat': 0.27, 'km': 0.35, 'kcat_km': 0.78,  'titer': 388},
    'K22M': {'ki': None, 'kcat': 0.09, 'km': 1.3,  'kcat_km': 0.12,  'titer': 22},
    'K22Y': {'ki': None, 'kcat': 0.09, 'km': 1.3,  'kcat_km': 0.12,  'titer': 22},
    'R74G': {'ki': 110,  'kcat': 0.14, 'km': 3.4,  'kcat_km': 0.04,  'titer': 975},
    'R74H': {'ki': None, 'kcat': 0.33, 'km': 0.75, 'kcat_km': 0.44,  'titer': 770},
    'I145A': {'ki': None, 'kcat': 0.029,'km': 2.0,  'kcat_km': 0.01,  'titer': 623},
    'I145F': {'ki': None, 'kcat': 0.28, 'km': 1.36, 'kcat_km': 0.20,  'titer': None},
    'R147K': {'ki': None, 'kcat': 0.149,'km': 0.5,  'kcat_km': 0.32,  'titer': 793},
    'S186C': {'ki': None, 'kcat': 0.07, 'km': 0.8,  'kcat_km': 0.08,  'titer': 596},
    'S208E': {'ki': None, 'kcat': None, 'km': None, 'kcat_km': None,  'titer': None},
    'T209D': {'ki': None, 'kcat': 0.13, 'km': 0.99, 'kcat_km': 0.13,  'titer': None},
    'M212Q': {'ki': None, 'kcat': 0.35, 'km': 0.7,  'kcat_km': 0.50,  'titer': 601},
    'I226V': {'ki': None, 'kcat': 0.16, 'km': 0.34, 'kcat_km': 0.46,  'titer': 633},
    'V230E': {'ki': 10,   'kcat': 0.07, 'km': 0.8,  'kcat_km': 0.08,  'titer': 278},
    'R74G-R147K': {'ki': None, 'kcat': 0.22, 'km': 0.53, 'kcat_km': 0.42, 'titer': 909},
    'HKQ':  {'ki': 80,   'kcat': 0.16, 'km': 0.40, 'kcat_km': 0.40,  'titer': 1079},
    'GKQ':  {'ki': 11,   'kcat': 0.22, 'km': 0.50, 'kcat_km': 0.50,  'titer': 8},
    'R74G-R147K-Q140L': {'ki': None, 'kcat': 0.06, 'km': 2.0, 'kcat_km': 0.03, 'titer': 401},
}

VARIANT_ORDER = list(VARIANT_NAMES.keys())
CATALYTIC = {18, 19, 20, 22, 120, 121, 153, 155, 158, 208, 302}
CLOSE_THRESHOLD = 15.0
INTERMEDIATE_THRESHOLD = 25.0


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


def compute_l2_metrics(chains):
    if 'L' not in chains or 'S' not in chains or 'A' not in chains:
        return None

    l1_coords = np.array([a['xyz'] for a in chains['L']])
    s_coords = np.array([a['xyz'] for a in chains['S']])
    l1_com = l1_coords.mean(axis=0)
    s_com = s_coords.mean(axis=0)
    l2_l1_dist = float(np.linalg.norm(s_com - l1_com))

    ca_coords = []
    ca_resnums = []
    for a in chains['A']:
        if a['atom'] == 'CA':
            ca_coords.append(a['xyz'])
            ca_resnums.append(a['resnum'])
    ca_coords = np.array(ca_coords)

    dists_to_l2 = np.linalg.norm(ca_coords - s_com, axis=1)
    closest_idx = np.argsort(dists_to_l2)[:10]
    closest_residues = [(ca_resnums[i], float(dists_to_l2[i])) for i in closest_idx]

    s_heavy = [a for a in chains['S'] if a['is_heavy']]
    prot_heavy = [a for a in chains['A'] if a['is_heavy']]
    contacts_4a = {}
    for pa in prot_heavy:
        for la in s_heavy:
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
                        for la in s_heavy for ca in cat_atoms)
            l2_to_catalytic_dists[rn] = float(min_d)

    n_catalytic_contacted = sum(1 for d in l2_to_catalytic_dists.values() if d < 4.0)

    if l2_l1_dist > CLOSE_THRESHOLD:
        if l2_l1_dist > INTERMEDIATE_THRESHOLD:
            binding_mode = 'distant'
        else:
            binding_mode = 'intermediate'
    elif n_catalytic_contacted >= 5:
        binding_mode = 'active_site'
    else:
        binding_mode = 'entrance'

    # Also compute L1 metrics (distance to active site)
    # Use catalytic residue CA positions as active site reference
    active_site_cas = []
    for a in chains['A']:
        if a['atom'] == 'CA' and a['resnum'] in {19, 22, 120, 121, 153}:
            active_site_cas.append(a['xyz'])
    if active_site_cas:
        as_center = np.mean(active_site_cas, axis=0)
        l1_to_active_site = float(np.linalg.norm(l1_com - as_center))
    else:
        l1_to_active_site = None

    return {
        'l2_l1_dist': l2_l1_dist,
        'closest_residues': closest_residues,
        'contacts_4a': {str(k): v for k, v in contacts_4a.items()},
        'n_contacts_4a': len(contacts_4a),
        'n_catalytic_contacted': n_catalytic_contacted,
        'catalytic_dists': {str(k): round(v, 2) for k, v in l2_to_catalytic_dists.items()},
        'binding_mode': binding_mode,
        'l1_to_active_site': l1_to_active_site,
    }


def find_all_pdbs(base_dir, dirname):
    """Find all PDB models from the merged directory."""
    pdbs = []
    pred_dir = os.path.join(base_dir, dirname, f'boltz_results_{dirname}',
                            'predictions', dirname)
    if os.path.isdir(pred_dir):
        pdbs.extend(sorted(glob.glob(os.path.join(pred_dir, f'{dirname}_model_*.pdb'))))
    return pdbs


def get_confidence(conf_path, s_chain_idx):
    """Extract confidence metrics with correct chain index for S."""
    if not os.path.exists(conf_path):
        return {}
    with open(conf_path) as f:
        conf_data = json.load(f)
    s_idx = str(s_chain_idx)
    return {
        'confidence_score': conf_data.get('confidence_score'),
        'iptm': conf_data.get('iptm'),
        'ligand_iptm': conf_data.get('ligand_iptm'),
        'A_S_iptm_from_S': conf_data.get('pair_chains_iptm', {}).get(s_idx, {}).get('0'),
        'A_S_iptm_from_A': conf_data.get('pair_chains_iptm', {}).get('0', {}).get(s_idx),
        'S_ptm': conf_data.get('chains_ptm', {}).get(s_idx),
        'A_L1_iptm_from_L1': conf_data.get('pair_chains_iptm', {}).get('1', {}).get('0'),
        'L1_ptm': conf_data.get('chains_ptm', {}).get('1'),
    }


def scan_all_models(base_dir, s_chain_idx):
    all_data = {}
    for short, dirname in VARIANT_NAMES.items():
        pdbs = find_all_pdbs(base_dir, dirname)
        if not pdbs:
            print(f'  WARNING: no PDBs found for {short}')
            continue

        models = []
        for pdb_path in pdbs:
            chains = parse_pdb_minimal(pdb_path)
            metrics = compute_l2_metrics(chains)
            if not metrics:
                continue

            model_name = os.path.basename(pdb_path).replace('.pdb', '')
            conf_path = os.path.join(os.path.dirname(pdb_path),
                                      f'confidence_{model_name}.json')
            conf = get_confidence(conf_path, s_chain_idx)

            source = 'seed' if '_seed' in pdb_path else 'main'
            model_idx = model_name.split('_')[-1]
            models.append({
                'model_id': f'{source}_m{model_idx}',
                'pdb_path': pdb_path,
                **metrics,
                **conf,
            })

        vdata = VARIANT_DATA.get(short, {})
        all_data[short] = {
            'variant': short,
            'ki': vdata.get('ki'),
            'titer': vdata.get('titer'),
            'kcat_km': vdata.get('kcat_km'),
            'n_models': len(models),
            'models': models,
        }
        print(f'  {short}: {len(models)} models')
    return all_data


def analyze_distributions(all_data):
    stats_table = {}
    for short in VARIANT_ORDER:
        vd = all_data.get(short)
        if not vd or not vd['models']:
            continue
        dists = [m['l2_l1_dist'] for m in vd['models']]
        n = len(dists)
        modes = Counter(m['binding_mode'] for m in vd['models'])
        n_close = sum(1 for d in dists if d < CLOSE_THRESHOLD)
        n_entrance = modes.get('entrance', 0)
        n_active_site = modes.get('active_site', 0)

        a_s_iptms = [m['A_S_iptm_from_S'] for m in vd['models']
                     if m.get('A_S_iptm_from_S') is not None]

        entrance_models = [m for m in vd['models'] if m['binding_mode'] == 'entrance']
        mean_entrance_contacts = (
            np.mean([m['n_contacts_4a'] for m in entrance_models])
            if entrance_models else 0.0
        )

        l1_dists = [m['l1_to_active_site'] for m in vd['models']
                    if m.get('l1_to_active_site') is not None]

        stats_table[short] = {
            'n_models': n,
            'ki': vd['ki'],
            'titer': vd['titer'],
            'kcat_km': vd.get('kcat_km'),
            'mean_dist': float(np.mean(dists)),
            'median_dist': float(np.median(dists)),
            'std_dist': float(np.std(dists)),
            'min_dist': float(np.min(dists)),
            'max_dist': float(np.max(dists)),
            'n_close': n_close,
            'frac_close': n_close / n,
            'n_entrance': n_entrance,
            'n_active_site': n_active_site,
            'frac_entrance': n_entrance / n,
            'frac_active_site': n_active_site / n,
            'mean_entrance_contacts': float(mean_entrance_contacts),
            'mean_A_S_iptm': float(np.mean(a_s_iptms)) if a_s_iptms else None,
            'binding_modes': dict(modes),
            'mean_l1_to_active_site': float(np.mean(l1_dists)) if l1_dists else None,
        }
    return stats_table


def compute_correlations(stats_table):
    corrs = {}
    variants_with_ki = [v for v in VARIANT_ORDER
                        if v in stats_table and stats_table[v].get('ki') is not None]
    if len(variants_with_ki) >= 3:
        ki_vals = [stats_table[v]['ki'] for v in variants_with_ki]
        for metric, label in [('frac_close', 'Ki_vs_frac_close'),
                              ('frac_entrance', 'Ki_vs_frac_entrance'),
                              ('mean_dist', 'Ki_vs_mean_dist')]:
            vals = [stats_table[v][metric] for v in variants_with_ki]
            r, p = stats.pearsonr(ki_vals, vals)
            rho, ps = stats.spearmanr(ki_vals, vals)
            corrs[label] = {'r': float(r), 'p': float(p),
                            'rho': float(rho), 'p_spearman': float(ps),
                            'n': len(ki_vals), 'variants': variants_with_ki}

    variants_with_titer = [v for v in VARIANT_ORDER
                           if v in stats_table and stats_table[v].get('titer') is not None]
    if len(variants_with_titer) >= 3:
        titer_vals = [stats_table[v]['titer'] for v in variants_with_titer]
        for metric, label in [('frac_close', 'titer_vs_frac_close'),
                              ('frac_entrance', 'titer_vs_frac_entrance')]:
            vals = [stats_table[v][metric] for v in variants_with_titer]
            r, p = stats.pearsonr(titer_vals, vals)
            corrs[label] = {'r': float(r), 'p': float(p),
                            'n': len(variants_with_titer)}
    return corrs


def consensus_contacts(all_data):
    all_close_contacts = Counter()
    n_total_close = 0
    for short in VARIANT_ORDER:
        vd = all_data.get(short)
        if not vd or not vd['models']:
            continue
        close_models = [m for m in vd['models'] if m['l2_l1_dist'] < CLOSE_THRESHOLD]
        n_total_close += len(close_models)
        for m in close_models:
            for rn_str in m['contacts_4a']:
                all_close_contacts[int(rn_str)] += 1
    return dict(all_close_contacts.most_common(30)), n_total_close


def group_comparison(stats_table):
    high_titer = [v for v in VARIANT_ORDER
                  if v in stats_table and stats_table[v].get('titer') is not None
                  and stats_table[v]['titer'] >= 700]
    low_titer = [v for v in VARIANT_ORDER
                 if v in stats_table and stats_table[v].get('titer') is not None
                 and stats_table[v]['titer'] < 400]
    results = {
        'high_titer_variants': high_titer,
        'low_titer_variants': low_titer,
    }
    for metric in ['frac_entrance', 'frac_close', 'mean_dist']:
        high_vals = [stats_table[v][metric] for v in high_titer]
        low_vals = [stats_table[v][metric] for v in low_titer]
        if len(high_vals) >= 2 and len(low_vals) >= 2:
            u_stat, p_val = stats.mannwhitneyu(high_vals, low_vals, alternative='two-sided')
            results[metric] = {
                'high_mean': float(np.mean(high_vals)),
                'low_mean': float(np.mean(low_vals)),
                'U': float(u_stat), 'p': float(p_val),
            }
    return results


def run_analysis(condition):
    if condition == 'no_cofactor':
        base_dir = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap_no_cofactor'
        output_dir = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis/analysis_no_cofactor'
        s_chain_idx = 2  # A=0, L=1, S=2
    elif condition == 'no_constraint':
        base_dir = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap_no_constraint'
        output_dir = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis/analysis_no_constraint'
        s_chain_idx = 2  # A=0, L=1, S=2, T=3, M=4
    else:
        raise ValueError(f'Unknown condition: {condition}')

    print(f'\n=== Analyzing {condition} (S chain idx = {s_chain_idx}) ===')
    print(f'Base: {base_dir}')
    print(f'Output: {output_dir}')

    all_data = scan_all_models(base_dir, s_chain_idx)
    stats_table = analyze_distributions(all_data)
    corrs = compute_correlations(stats_table)
    contacts, n_close_total = consensus_contacts(all_data)
    groups = group_comparison(stats_table)

    results = {
        'condition': condition,
        'base_dir': base_dir,
        'total_variants': len(all_data),
        'total_models': sum(v['n_models'] for v in all_data.values()),
        'stats': stats_table,
        'correlations': corrs,
        'consensus_contacts': contacts,
        'n_close_models_total': n_close_total,
        'group_comparison': groups,
    }

    # Make serializable
    for v in stats_table.values():
        for k, val in v.items():
            if isinstance(val, (np.floating, np.integer)):
                v[k] = float(val) if isinstance(val, np.floating) else int(val)

    out_path = os.path.join(output_dir, 'analysis_results.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\nSaved results to {out_path}')
    print(f'Total models analyzed: {results["total_models"]}')
    print(f'Total close-binding models: {n_close_total}')

    return results


if __name__ == '__main__':
    condition = sys.argv[1] if len(sys.argv) > 1 else 'no_cofactor'
    run_analysis(condition)
