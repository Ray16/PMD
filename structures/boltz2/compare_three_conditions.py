#!/usr/bin/env python3
"""Compare Boltz2 cofolding results across three conditions:
1. Constrained (cofactors + constraints on cofactors/substrates)
2. No constraint (cofactors present but no constraints)
3. No cofactor (no ATP/Mg2+, no constraints)

Chains: constrained/no_constraint have A(protein), L(substrate), S(2nd substrate), M(Mg), T(ATP)
        no_cofactor has A(protein), L(substrate), S(2nd substrate)
"""
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

BASES = {
    'constrained': '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap',
    'no_constraint': '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap_no_constraint',
    'no_cofactor': '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap_no_cofactor',
}

OUTPUT_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2'

VARIANT_DIRS = {
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

VARIANT_ORDER = list(VARIANT_DIRS.keys())

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
    n_catalytic_near = sum(1 for d in l2_to_catalytic_dists.values() if d < 6.0)

    if l2_l1_dist > CLOSE_THRESHOLD:
        if l2_l1_dist > INTERMEDIATE_THRESHOLD:
            binding_mode = 'distant'
        else:
            binding_mode = 'intermediate'
    elif n_catalytic_contacted >= 5:
        binding_mode = 'active_site'
    else:
        binding_mode = 'entrance'

    # Also compute protein-substrate (L) and protein-S distances
    prot_com = np.array([a['xyz'] for a in chains['A']]).mean(axis=0)
    l_prot_dist = float(np.linalg.norm(l1_com - prot_com))
    s_prot_dist = float(np.linalg.norm(s_com - prot_com))

    return {
        'l2_l1_dist': l2_l1_dist,
        'closest_residues': closest_residues,
        'contacts_4a': {str(k): v for k, v in contacts_4a.items()},
        'n_contacts_4a': len(contacts_4a),
        'n_catalytic_contacted': n_catalytic_contacted,
        'n_catalytic_near': n_catalytic_near,
        'catalytic_dists': {str(k): round(v, 2) for k, v in l2_to_catalytic_dists.items()},
        'binding_mode': binding_mode,
    }


def find_all_pdbs(base, dirname):
    """Find PDBs from main dir and all seed dirs."""
    pdbs = []
    main_pred = os.path.join(base, dirname, f'boltz_results_{dirname}',
                              'predictions', dirname)
    if os.path.isdir(main_pred):
        pdbs.extend(sorted(glob.glob(os.path.join(main_pred, f'{dirname}_model_*.pdb'))))

    # Also search seed directories
    for seed_dir in sorted(glob.glob(os.path.join(base, f'{dirname}_seed*'))):
        seed_pred = os.path.join(seed_dir, f'boltz_results_{dirname}',
                                  'predictions', dirname)
        if os.path.isdir(seed_pred):
            pdbs.extend(sorted(glob.glob(os.path.join(seed_pred, f'{dirname}_model_*.pdb'))))
    return pdbs


def get_confidence(pdb_path, dirname):
    """Get confidence JSON corresponding to a PDB file."""
    pred_dir = os.path.dirname(pdb_path)
    model_name = os.path.basename(pdb_path).replace('.pdb', '')
    conf_path = os.path.join(pred_dir, f'confidence_{model_name}.json')
    if os.path.exists(conf_path):
        with open(conf_path) as f:
            return json.load(f)
    return None


def scan_condition(condition_name, base):
    all_data = {}
    for short, dirname in VARIANT_DIRS.items():
        pdbs = find_all_pdbs(base, dirname)
        models = []
        for pdb_path in pdbs:
            chains = parse_pdb_minimal(pdb_path)
            metrics = compute_l2_metrics(chains)
            if not metrics:
                continue

            conf_data = get_confidence(pdb_path, dirname)
            conf = {}
            if conf_data:
                # Chain indices differ: no_cofactor has 3 chains (0=A, 1=L, 2=S)
                # constrained/no_constraint have 5 chains (0=A, 1=L, 2=S, 3=T/M, 4=M/T)
                conf = {
                    'confidence_score': conf_data.get('confidence_score'),
                    'iptm': conf_data.get('iptm'),
                    'ligand_iptm': conf_data.get('ligand_iptm'),
                    'ptm': conf_data.get('ptm'),
                    'complex_plddt': conf_data.get('complex_plddt'),
                }

            model_name = os.path.basename(pdb_path).replace('.pdb', '')
            models.append({
                'model_id': model_name,
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
    return all_data


def summarize(all_data):
    stats_table = {}
    for short in VARIANT_ORDER:
        vd = all_data.get(short)
        if not vd or not vd['models']:
            continue
        dists = [m['l2_l1_dist'] for m in vd['models']]
        n = len(dists)
        modes = Counter(m['binding_mode'] for m in vd['models'])
        n_entrance = modes.get('entrance', 0)
        n_active_site = modes.get('active_site', 0)
        n_close = sum(1 for d in dists if d < CLOSE_THRESHOLD)

        conf_scores = [m['confidence_score'] for m in vd['models']
                       if m.get('confidence_score') is not None]
        iptms = [m['iptm'] for m in vd['models'] if m.get('iptm') is not None]
        plddts = [m['complex_plddt'] for m in vd['models']
                  if m.get('complex_plddt') is not None]

        stats_table[short] = {
            'n_models': n,
            'ki': vd['ki'],
            'titer': vd['titer'],
            'kcat_km': vd['kcat_km'],
            'mean_dist': float(np.mean(dists)),
            'std_dist': float(np.std(dists)),
            'median_dist': float(np.median(dists)),
            'frac_entrance': n_entrance / n if n > 0 else 0,
            'frac_active_site': n_active_site / n if n > 0 else 0,
            'frac_close': n_close / n if n > 0 else 0,
            'n_entrance': n_entrance,
            'n_active_site': n_active_site,
            'binding_modes': dict(modes),
            'mean_confidence': float(np.mean(conf_scores)) if conf_scores else None,
            'mean_iptm': float(np.mean(iptms)) if iptms else None,
            'mean_plddt': float(np.mean(plddts)) if plddts else None,
            'mean_contacts': float(np.mean([m['n_contacts_4a'] for m in vd['models']])),
        }
    return stats_table


def print_comparison_table(all_stats):
    print('\n' + '=' * 120)
    print(f'{"Variant":>22} | {"Condition":>15} | {"N":>3} | {"FracEnt":>7} | {"FracAS":>7} | '
          f'{"FracClose":>9} | {"MeanDist":>8} | {"StdDist":>7} | {"Conf":>6} | {"iPTM":>6} | {"pLDDT":>6}')
    print('-' * 120)

    for v in VARIANT_ORDER:
        first = True
        for cond in ['constrained', 'no_constraint', 'no_cofactor']:
            s = all_stats[cond].get(v)
            if not s:
                continue
            label = v if first else ''
            first = False
            conf_str = f'{s["mean_confidence"]:.3f}' if s['mean_confidence'] is not None else 'N/A'
            iptm_str = f'{s["mean_iptm"]:.3f}' if s['mean_iptm'] is not None else 'N/A'
            plddt_str = f'{s["mean_plddt"]:.3f}' if s['mean_plddt'] is not None else 'N/A'
            print(f'{label:>22} | {cond:>15} | {s["n_models"]:>3} | '
                  f'{s["frac_entrance"]:>7.2f} | {s["frac_active_site"]:>7.2f} | '
                  f'{s["frac_close"]:>9.2f} | {s["mean_dist"]:>8.1f} | {s["std_dist"]:>7.1f} | '
                  f'{conf_str:>6} | {iptm_str:>6} | {plddt_str:>6}')
        if not first:
            print('-' * 120)


def compute_condition_correlations(stats_table, condition_name):
    corrs = {}
    variants_with_ki = [v for v in VARIANT_ORDER
                        if v in stats_table and stats_table[v].get('ki') is not None]
    if len(variants_with_ki) >= 3:
        ki_vals = [stats_table[v]['ki'] for v in variants_with_ki]
        for metric, label in [('frac_entrance', 'Ki_vs_frac_entrance'),
                              ('frac_close', 'Ki_vs_frac_close'),
                              ('mean_dist', 'Ki_vs_mean_dist')]:
            vals = [stats_table[v][metric] for v in variants_with_ki]
            r, p = stats.pearsonr(ki_vals, vals)
            rho, ps = stats.spearmanr(ki_vals, vals)
            corrs[label] = {'r': r, 'p': p, 'rho': rho, 'p_spearman': ps,
                            'n': len(ki_vals), 'variants': variants_with_ki}

    variants_with_titer = [v for v in VARIANT_ORDER
                           if v in stats_table and stats_table[v].get('titer') is not None]
    if len(variants_with_titer) >= 3:
        titer_vals = [stats_table[v]['titer'] for v in variants_with_titer]
        for metric, label in [('frac_entrance', 'titer_vs_frac_entrance'),
                              ('frac_close', 'titer_vs_frac_close')]:
            vals = [stats_table[v][metric] for v in variants_with_titer]
            r, p = stats.pearsonr(titer_vals, vals)
            rho, ps = stats.spearmanr(titer_vals, vals)
            corrs[label] = {'r': r, 'p': p, 'rho': rho, 'p_spearman': ps,
                            'n': len(variants_with_titer)}

    return corrs


def plot_comparison(all_stats, all_corrs):
    conditions = ['constrained', 'no_constraint', 'no_cofactor']
    cond_colors = {'constrained': '#1f77b4', 'no_constraint': '#ff7f0e', 'no_cofactor': '#2ca02c'}
    cond_labels = {'constrained': 'Constrained', 'no_constraint': 'No Constraint', 'no_cofactor': 'No Cofactor'}

    fig = plt.figure(figsize=(24, 28))
    gs = gridspec.GridSpec(4, 2, hspace=0.45, wspace=0.35)

    # --- Panel A: Entrance fraction comparison across conditions ---
    ax = fig.add_subplot(gs[0, 0])
    available = [v for v in VARIANT_ORDER
                 if any(v in all_stats[c] for c in conditions)]
    x = np.arange(len(available))
    width = 0.25
    for i, cond in enumerate(conditions):
        vals = []
        for v in available:
            s = all_stats[cond].get(v)
            vals.append(s['frac_entrance'] if s else 0)
        ax.bar(x + (i - 1) * width, vals, width, label=cond_labels[cond],
               color=cond_colors[cond], edgecolor='k', linewidth=0.5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Entrance-Binding Fraction')
    ax.set_title('A. Entrance-Binding Fraction by Variant & Condition', fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)

    # --- Panel B: Active-site fraction comparison ---
    ax = fig.add_subplot(gs[0, 1])
    for i, cond in enumerate(conditions):
        vals = []
        for v in available:
            s = all_stats[cond].get(v)
            vals.append(s['frac_active_site'] if s else 0)
        ax.bar(x + (i - 1) * width, vals, width, label=cond_labels[cond],
               color=cond_colors[cond], edgecolor='k', linewidth=0.5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Active-Site Fraction')
    ax.set_title('B. Active-Site Binding Fraction by Variant & Condition', fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)

    # --- Panel C: Mean S-L distance comparison ---
    ax = fig.add_subplot(gs[1, 0])
    for i, cond in enumerate(conditions):
        vals = []
        errs = []
        for v in available:
            s = all_stats[cond].get(v)
            if s:
                vals.append(s['mean_dist'])
                errs.append(s['std_dist'])
            else:
                vals.append(0)
                errs.append(0)
        ax.bar(x + (i - 1) * width, vals, width, yerr=errs,
               label=cond_labels[cond], color=cond_colors[cond],
               edgecolor='k', linewidth=0.5, alpha=0.8, capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Mean S-L Distance (Angstrom)')
    ax.set_title('C. Mean Secondary Substrate Distance', fontweight='bold')
    ax.legend(fontsize=8)
    ax.axhline(y=CLOSE_THRESHOLD, color='r', linestyle='--', alpha=0.5, label='Close threshold')

    # --- Panel D: Confidence score comparison ---
    ax = fig.add_subplot(gs[1, 1])
    for i, cond in enumerate(conditions):
        vals = []
        for v in available:
            s = all_stats[cond].get(v)
            if s and s['mean_confidence'] is not None:
                vals.append(s['mean_confidence'])
            else:
                vals.append(0)
        ax.bar(x + (i - 1) * width, vals, width, label=cond_labels[cond],
               color=cond_colors[cond], edgecolor='k', linewidth=0.5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Mean Confidence Score')
    ax.set_title('D. Boltz2 Confidence Score Comparison', fontweight='bold')
    ax.legend(fontsize=8)

    # --- Panel E: Ki correlation comparison ---
    ax = fig.add_subplot(gs[2, 0])
    ki_variants = [v for v in VARIANT_ORDER
                   if any(v in all_stats[c] and all_stats[c][v].get('ki') is not None
                          for c in conditions)]
    markers = {'constrained': 'o', 'no_constraint': 's', 'no_cofactor': '^'}
    for cond in conditions:
        ki_v = [v for v in ki_variants if v in all_stats[cond] and all_stats[cond][v].get('ki') is not None]
        if not ki_v:
            continue
        ki_vals = [all_stats[cond][v]['ki'] for v in ki_v]
        fe_vals = [all_stats[cond][v]['frac_entrance'] for v in ki_v]
        ax.scatter(fe_vals, ki_vals, c=cond_colors[cond], s=120,
                   marker=markers[cond], edgecolors='k', linewidth=1,
                   label=cond_labels[cond], zorder=3)
        for v, fe, ki in zip(ki_v, fe_vals, ki_vals):
            ax.annotate(v, (fe, ki), textcoords="offset points",
                       xytext=(5, 3), fontsize=7)

        corr_info = all_corrs[cond].get('Ki_vs_frac_entrance', {})
        if corr_info:
            z = np.polyfit(fe_vals, ki_vals, 1)
            xline = np.linspace(min(fe_vals) - 0.05, max(fe_vals) + 0.05, 100)
            ax.plot(xline, np.polyval(z, xline), '--', color=cond_colors[cond], alpha=0.4)

    legend_text = []
    for cond in conditions:
        ci = all_corrs[cond].get('Ki_vs_frac_entrance', {})
        if ci:
            legend_text.append(f'{cond_labels[cond]}: r={ci["r"]:.2f} (p={ci["p"]:.3f})')
    ax.set_xlabel('Entrance-Binding Fraction')
    ax.set_ylabel('Ki (mM)')
    ax.set_title('E. Ki vs Entrance Fraction\n' + ', '.join(legend_text),
                 fontweight='bold', fontsize=9)
    ax.legend(fontsize=8)

    # --- Panel F: Titer correlation comparison ---
    ax = fig.add_subplot(gs[2, 1])
    titer_variants = [v for v in VARIANT_ORDER
                      if any(v in all_stats[c] and all_stats[c][v].get('titer') is not None
                             for c in conditions)]
    for cond in conditions:
        tv = [v for v in titer_variants if v in all_stats[cond] and all_stats[cond][v].get('titer') is not None]
        if not tv:
            continue
        titer_vals = [all_stats[cond][v]['titer'] for v in tv]
        fe_vals = [all_stats[cond][v]['frac_entrance'] for v in tv]
        ax.scatter(fe_vals, titer_vals, c=cond_colors[cond], s=80,
                   marker=markers[cond], edgecolors='k', linewidth=0.8,
                   label=cond_labels[cond], zorder=3, alpha=0.8)

    legend_text = []
    for cond in conditions:
        ci = all_corrs[cond].get('titer_vs_frac_entrance', {})
        if ci:
            legend_text.append(f'{cond_labels[cond]}: r={ci["r"]:.2f} (p={ci["p"]:.3f})')
    ax.set_xlabel('Entrance-Binding Fraction')
    ax.set_ylabel('Titer (mg/L)')
    ax.set_title('F. Titer vs Entrance Fraction\n' + ', '.join(legend_text),
                 fontweight='bold', fontsize=9)
    ax.legend(fontsize=8)

    # --- Panel G: Binding mode composition heatmap ---
    ax = fig.add_subplot(gs[3, 0])
    mode_names = ['entrance', 'active_site', 'intermediate', 'distant']
    mode_colors_list = ['#ff7f0e', '#d62728', '#ffbb78', '#aec7e8']
    data_matrix = []
    row_labels = []
    for v in available:
        for cond in conditions:
            s = all_stats[cond].get(v)
            if not s:
                continue
            n = s['n_models']
            row = [s['binding_modes'].get(m, 0) / n if n > 0 else 0 for m in mode_names]
            data_matrix.append(row)
            row_labels.append(f'{v} ({cond_labels[cond][:4]})')

    if data_matrix:
        data_matrix = np.array(data_matrix)
        im = ax.imshow(data_matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=1)
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=5)
        ax.set_xticks(range(len(mode_names)))
        ax.set_xticklabels([m.replace('_', ' ').title() for m in mode_names], fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.6, label='Fraction')
        ax.set_title('G. Binding Mode Heatmap', fontweight='bold')

    # --- Panel H: Summary table ---
    ax = fig.add_subplot(gs[3, 1])
    ax.axis('off')
    headers = ['Variant', 'Cond', 'N', 'Frac Ent', 'Frac AS', 'Mean Dist', 'Conf', 'iPTM']
    cell_text = []
    for v in VARIANT_ORDER:
        for cond in conditions:
            s = all_stats[cond].get(v)
            if not s:
                continue
            conf_str = f'{s["mean_confidence"]:.3f}' if s['mean_confidence'] is not None else 'N/A'
            iptm_str = f'{s["mean_iptm"]:.3f}' if s['mean_iptm'] is not None else 'N/A'
            cell_text.append([
                v, cond[:6], str(s['n_models']),
                f'{s["frac_entrance"]:.2f}', f'{s["frac_active_site"]:.2f}',
                f'{s["mean_dist"]:.1f}', conf_str, iptm_str,
            ])

    if cell_text:
        table = ax.table(cellText=cell_text, colLabels=headers,
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(5)
        table.scale(1.0, 1.1)
        for key, cell in table.get_celld().items():
            if key[0] == 0:
                cell.set_facecolor('#e6e6e6')
                cell.set_text_props(fontweight='bold')
        ax.set_title('H. Summary Table', fontweight='bold', fontsize=10)

    fig.suptitle('Boltz2 Cofolding: Three-Condition Comparison\n'
                 '(Constrained vs No Constraint vs No Cofactor)',
                 fontsize=14, fontweight='bold', y=1.01)

    path = os.path.join(OUTPUT_DIR, 'three_condition_comparison.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {path}')


def plot_pairwise_condition(all_stats):
    """Pairwise scatter of entrance fraction between conditions."""
    conditions = ['constrained', 'no_constraint', 'no_cofactor']
    cond_labels = {'constrained': 'Constrained', 'no_constraint': 'No Constraint', 'no_cofactor': 'No Cofactor'}
    pairs = [('constrained', 'no_constraint'), ('constrained', 'no_cofactor'),
             ('no_constraint', 'no_cofactor')]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for idx, (c1, c2) in enumerate(pairs):
        ax = axes[idx]
        common = [v for v in VARIANT_ORDER
                  if v in all_stats[c1] and v in all_stats[c2]
                  and all_stats[c1][v]['n_models'] > 0 and all_stats[c2][v]['n_models'] > 0]

        x_vals = [all_stats[c1][v]['frac_entrance'] for v in common]
        y_vals = [all_stats[c2][v]['frac_entrance'] for v in common]

        for v, xv, yv in zip(common, x_vals, y_vals):
            ki = VARIANT_DATA[v].get('ki')
            if ki is not None:
                color = '#d62728'
                marker = 'D'
            else:
                color = '#1f77b4'
                marker = 'o'
            ax.scatter(xv, yv, c=color, s=80, marker=marker, edgecolors='k',
                      linewidth=0.8, zorder=3)
            ax.annotate(v, (xv, yv), textcoords="offset points",
                       xytext=(4, 3), fontsize=7)

        lims = [-0.05, 1.05]
        ax.plot(lims, lims, 'k--', alpha=0.3)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel(f'{cond_labels[c1]} Entrance Frac')
        ax.set_ylabel(f'{cond_labels[c2]} Entrance Frac')

        if len(common) >= 3:
            r, p = stats.pearsonr(x_vals, y_vals)
            ax.set_title(f'{cond_labels[c1]} vs {cond_labels[c2]}\nr={r:.2f} (p={p:.3f}, n={len(common)})',
                         fontweight='bold', fontsize=9)
        else:
            ax.set_title(f'{cond_labels[c1]} vs {cond_labels[c2]}', fontweight='bold', fontsize=9)

    fig.suptitle('Pairwise Comparison: Entrance Fraction Across Conditions\n(diamonds = Ki measured)',
                 fontweight='bold', fontsize=12)
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'pairwise_entrance_comparison.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {path}')


def plot_delta_analysis(all_stats):
    """Analyze how much each metric changes between conditions."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    common = [v for v in VARIANT_ORDER
              if all(v in all_stats[c] and all_stats[c][v]['n_models'] > 0
                     for c in ['constrained', 'no_constraint', 'no_cofactor'])]

    metrics = [
        ('frac_entrance', 'Entrance Fraction'),
        ('frac_active_site', 'Active-Site Fraction'),
        ('mean_dist', 'Mean S-L Distance'),
        ('mean_confidence', 'Mean Confidence'),
    ]

    for idx, (metric, label) in enumerate(metrics):
        ax = axes[idx // 2][idx % 2]

        x_pos = np.arange(len(common))
        constrained_vals = []
        no_constraint_vals = []
        no_cofactor_vals = []

        for v in common:
            cv = all_stats['constrained'][v].get(metric, 0) or 0
            ncv = all_stats['no_constraint'][v].get(metric, 0) or 0
            nfv = all_stats['no_cofactor'][v].get(metric, 0) or 0
            constrained_vals.append(cv)
            no_constraint_vals.append(ncv)
            no_cofactor_vals.append(nfv)

        width = 0.25
        ax.bar(x_pos - width, constrained_vals, width, label='Constrained',
               color='#1f77b4', edgecolor='k', linewidth=0.5, alpha=0.8)
        ax.bar(x_pos, no_constraint_vals, width, label='No Constraint',
               color='#ff7f0e', edgecolor='k', linewidth=0.5, alpha=0.8)
        ax.bar(x_pos + width, no_cofactor_vals, width, label='No Cofactor',
               color='#2ca02c', edgecolor='k', linewidth=0.5, alpha=0.8)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(common, rotation=55, ha='right', fontsize=7)
        ax.set_ylabel(label)
        ax.set_title(label, fontweight='bold')
        ax.legend(fontsize=7)

    fig.suptitle('Metric Comparison Across Three Conditions',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'delta_analysis.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {path}')


def main():
    print('=' * 80)
    print('  Three-Condition Boltz2 Comparison')
    print('  1. Constrained (cofactors + constraints)')
    print('  2. No constraint (cofactors, no constraints)')
    print('  3. No cofactor (no ATP/Mg2+)')
    print('=' * 80)

    all_data = {}
    all_stats = {}
    all_corrs = {}

    for cond, base in BASES.items():
        print(f'\n--- Scanning {cond} ---')
        all_data[cond] = scan_condition(cond, base)
        all_stats[cond] = summarize(all_data[cond])
        all_corrs[cond] = compute_condition_correlations(all_stats[cond], cond)

        n_total = sum(all_stats[cond][v]['n_models'] for v in VARIANT_ORDER if v in all_stats[cond])
        n_variants = sum(1 for v in VARIANT_ORDER if v in all_stats[cond] and all_stats[cond][v]['n_models'] > 0)
        print(f'  {n_variants} variants, {n_total} total models')

    print_comparison_table(all_stats)

    print('\n' + '=' * 80)
    print('  CORRELATION COMPARISON')
    print('=' * 80)
    for corr_key in ['Ki_vs_frac_entrance', 'Ki_vs_frac_close', 'Ki_vs_mean_dist',
                     'titer_vs_frac_entrance', 'titer_vs_frac_close']:
        print(f'\n  {corr_key}:')
        for cond in ['constrained', 'no_constraint', 'no_cofactor']:
            ci = all_corrs[cond].get(corr_key)
            if ci:
                print(f'    {cond:>15}: r={ci["r"]:.3f} (p={ci["p"]:.4f}), '
                      f'rho={ci["rho"]:.3f} (p={ci["p_spearman"]:.4f}), n={ci["n"]}')
            else:
                print(f'    {cond:>15}: N/A')

    print('\n' + '=' * 80)
    print('  CONDITION EFFECT SUMMARY')
    print('=' * 80)

    for metric, label in [('frac_entrance', 'Entrance Fraction'),
                          ('frac_active_site', 'Active-Site Fraction'),
                          ('mean_dist', 'Mean S-L Distance'),
                          ('mean_confidence', 'Mean Confidence')]:
        print(f'\n  {label}:')
        for cond in ['constrained', 'no_constraint', 'no_cofactor']:
            vals = [all_stats[cond][v][metric] for v in VARIANT_ORDER
                    if v in all_stats[cond] and all_stats[cond][v].get(metric) is not None]
            if vals:
                print(f'    {cond:>15}: mean={np.mean(vals):.3f}, '
                      f'std={np.std(vals):.3f}, '
                      f'range=[{np.min(vals):.3f}, {np.max(vals):.3f}]')

    # Pairwise condition changes
    print('\n  Per-variant entrance fraction changes:')
    common = [v for v in VARIANT_ORDER
              if all(v in all_stats[c] and all_stats[c][v]['n_models'] > 0
                     for c in ['constrained', 'no_constraint', 'no_cofactor'])]
    print(f'  (Variants in all three conditions: {len(common)})')
    if common:
        c_vals = np.array([all_stats['constrained'][v]['frac_entrance'] for v in common])
        nc_vals = np.array([all_stats['no_constraint'][v]['frac_entrance'] for v in common])
        nf_vals = np.array([all_stats['no_cofactor'][v]['frac_entrance'] for v in common])

        delta_nc = nc_vals - c_vals
        delta_nf = nf_vals - c_vals
        delta_nc_nf = nf_vals - nc_vals

        print(f'    Constrained -> No Constraint: mean delta = {np.mean(delta_nc):+.3f} '
              f'(std={np.std(delta_nc):.3f})')
        print(f'    Constrained -> No Cofactor:   mean delta = {np.mean(delta_nf):+.3f} '
              f'(std={np.std(delta_nf):.3f})')
        print(f'    No Constraint -> No Cofactor: mean delta = {np.mean(delta_nc_nf):+.3f} '
              f'(std={np.std(delta_nc_nf):.3f})')

        if len(common) >= 3:
            t_nc, p_nc = stats.ttest_rel(c_vals, nc_vals)
            t_nf, p_nf = stats.ttest_rel(c_vals, nf_vals)
            print(f'    Paired t-test (constrained vs no_constraint): t={t_nc:.2f}, p={p_nc:.4f}')
            print(f'    Paired t-test (constrained vs no_cofactor):   t={t_nf:.2f}, p={p_nf:.4f}')

    print('\n--- Generating figures ---')
    plot_comparison(all_stats, all_corrs)
    plot_pairwise_condition(all_stats)
    plot_delta_analysis(all_stats)

    # Save JSON
    output = {
        'conditions': {},
        'correlations': {},
    }
    for cond in ['constrained', 'no_constraint', 'no_cofactor']:
        output['conditions'][cond] = {}
        for v, s in all_stats[cond].items():
            output['conditions'][cond][v] = {
                k: (float(val) if isinstance(val, (np.floating, np.integer)) else val)
                for k, val in s.items()
            }
        output['correlations'][cond] = {}
        for k, ci in all_corrs[cond].items():
            output['correlations'][cond][k] = {
                kk: (float(vv) if isinstance(vv, (np.floating, float)) else vv)
                for kk, vv in ci.items()
            }

    json_path = os.path.join(OUTPUT_DIR, 'three_condition_comparison.json')
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f'Saved {json_path}')

    print('\n--- DONE ---')


if __name__ == '__main__':
    main()
