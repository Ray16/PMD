#!/usr/bin/env python3
"""Compare Boltz2 cofolding results across four conditions:
1. Constrained (cofactors + constraints on cofactors/substrates)
2. No constraint (cofactors present but no constraints)
3. No cofactor (no ATP/Mg2+, no constraints)
4. No ATP (Mg2+ present, no ATP, L1 pocket-constrained, Mg-D302 constrained)

Chains: constrained/no_constraint have A(protein), L(substrate), S(2nd substrate), M(Mg), T(ATP)
        no_cofactor has A(protein), L(substrate), S(2nd substrate)
        no_atp has A(protein), L(substrate), S(2nd substrate), M(Mg)
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
    'no_atp': '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap_no_atp',
}

OUTPUT_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'

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

CONDITIONS = ['constrained', 'no_constraint', 'no_cofactor', 'no_atp']
COND_COLORS = {
    'constrained': '#1f77b4',
    'no_constraint': '#ff7f0e',
    'no_cofactor': '#2ca02c',
    'no_atp': '#d62728',
}
COND_LABELS = {
    'constrained': 'Constrained',
    'no_constraint': 'No Constraint',
    'no_cofactor': 'No Cofactor',
    'no_atp': 'No ATP (Mg only)',
}
COND_MARKERS = {
    'constrained': 'o',
    'no_constraint': 's',
    'no_cofactor': '^',
    'no_atp': 'D',
}


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

    # Mg2+ distance to S (for no_atp condition)
    mg_s_dist = None
    if 'M' in chains:
        mg_coords = np.array([a['xyz'] for a in chains['M']])
        mg_com = mg_coords.mean(axis=0)
        mg_s_dist = float(np.linalg.norm(mg_com - s_com))

    # L1 distance to active site center
    active_site_cas = []
    for a in chains['A']:
        if a['atom'] == 'CA' and a['resnum'] in {19, 22, 120, 121, 153}:
            active_site_cas.append(a['xyz'])
    l1_to_active_site = None
    if active_site_cas:
        as_center = np.mean(active_site_cas, axis=0)
        l1_to_active_site = float(np.linalg.norm(l1_com - as_center))

    return {
        'l2_l1_dist': l2_l1_dist,
        'closest_residues': closest_residues,
        'contacts_4a': {str(k): v for k, v in contacts_4a.items()},
        'n_contacts_4a': len(contacts_4a),
        'n_catalytic_contacted': n_catalytic_contacted,
        'n_catalytic_near': n_catalytic_near,
        'catalytic_dists': {str(k): round(v, 2) for k, v in l2_to_catalytic_dists.items()},
        'binding_mode': binding_mode,
        'mg_s_dist': mg_s_dist,
        'l1_to_active_site': l1_to_active_site,
    }


def find_all_pdbs(base, dirname):
    pdbs = []
    main_pred = os.path.join(base, dirname, f'boltz_results_{dirname}',
                              'predictions', dirname)
    if os.path.isdir(main_pred):
        pdbs.extend(sorted(glob.glob(os.path.join(main_pred, f'{dirname}_model_*.pdb'))))
    for seed_dir in sorted(glob.glob(os.path.join(base, f'{dirname}_seed*'))):
        seed_pred = os.path.join(seed_dir, f'boltz_results_{dirname}',
                                  'predictions', dirname)
        if os.path.isdir(seed_pred):
            pdbs.extend(sorted(glob.glob(os.path.join(seed_pred, f'{dirname}_model_*.pdb'))))
    return pdbs


def get_confidence(pdb_path):
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

            conf_data = get_confidence(pdb_path)
            conf = {}
            if conf_data:
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

        mg_dists = [m['mg_s_dist'] for m in vd['models']
                    if m.get('mg_s_dist') is not None]

        entrance_models = [m for m in vd['models'] if m['binding_mode'] == 'entrance']
        mean_entrance_contacts = (
            float(np.mean([m['n_contacts_4a'] for m in entrance_models]))
            if entrance_models else 0.0
        )

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
            'mean_entrance_contacts': mean_entrance_contacts,
            'mean_mg_s_dist': float(np.mean(mg_dists)) if mg_dists else None,
        }
    return stats_table


def compute_condition_correlations(stats_table):
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


def print_comparison_table(all_stats):
    print('\n' + '=' * 140)
    print(f'{"Variant":>22} | {"Condition":>18} | {"N":>3} | {"FracEnt":>7} | {"FracAS":>7} | '
          f'{"FracClose":>9} | {"MeanDist":>8} | {"StdDist":>7} | {"Conf":>6} | {"iPTM":>6} | {"pLDDT":>6} | {"MgSDist":>7}')
    print('-' * 140)

    for v in VARIANT_ORDER:
        first = True
        for cond in CONDITIONS:
            s = all_stats[cond].get(v)
            if not s:
                continue
            label = v if first else ''
            first = False
            conf_str = f'{s["mean_confidence"]:.3f}' if s['mean_confidence'] is not None else 'N/A'
            iptm_str = f'{s["mean_iptm"]:.3f}' if s['mean_iptm'] is not None else 'N/A'
            plddt_str = f'{s["mean_plddt"]:.3f}' if s['mean_plddt'] is not None else 'N/A'
            mg_str = f'{s["mean_mg_s_dist"]:.1f}' if s.get('mean_mg_s_dist') is not None else 'N/A'
            print(f'{label:>22} | {cond:>18} | {s["n_models"]:>3} | '
                  f'{s["frac_entrance"]:>7.2f} | {s["frac_active_site"]:>7.2f} | '
                  f'{s["frac_close"]:>9.2f} | {s["mean_dist"]:>8.1f} | {s["std_dist"]:>7.1f} | '
                  f'{conf_str:>6} | {iptm_str:>6} | {plddt_str:>6} | {mg_str:>7}')
        if not first:
            print('-' * 140)


def plot_four_condition_comparison(all_stats, all_corrs):
    fig = plt.figure(figsize=(28, 36))
    gs = gridspec.GridSpec(5, 2, hspace=0.45, wspace=0.35)

    available = [v for v in VARIANT_ORDER
                 if any(v in all_stats[c] for c in CONDITIONS)]
    x = np.arange(len(available))
    width = 0.20

    # --- Panel A: Entrance fraction ---
    ax = fig.add_subplot(gs[0, 0])
    for i, cond in enumerate(CONDITIONS):
        vals = [all_stats[cond].get(v, {}).get('frac_entrance', 0) for v in available]
        ax.bar(x + (i - 1.5) * width, vals, width, label=COND_LABELS[cond],
               color=COND_COLORS[cond], edgecolor='k', linewidth=0.5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Entrance-Binding Fraction')
    ax.set_title('A. Entrance-Binding Fraction by Variant & Condition', fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_ylim(0, 1.05)

    # --- Panel B: Active-site fraction ---
    ax = fig.add_subplot(gs[0, 1])
    for i, cond in enumerate(CONDITIONS):
        vals = [all_stats[cond].get(v, {}).get('frac_active_site', 0) for v in available]
        ax.bar(x + (i - 1.5) * width, vals, width, label=COND_LABELS[cond],
               color=COND_COLORS[cond], edgecolor='k', linewidth=0.5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Active-Site Fraction')
    ax.set_title('B. Active-Site Binding Fraction by Variant & Condition', fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_ylim(0, 1.05)

    # --- Panel C: Mean S-L distance ---
    ax = fig.add_subplot(gs[1, 0])
    for i, cond in enumerate(CONDITIONS):
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
        ax.bar(x + (i - 1.5) * width, vals, width, yerr=errs,
               label=COND_LABELS[cond], color=COND_COLORS[cond],
               edgecolor='k', linewidth=0.5, alpha=0.8, capsize=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Mean S-L Distance (Angstrom)')
    ax.set_title('C. Mean Secondary Substrate Distance', fontweight='bold')
    ax.legend(fontsize=7)
    ax.axhline(y=CLOSE_THRESHOLD, color='gray', linestyle='--', alpha=0.5)

    # --- Panel D: Confidence score ---
    ax = fig.add_subplot(gs[1, 1])
    for i, cond in enumerate(CONDITIONS):
        vals = []
        for v in available:
            s = all_stats[cond].get(v)
            if s and s['mean_confidence'] is not None:
                vals.append(s['mean_confidence'])
            else:
                vals.append(0)
        ax.bar(x + (i - 1.5) * width, vals, width, label=COND_LABELS[cond],
               color=COND_COLORS[cond], edgecolor='k', linewidth=0.5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=55, ha='right', fontsize=7)
    ax.set_ylabel('Mean Confidence Score')
    ax.set_title('D. Boltz2 Confidence Score Comparison', fontweight='bold')
    ax.legend(fontsize=7)

    # --- Panel E: Ki correlation ---
    ax = fig.add_subplot(gs[2, 0])
    for cond in CONDITIONS:
        ki_v = [v for v in VARIANT_ORDER
                if v in all_stats[cond] and all_stats[cond][v].get('ki') is not None]
        if not ki_v:
            continue
        ki_vals = [all_stats[cond][v]['ki'] for v in ki_v]
        fe_vals = [all_stats[cond][v]['frac_entrance'] for v in ki_v]
        ax.scatter(fe_vals, ki_vals, c=COND_COLORS[cond], s=100,
                   marker=COND_MARKERS[cond], edgecolors='k', linewidth=0.8,
                   label=COND_LABELS[cond], zorder=3)
        for v, fe, ki in zip(ki_v, fe_vals, ki_vals):
            ax.annotate(v, (fe, ki), textcoords="offset points",
                       xytext=(5, 3), fontsize=6)
        if len(ki_v) >= 3:
            z = np.polyfit(fe_vals, ki_vals, 1)
            xline = np.linspace(min(fe_vals) - 0.05, max(fe_vals) + 0.05, 100)
            ax.plot(xline, np.polyval(z, xline), '--', color=COND_COLORS[cond], alpha=0.3)

    legend_text = []
    for cond in CONDITIONS:
        ci = all_corrs[cond].get('Ki_vs_frac_entrance', {})
        if ci:
            legend_text.append(f'{COND_LABELS[cond]}: r={ci["r"]:.2f} (p={ci["p"]:.3f})')
    ax.set_xlabel('Entrance-Binding Fraction')
    ax.set_ylabel('Ki (mM)')
    title_extra = '\n'.join(legend_text) if legend_text else ''
    ax.set_title(f'E. Ki vs Entrance Fraction\n{title_extra}',
                 fontweight='bold', fontsize=8)
    ax.legend(fontsize=7)

    # --- Panel F: Titer correlation ---
    ax = fig.add_subplot(gs[2, 1])
    for cond in CONDITIONS:
        tv = [v for v in VARIANT_ORDER
              if v in all_stats[cond] and all_stats[cond][v].get('titer') is not None]
        if not tv:
            continue
        titer_vals = [all_stats[cond][v]['titer'] for v in tv]
        fe_vals = [all_stats[cond][v]['frac_entrance'] for v in tv]
        ax.scatter(fe_vals, titer_vals, c=COND_COLORS[cond], s=70,
                   marker=COND_MARKERS[cond], edgecolors='k', linewidth=0.8,
                   label=COND_LABELS[cond], zorder=3, alpha=0.8)

    legend_text = []
    for cond in CONDITIONS:
        ci = all_corrs[cond].get('titer_vs_frac_entrance', {})
        if ci:
            legend_text.append(f'{COND_LABELS[cond]}: r={ci["r"]:.2f} (p={ci["p"]:.3f})')
    ax.set_xlabel('Entrance-Binding Fraction')
    ax.set_ylabel('Titer (mg/L)')
    title_extra = '\n'.join(legend_text) if legend_text else ''
    ax.set_title(f'F. Titer vs Entrance Fraction\n{title_extra}',
                 fontweight='bold', fontsize=8)
    ax.legend(fontsize=7)

    # --- Panel G: no_atp vs constrained entrance scatter ---
    ax = fig.add_subplot(gs[3, 0])
    common = [v for v in VARIANT_ORDER
              if v in all_stats['constrained'] and v in all_stats['no_atp']
              and all_stats['constrained'][v]['n_models'] > 0
              and all_stats['no_atp'][v]['n_models'] > 0]
    x_vals = [all_stats['constrained'][v]['frac_entrance'] for v in common]
    y_vals = [all_stats['no_atp'][v]['frac_entrance'] for v in common]
    for v, xv, yv in zip(common, x_vals, y_vals):
        ki = VARIANT_DATA[v].get('ki')
        color = '#d62728' if ki is not None else '#1f77b4'
        marker = 'D' if ki is not None else 'o'
        ax.scatter(xv, yv, c=color, s=80, marker=marker, edgecolors='k',
                  linewidth=0.8, zorder=3)
        ax.annotate(v, (xv, yv), textcoords="offset points",
                   xytext=(4, 3), fontsize=7)
    ax.plot([-0.05, 1.05], [-0.05, 1.05], 'k--', alpha=0.3)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('Constrained: Entrance Fraction')
    ax.set_ylabel('No ATP (Mg only): Entrance Fraction')
    if len(common) >= 3:
        r, p = stats.pearsonr(x_vals, y_vals)
        ax.set_title(f'G. Constrained vs No ATP\nr={r:.2f} (p={p:.3f}, n={len(common)})',
                     fontweight='bold', fontsize=9)
    else:
        ax.set_title('G. Constrained vs No ATP', fontweight='bold')

    # --- Panel H: no_atp vs no_cofactor entrance scatter ---
    ax = fig.add_subplot(gs[3, 1])
    common2 = [v for v in VARIANT_ORDER
               if v in all_stats['no_cofactor'] and v in all_stats['no_atp']
               and all_stats['no_cofactor'][v]['n_models'] > 0
               and all_stats['no_atp'][v]['n_models'] > 0]
    x_vals2 = [all_stats['no_cofactor'][v]['frac_entrance'] for v in common2]
    y_vals2 = [all_stats['no_atp'][v]['frac_entrance'] for v in common2]
    for v, xv, yv in zip(common2, x_vals2, y_vals2):
        ki = VARIANT_DATA[v].get('ki')
        color = '#d62728' if ki is not None else '#1f77b4'
        marker = 'D' if ki is not None else 'o'
        ax.scatter(xv, yv, c=color, s=80, marker=marker, edgecolors='k',
                  linewidth=0.8, zorder=3)
        ax.annotate(v, (xv, yv), textcoords="offset points",
                   xytext=(4, 3), fontsize=7)
    ax.plot([-0.05, 1.05], [-0.05, 1.05], 'k--', alpha=0.3)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('No Cofactor: Entrance Fraction')
    ax.set_ylabel('No ATP (Mg only): Entrance Fraction')
    if len(common2) >= 3:
        r, p = stats.pearsonr(x_vals2, y_vals2)
        ax.set_title(f'H. No Cofactor vs No ATP\nr={r:.2f} (p={p:.3f}, n={len(common2)})',
                     fontweight='bold', fontsize=9)
    else:
        ax.set_title('H. No Cofactor vs No ATP', fontweight='bold')

    # --- Panel I: Binding mode heatmap ---
    ax = fig.add_subplot(gs[4, 0])
    mode_names = ['entrance', 'active_site', 'intermediate', 'distant']
    data_matrix = []
    row_labels = []
    for v in available:
        for cond in CONDITIONS:
            s = all_stats[cond].get(v)
            if not s:
                continue
            n = s['n_models']
            row = [s['binding_modes'].get(m, 0) / n if n > 0 else 0 for m in mode_names]
            data_matrix.append(row)
            row_labels.append(f'{v} ({COND_LABELS[cond][:5]})')

    if data_matrix:
        data_matrix = np.array(data_matrix)
        im = ax.imshow(data_matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=1)
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=4)
        ax.set_xticks(range(len(mode_names)))
        ax.set_xticklabels([m.replace('_', ' ').title() for m in mode_names], fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.5, label='Fraction')
        ax.set_title('I. Binding Mode Heatmap (4 conditions)', fontweight='bold')

    # --- Panel J: Mg-S distance (no_atp specific) ---
    ax = fig.add_subplot(gs[4, 1])
    mg_variants = [v for v in VARIANT_ORDER
                   if v in all_stats['no_atp'] and all_stats['no_atp'][v].get('mean_mg_s_dist') is not None]
    if mg_variants:
        mg_dists = [all_stats['no_atp'][v]['mean_mg_s_dist'] for v in mg_variants]
        fe_vals = [all_stats['no_atp'][v]['frac_entrance'] for v in mg_variants]
        colors = []
        for v in mg_variants:
            ki = VARIANT_DATA[v].get('ki')
            colors.append('#d62728' if ki is not None else '#1f77b4')
        ax.scatter(mg_dists, fe_vals, c=colors, s=80, edgecolors='k', linewidth=0.8, zorder=3)
        for v, md, fe in zip(mg_variants, mg_dists, fe_vals):
            ax.annotate(v, (md, fe), textcoords="offset points", xytext=(4, 3), fontsize=7)
        ax.set_xlabel('Mean Mg-S Distance (Angstrom)')
        ax.set_ylabel('Entrance Fraction')
        if len(mg_variants) >= 3:
            r, p = stats.pearsonr(mg_dists, fe_vals)
            ax.set_title(f'J. Mg-S Distance vs Entrance Fraction (No ATP)\nr={r:.2f} (p={p:.3f})',
                         fontweight='bold', fontsize=9)
        else:
            ax.set_title('J. Mg-S Distance vs Entrance Fraction (No ATP)', fontweight='bold')
    else:
        ax.set_title('J. No Mg-S data available', fontweight='bold')
        ax.axis('off')

    fig.suptitle('Boltz2 Cofolding: Four-Condition Comparison\n'
                 '(Constrained vs No Constraint vs No Cofactor vs No ATP)',
                 fontsize=14, fontweight='bold', y=1.01)

    path = os.path.join(OUTPUT_DIR, 'four_condition_comparison.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {path}')


def plot_no_atp_focus(all_stats, all_data):
    """Focused analysis on no_atp condition: what does Mg2+ alone do?"""
    fig, axes = plt.subplots(2, 3, figsize=(24, 14))

    # Panel 1: no_atp entrance fraction vs constrained delta
    ax = axes[0, 0]
    common = [v for v in VARIANT_ORDER
              if all(v in all_stats[c] and all_stats[c][v]['n_models'] > 0
                     for c in ['constrained', 'no_atp'])]
    if common:
        c_vals = [all_stats['constrained'][v]['frac_entrance'] for v in common]
        na_vals = [all_stats['no_atp'][v]['frac_entrance'] for v in common]
        delta = [na - c for na, c in zip(na_vals, c_vals)]
        colors = ['#d62728' if VARIANT_DATA[v].get('ki') is not None else '#1f77b4' for v in common]
        x_pos = np.arange(len(common))
        ax.bar(x_pos, delta, color=colors, edgecolor='k', linewidth=0.5, alpha=0.8)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(common, rotation=55, ha='right', fontsize=7)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_ylabel('Delta Entrance Fraction\n(No ATP - Constrained)')
        ax.set_title('Effect of Removing ATP\n(positive = more entrance binding)', fontweight='bold', fontsize=9)

    # Panel 2: no_atp vs no_cofactor delta (isolating Mg2+ effect)
    ax = axes[0, 1]
    common2 = [v for v in VARIANT_ORDER
               if all(v in all_stats[c] and all_stats[c][v]['n_models'] > 0
                      for c in ['no_cofactor', 'no_atp'])]
    if common2:
        nf_vals = [all_stats['no_cofactor'][v]['frac_entrance'] for v in common2]
        na_vals = [all_stats['no_atp'][v]['frac_entrance'] for v in common2]
        delta = [na - nf for na, nf in zip(na_vals, nf_vals)]
        colors = ['#d62728' if VARIANT_DATA[v].get('ki') is not None else '#1f77b4' for v in common2]
        x_pos = np.arange(len(common2))
        ax.bar(x_pos, delta, color=colors, edgecolor='k', linewidth=0.5, alpha=0.8)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(common2, rotation=55, ha='right', fontsize=7)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_ylabel('Delta Entrance Fraction\n(No ATP - No Cofactor)')
        ax.set_title('Mg2+ Effect (No ATP vs No Cofactor)\n(positive = Mg increases entrance binding)',
                     fontweight='bold', fontsize=9)

    # Panel 3: Condition ranking per variant
    ax = axes[0, 2]
    all_four = [v for v in VARIANT_ORDER
                if all(v in all_stats[c] and all_stats[c][v]['n_models'] > 0 for c in CONDITIONS)]
    if all_four:
        for i, cond in enumerate(CONDITIONS):
            vals = [all_stats[cond][v]['frac_entrance'] for v in all_four]
            ax.plot(range(len(all_four)), vals, '-o', color=COND_COLORS[cond],
                    label=COND_LABELS[cond], markersize=5, alpha=0.8)
        ax.set_xticks(range(len(all_four)))
        ax.set_xticklabels(all_four, rotation=55, ha='right', fontsize=7)
        ax.set_ylabel('Entrance Fraction')
        ax.set_title('Entrance Fraction Across All 4 Conditions', fontweight='bold', fontsize=9)
        ax.legend(fontsize=7)
        ax.set_ylim(-0.05, 1.05)

    # Panel 4: S-L distance distributions per condition (box plots)
    ax = axes[1, 0]
    box_data = {cond: [] for cond in CONDITIONS}
    for cond in CONDITIONS:
        for v in VARIANT_ORDER:
            vd = all_data.get(cond, {}).get(v)
            if vd and vd['models']:
                for m in vd['models']:
                    box_data[cond].append(m['l2_l1_dist'])
    bp_data = [box_data[c] for c in CONDITIONS if box_data[c]]
    bp_labels = [COND_LABELS[c] for c in CONDITIONS if box_data[c]]
    bp_colors = [COND_COLORS[c] for c in CONDITIONS if box_data[c]]
    if bp_data:
        bp = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], bp_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax.set_ylabel('S-L Distance (Angstrom)')
        ax.set_title('S-L Distance Distribution by Condition', fontweight='bold', fontsize=9)
        ax.axhline(y=CLOSE_THRESHOLD, color='gray', linestyle='--', alpha=0.5)

    # Panel 5: Confidence comparison boxplot
    ax = axes[1, 1]
    conf_data = {cond: [] for cond in CONDITIONS}
    for cond in CONDITIONS:
        for v in VARIANT_ORDER:
            vd = all_data.get(cond, {}).get(v)
            if vd and vd['models']:
                for m in vd['models']:
                    if m.get('confidence_score') is not None:
                        conf_data[cond].append(m['confidence_score'])
    bp_data = [conf_data[c] for c in CONDITIONS if conf_data[c]]
    bp_labels = [COND_LABELS[c] for c in CONDITIONS if conf_data[c]]
    bp_colors = [COND_COLORS[c] for c in CONDITIONS if conf_data[c]]
    if bp_data:
        bp = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], bp_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax.set_ylabel('Confidence Score')
        ax.set_title('Confidence Score by Condition', fontweight='bold', fontsize=9)

    # Panel 6: Per-variant Mg-S distance (no_atp only)
    ax = axes[1, 2]
    na_data = all_data.get('no_atp', {})
    mg_variants = []
    mg_dist_lists = []
    for v in VARIANT_ORDER:
        vd = na_data.get(v)
        if vd and vd['models']:
            dists = [m['mg_s_dist'] for m in vd['models'] if m.get('mg_s_dist') is not None]
            if dists:
                mg_variants.append(v)
                mg_dist_lists.append(dists)
    if mg_dist_lists:
        bp = ax.boxplot(mg_dist_lists, labels=mg_variants, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('#d62728')
            patch.set_alpha(0.5)
        ax.set_xticklabels(mg_variants, rotation=55, ha='right', fontsize=7)
        ax.set_ylabel('Mg-S Distance (Angstrom)')
        ax.set_title('Mg2+ to S Distance per Variant (No ATP)', fontweight='bold', fontsize=9)

    fig.suptitle('No ATP Condition: Focused Analysis\n'
                 'Isolating Mg2+ contribution to competitive binding',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'no_atp_focused_analysis.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {path}')


def plot_pairwise_all(all_stats):
    """Pairwise scatter of entrance fraction between all condition pairs."""
    pairs = []
    for i, c1 in enumerate(CONDITIONS):
        for c2 in CONDITIONS[i+1:]:
            pairs.append((c1, c2))

    n_pairs = len(pairs)
    cols = 3
    rows = (n_pairs + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(7 * cols, 6 * rows))
    axes = axes.flatten() if n_pairs > 1 else [axes]

    for idx, (c1, c2) in enumerate(pairs):
        ax = axes[idx]
        common = [v for v in VARIANT_ORDER
                  if v in all_stats[c1] and v in all_stats[c2]
                  and all_stats[c1][v]['n_models'] > 0 and all_stats[c2][v]['n_models'] > 0]

        x_vals = [all_stats[c1][v]['frac_entrance'] for v in common]
        y_vals = [all_stats[c2][v]['frac_entrance'] for v in common]

        for v, xv, yv in zip(common, x_vals, y_vals):
            ki = VARIANT_DATA[v].get('ki')
            color = '#d62728' if ki is not None else '#1f77b4'
            marker = 'D' if ki is not None else 'o'
            ax.scatter(xv, yv, c=color, s=80, marker=marker, edgecolors='k',
                      linewidth=0.8, zorder=3)
            ax.annotate(v, (xv, yv), textcoords="offset points",
                       xytext=(4, 3), fontsize=6)

        ax.plot([-0.05, 1.05], [-0.05, 1.05], 'k--', alpha=0.3)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel(f'{COND_LABELS[c1]} Entrance Frac')
        ax.set_ylabel(f'{COND_LABELS[c2]} Entrance Frac')

        if len(common) >= 3:
            r, p = stats.pearsonr(x_vals, y_vals)
            ax.set_title(f'{COND_LABELS[c1]} vs {COND_LABELS[c2]}\nr={r:.2f} (p={p:.3f}, n={len(common)})',
                         fontweight='bold', fontsize=9)
        else:
            ax.set_title(f'{COND_LABELS[c1]} vs {COND_LABELS[c2]}', fontweight='bold', fontsize=9)

    for idx in range(len(pairs), len(axes)):
        axes[idx].axis('off')

    fig.suptitle('Pairwise Comparison: Entrance Fraction Across All Conditions\n(diamonds = Ki measured)',
                 fontweight='bold', fontsize=12)
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'pairwise_entrance_four_conditions.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {path}')


def main():
    print('=' * 80)
    print('  Four-Condition Boltz2 Comparison')
    print('  1. Constrained (cofactors + constraints)')
    print('  2. No constraint (cofactors, no constraints)')
    print('  3. No cofactor (no ATP/Mg2+)')
    print('  4. No ATP (Mg2+ only, L1 constrained, Mg-D302 constrained)')
    print('=' * 80)

    all_data = {}
    all_stats = {}
    all_corrs = {}

    for cond, base in BASES.items():
        print(f'\n--- Scanning {cond} ---')
        all_data[cond] = scan_condition(cond, base)
        all_stats[cond] = summarize(all_data[cond])
        all_corrs[cond] = compute_condition_correlations(all_stats[cond])

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
        for cond in CONDITIONS:
            ci = all_corrs[cond].get(corr_key)
            if ci:
                print(f'    {cond:>18}: r={ci["r"]:.3f} (p={ci["p"]:.4f}), '
                      f'rho={ci["rho"]:.3f} (p={ci["p_spearman"]:.4f}), n={ci["n"]}')
            else:
                print(f'    {cond:>18}: N/A')

    print('\n' + '=' * 80)
    print('  CONDITION EFFECT SUMMARY')
    print('=' * 80)

    for metric, label in [('frac_entrance', 'Entrance Fraction'),
                          ('frac_active_site', 'Active-Site Fraction'),
                          ('mean_dist', 'Mean S-L Distance'),
                          ('mean_confidence', 'Mean Confidence')]:
        print(f'\n  {label}:')
        for cond in CONDITIONS:
            vals = [all_stats[cond][v][metric] for v in VARIANT_ORDER
                    if v in all_stats[cond] and all_stats[cond][v].get(metric) is not None]
            if vals:
                print(f'    {cond:>18}: mean={np.mean(vals):.3f}, '
                      f'std={np.std(vals):.3f}, '
                      f'range=[{np.min(vals):.3f}, {np.max(vals):.3f}]')

    # Pairwise condition changes for entrance fraction
    print('\n  Per-variant entrance fraction changes (pairwise):')
    for c1, c2 in [('constrained', 'no_atp'), ('no_cofactor', 'no_atp'),
                    ('constrained', 'no_cofactor'), ('constrained', 'no_constraint')]:
        common = [v for v in VARIANT_ORDER
                  if all(v in all_stats[c] and all_stats[c][v]['n_models'] > 0
                         for c in [c1, c2])]
        if len(common) >= 3:
            c1_vals = np.array([all_stats[c1][v]['frac_entrance'] for v in common])
            c2_vals = np.array([all_stats[c2][v]['frac_entrance'] for v in common])
            delta = c2_vals - c1_vals
            t, p = stats.ttest_rel(c1_vals, c2_vals)
            r, rp = stats.pearsonr(list(c1_vals), list(c2_vals))
            print(f'    {c1} -> {c2}: mean delta = {np.mean(delta):+.3f} '
                  f'(std={np.std(delta):.3f}), paired t={t:.2f} p={p:.4f}, '
                  f'correlation r={r:.3f} (n={len(common)})')

    # Mg2+ specific analysis
    print('\n' + '=' * 80)
    print('  Mg2+ PROXIMITY ANALYSIS (No ATP condition)')
    print('=' * 80)
    for v in VARIANT_ORDER:
        s = all_stats['no_atp'].get(v)
        if s and s.get('mean_mg_s_dist') is not None:
            print(f'    {v:>22}: mean Mg-S dist = {s["mean_mg_s_dist"]:.1f} A, '
                  f'entrance frac = {s["frac_entrance"]:.2f}')

    print('\n--- Generating figures ---')
    plot_four_condition_comparison(all_stats, all_corrs)
    plot_no_atp_focus(all_stats, all_data)
    plot_pairwise_all(all_stats)

    # Save JSON
    output = {
        'conditions': {},
        'correlations': {},
    }
    for cond in CONDITIONS:
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

    json_path = os.path.join(OUTPUT_DIR, 'four_condition_comparison.json')
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f'Saved {json_path}')

    print('\n--- DONE ---')


if __name__ == '__main__':
    main()
