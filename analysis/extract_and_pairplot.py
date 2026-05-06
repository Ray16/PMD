#!/usr/bin/env python3
"""Extract all Boltz2 cofolding metrics from confidence JSONs across 4 conditions,
merge with experimental data, save to CSV, and generate pair plots."""

import json
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

BOLTZ_BASE = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2'
ANALYSIS_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'

CONDITION_DIRS = {
    'constrained':    'output_2mvap',
    'no_constraint':  'output_2mvap_no_constraint',
    'no_atp':         'output_2mvap_no_atp',
    'no_cofactor':    'output_2mvap_no_cofactor',
}

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

# Chain index mapping: 0=A(protein), 1=L(MVAP substrate), 2=L2(second MVAP/S), 3=T(ATP), 4=M(Mg2+)
CHAIN_NAMES = {0: 'A_protein', 1: 'L_substrate', 2: 'S_substrate2', 3: 'T_atp', 4: 'M_mg'}

VARIANT_DATA = {
    'WT':   {'ki': 18,   'kcat': 0.15, 'km': 2.3,  'kcat_km': 0.066, 'titer': 475,  'growth_rate': 0.39},
    'Y19H': {'ki': None, 'kcat': 0.27, 'km': 0.35, 'kcat_km': 0.78,  'titer': 388,  'growth_rate': 0.22},
    'K22M': {'ki': None, 'kcat': 0.09, 'km': 1.3,  'kcat_km': 0.12,  'titer': 22,   'growth_rate': None},
    'K22Y': {'ki': None, 'kcat': 0.09, 'km': 1.3,  'kcat_km': 0.12,  'titer': 22,   'growth_rate': 0.13},
    'R74G': {'ki': 110,  'kcat': 0.14, 'km': 3.4,  'kcat_km': 0.04,  'titer': 975,  'growth_rate': 0.81},
    'R74H': {'ki': None, 'kcat': 0.33, 'km': 0.75, 'kcat_km': 0.44,  'titer': 770,  'growth_rate': 0.68},
    'I145A': {'ki': None, 'kcat': 0.029,'km': 2.0,  'kcat_km': 0.01,  'titer': 623,  'growth_rate': 0.17},
    'I145F': {'ki': None, 'kcat': 0.28, 'km': 1.36, 'kcat_km': 0.20,  'titer': None, 'growth_rate': None},
    'R147K': {'ki': None, 'kcat': 0.149,'km': 0.5,  'kcat_km': 0.32,  'titer': 793,  'growth_rate': 0.56},
    'S186C': {'ki': None, 'kcat': 0.07, 'km': 0.8,  'kcat_km': 0.08,  'titer': 596,  'growth_rate': 0.30},
    'S208E': {'ki': None, 'kcat': None, 'km': None, 'kcat_km': None,  'titer': None, 'growth_rate': None},
    'T209D': {'ki': None, 'kcat': 0.13, 'km': 0.99, 'kcat_km': 0.13,  'titer': None, 'growth_rate': None},
    'M212Q': {'ki': None, 'kcat': 0.35, 'km': 0.7,  'kcat_km': 0.50,  'titer': 601,  'growth_rate': 0.09},
    'I226V': {'ki': None, 'kcat': 0.16, 'km': 0.34, 'kcat_km': 0.46,  'titer': 633,  'growth_rate': 0.35},
    'V230E': {'ki': 10,   'kcat': 0.07, 'km': 0.8,  'kcat_km': 0.08,  'titer': 278,  'growth_rate': 0.21},
    'R74G-R147K': {'ki': None, 'kcat': 0.22, 'km': 0.53, 'kcat_km': 0.42, 'titer': 909, 'growth_rate': 0.84},
    'HKQ':  {'ki': 80,   'kcat': 0.16, 'km': 0.40, 'kcat_km': 0.40,  'titer': 1079, 'growth_rate': 0.79},
    'GKQ':  {'ki': 11,   'kcat': 0.22, 'km': 0.50, 'kcat_km': 0.50,  'titer': 8,    'growth_rate': None},
    'R74G-R147K-Q140L': {'ki': None, 'kcat': 0.06, 'km': 2.0, 'kcat_km': 0.03, 'titer': 401, 'growth_rate': None},
}

# Also load structural metrics from the four_condition_comparison.json
with open(os.path.join(ANALYSIS_DIR, 'four_condition_comparison.json')) as f:
    four_cond = json.load(f)


def extract_confidence_metrics(condition, variant_key, variant_dir_name):
    """Extract per-model confidence metrics from Boltz2 confidence JSONs."""
    cond_dir = CONDITION_DIRS[condition]
    pred_dir = os.path.join(
        BOLTZ_BASE, cond_dir, variant_dir_name,
        f'boltz_results_{variant_dir_name}', 'predictions', variant_dir_name
    )

    if not os.path.isdir(pred_dir):
        return []

    conf_files = sorted(glob.glob(os.path.join(pred_dir, 'confidence_*.json')))
    rows = []

    for cf in conf_files:
        model_num = int(os.path.basename(cf).split('_model_')[1].replace('.json', ''))
        with open(cf) as f:
            cd = json.load(f)

        row = {
            'variant': variant_key,
            'condition': condition,
            'model': model_num,
            'confidence_score': cd['confidence_score'],
            'ptm': cd['ptm'],
            'iptm': cd['iptm'],
            'ligand_iptm': cd['ligand_iptm'],
            'protein_iptm': cd['protein_iptm'],
            'complex_plddt': cd['complex_plddt'],
            'complex_iplddt': cd['complex_iplddt'],
            'complex_pde': cd['complex_pde'],
            'complex_ipde': cd['complex_ipde'],
        }

        # Chain-level pTM
        for chain_idx_str, ptm_val in cd['chains_ptm'].items():
            chain_idx = int(chain_idx_str)
            chain_name = CHAIN_NAMES.get(chain_idx, f'chain{chain_idx}')
            row[f'chain_ptm_{chain_name}'] = ptm_val

        # Pairwise chain ipTM (extract key pairs)
        pair_iptm = cd['pair_chains_iptm']
        for i_str, targets in pair_iptm.items():
            i = int(i_str)
            i_name = CHAIN_NAMES.get(i, f'chain{i}')
            for j_str, val in targets.items():
                j = int(j_str)
                if i <= j:
                    j_name = CHAIN_NAMES.get(j, f'chain{j}')
                    row[f'pair_iptm_{i_name}__{j_name}'] = val

        rows.append(row)

    return rows


def main():
    print("Extracting per-model confidence metrics from all conditions...")
    all_rows = []

    for condition in CONDITION_DIRS:
        print(f"  Condition: {condition}")
        for variant_key, variant_dir_name in VARIANT_NAMES.items():
            rows = extract_confidence_metrics(condition, variant_key, variant_dir_name)
            if rows:
                print(f"    {variant_key}: {len(rows)} models")
            else:
                print(f"    {variant_key}: NO DATA")
            all_rows.extend(rows)

    df_models = pd.DataFrame(all_rows)
    print(f"\nTotal per-model rows: {len(df_models)}")
    print(f"Columns: {list(df_models.columns)}")

    # Save per-model CSV
    per_model_csv = os.path.join(ANALYSIS_DIR, 'boltz2_per_model_metrics.csv')
    df_models.to_csv(per_model_csv, index=False)
    print(f"\nSaved per-model metrics to {per_model_csv}")

    # Aggregate per variant-condition (mean across models)
    group_cols = ['variant', 'condition']
    numeric_cols = [c for c in df_models.columns if c not in group_cols + ['model']]
    df_agg = df_models.groupby(group_cols)[numeric_cols].mean().reset_index()

    # Add experimental data
    df_agg['ki'] = df_agg['variant'].map(lambda v: VARIANT_DATA[v]['ki'])
    df_agg['titer'] = df_agg['variant'].map(lambda v: VARIANT_DATA[v]['titer'])
    df_agg['kcat'] = df_agg['variant'].map(lambda v: VARIANT_DATA[v]['kcat'])
    df_agg['km'] = df_agg['variant'].map(lambda v: VARIANT_DATA[v]['km'])
    df_agg['kcat_km'] = df_agg['variant'].map(lambda v: VARIANT_DATA[v]['kcat_km'])
    df_agg['growth_rate'] = df_agg['variant'].map(lambda v: VARIANT_DATA[v]['growth_rate'])

    # Add structural metrics from four_condition_comparison.json
    structural_cols = ['mean_dist', 'std_dist', 'median_dist', 'frac_entrance',
                       'frac_active_site', 'frac_close', 'mean_contacts',
                       'mean_entrance_contacts', 'mean_mg_s_dist']
    for col in structural_cols:
        df_agg[col] = df_agg.apply(
            lambda r: four_cond['conditions'].get(r['condition'], {}).get(r['variant'], {}).get(col),
            axis=1
        )

    agg_csv = os.path.join(ANALYSIS_DIR, 'boltz2_aggregated_metrics.csv')
    df_agg.to_csv(agg_csv, index=False)
    print(f"Saved aggregated metrics to {agg_csv}")
    print(f"Aggregated shape: {df_agg.shape}")
    print(f"Aggregated columns: {list(df_agg.columns)}")

    # --- PAIR PLOTS ---
    # Define the cofolding/structural metrics to plot
    cofolding_metrics = [
        'confidence_score', 'ptm', 'iptm', 'ligand_iptm',
        'complex_plddt', 'complex_iplddt', 'complex_pde', 'complex_ipde',
        'chain_ptm_A_protein', 'chain_ptm_L_substrate', 'chain_ptm_S_substrate2', 'chain_ptm_T_atp',
        'pair_iptm_A_protein__L_substrate', 'pair_iptm_A_protein__S_substrate2',
        'pair_iptm_A_protein__T_atp', 'pair_iptm_L_substrate__S_substrate2',
        'pair_iptm_L_substrate__T_atp', 'pair_iptm_S_substrate2__T_atp',
        'mean_dist', 'median_dist', 'frac_entrance', 'frac_active_site',
        'frac_close', 'mean_contacts', 'mean_entrance_contacts', 'mean_mg_s_dist',
    ]

    experimental_targets = ['titer', 'ki', 'kcat', 'km', 'kcat_km', 'growth_rate']

    # Available cofolding metrics (some pair_iptm cols may not exist in all conditions)
    available_metrics = [m for m in cofolding_metrics if m in df_agg.columns]

    for condition in CONDITION_DIRS:
        print(f"\nGenerating pair plots for {condition}...")
        df_cond = df_agg[df_agg['condition'] == condition].copy()

        # --- 1. Pair plot: cofolding metrics vs titer ---
        df_titer = df_cond.dropna(subset=['titer'])
        if len(df_titer) >= 4:
            plot_pairwise_correlations(
                df_titer, available_metrics, 'titer',
                f'Cofolding Metrics vs Cell Growth Titer ({condition})',
                os.path.join(ANALYSIS_DIR, f'pairplot_titer_{condition}.png'),
                label_col='variant'
            )

        # --- 2. Pair plot: cofolding metrics vs Ki ---
        df_ki = df_cond.dropna(subset=['ki'])
        if len(df_ki) >= 3:
            plot_pairwise_correlations(
                df_ki, available_metrics, 'ki',
                f'Cofolding Metrics vs Ki ({condition})',
                os.path.join(ANALYSIS_DIR, f'pairplot_ki_{condition}.png'),
                label_col='variant'
            )

        # --- 3. Pair plot: cofolding metrics vs kcat_km ---
        df_kcatkm = df_cond.dropna(subset=['kcat_km'])
        if len(df_kcatkm) >= 4:
            plot_pairwise_correlations(
                df_kcatkm, available_metrics, 'kcat_km',
                f'Cofolding Metrics vs kcat/Km ({condition})',
                os.path.join(ANALYSIS_DIR, f'pairplot_kcat_km_{condition}.png'),
                label_col='variant'
            )

        # --- 4. Pair plot: cofolding metrics vs growth_rate ---
        df_gr = df_cond.dropna(subset=['growth_rate'])
        if len(df_gr) >= 4:
            plot_pairwise_correlations(
                df_gr, available_metrics, 'growth_rate',
                f'Cofolding Metrics vs Growth Rate ({condition})',
                os.path.join(ANALYSIS_DIR, f'pairplot_growth_rate_{condition}.png'),
                label_col='variant'
            )

    # --- 4. Cross-condition comparison: one big pair plot per experimental target ---
    for target in ['titer', 'ki', 'growth_rate']:
        df_target = df_agg.dropna(subset=[target])
        if len(df_target) >= 4:
            plot_pairwise_correlations_by_condition(
                df_target, available_metrics, target,
                f'Cofolding Metrics vs {target.upper()} (All Conditions)',
                os.path.join(ANALYSIS_DIR, f'pairplot_{target}_all_conditions.png')
            )

    # --- 5. Correlation heatmap across all metrics x conditions ---
    plot_correlation_heatmap(df_agg, available_metrics, experimental_targets,
                            os.path.join(ANALYSIS_DIR, 'correlation_heatmap.png'))

    print("\nDone! All outputs saved to", ANALYSIS_DIR)


def plot_pairwise_correlations(df, metrics, target, title, outpath, label_col='variant'):
    """Scatter plot grid: each cofolding metric vs a single experimental target."""
    valid_metrics = [m for m in metrics if m in df.columns and df[m].notna().sum() >= 3
                     and df[m].dropna().std() > 1e-12]
    n = len(valid_metrics)
    if n == 0:
        return

    ncols = 5
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
    axes = np.array(axes).flatten()

    for idx, metric in enumerate(valid_metrics):
        ax = axes[idx]
        x = df[metric].values
        y = df[target].values
        mask = np.isfinite(x.astype(float)) & np.isfinite(y.astype(float))
        x_clean, y_clean = x[mask].astype(float), y[mask].astype(float)

        ax.scatter(x_clean, y_clean, s=40, alpha=0.8, edgecolors='k', linewidth=0.5)

        if label_col in df.columns:
            for _, row in df[mask].iterrows():
                ax.annotate(row[label_col], (float(row[metric]), float(row[target])),
                            fontsize=5, alpha=0.7, ha='center', va='bottom')

        if len(x_clean) >= 3 and np.std(x_clean) > 1e-12 and np.std(y_clean) > 1e-12:
            r, p = stats.pearsonr(x_clean, y_clean)
            rho, p_s = stats.spearmanr(x_clean, y_clean)
            ax.set_title(f'{metric}\nr={r:.2f} (p={p:.3f})\nρ={rho:.2f} (p={p_s:.3f})', fontsize=7)

            try:
                z = np.polyfit(x_clean, y_clean, 1)
                xline = np.linspace(x_clean.min(), x_clean.max(), 100)
                ax.plot(xline, np.polyval(z, xline), 'r--', alpha=0.5, linewidth=1)
            except np.linalg.LinAlgError:
                pass
        else:
            ax.set_title(metric, fontsize=7)

        ax.set_xlabel(metric, fontsize=6)
        ax.set_ylabel(target, fontsize=7)
        ax.tick_params(labelsize=6)

    for idx in range(n, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(title, fontsize=12, fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_pairwise_correlations_by_condition(df, metrics, target, title, outpath):
    """Scatter plot grid colored by condition."""
    valid_metrics = [m for m in metrics if m in df.columns and df[m].notna().sum() >= 3
                     and df[m].dropna().std() > 1e-12]
    n = len(valid_metrics)
    if n == 0:
        return

    ncols = 5
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
    axes = np.array(axes).flatten()

    conditions = df['condition'].unique()
    colors = dict(zip(conditions, sns.color_palette('Set2', len(conditions))))

    for idx, metric in enumerate(valid_metrics):
        ax = axes[idx]
        for cond in conditions:
            dc = df[df['condition'] == cond]
            x = dc[metric].values.astype(float)
            y = dc[target].values.astype(float)
            mask = np.isfinite(x) & np.isfinite(y)
            ax.scatter(x[mask], y[mask], s=30, alpha=0.7, label=cond,
                       color=colors[cond], edgecolors='k', linewidth=0.3)

        all_x = df[metric].values.astype(float)
        all_y = df[target].values.astype(float)
        mask = np.isfinite(all_x) & np.isfinite(all_y)
        if mask.sum() >= 3 and np.std(all_x[mask]) > 1e-12 and np.std(all_y[mask]) > 1e-12:
            r, p = stats.pearsonr(all_x[mask], all_y[mask])
            rho, p_s = stats.spearmanr(all_x[mask], all_y[mask])
            ax.set_title(f'{metric}\nr={r:.2f} (p={p:.3f}) ρ={rho:.2f} (p={p_s:.3f})', fontsize=7)
        else:
            ax.set_title(metric, fontsize=7)

        ax.set_xlabel(metric, fontsize=6)
        ax.set_ylabel(target, fontsize=7)
        ax.tick_params(labelsize=6)

    for idx in range(n, len(axes)):
        axes[idx].set_visible(False)

    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[c],
                           markersize=8, label=c) for c in conditions]
    fig.legend(handles=handles, loc='upper right', fontsize=8)
    fig.suptitle(title, fontsize=12, fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_correlation_heatmap(df_agg, metrics, targets, outpath):
    """Heatmap of Spearman correlations: metrics × targets, split by condition."""
    conditions = df_agg['condition'].unique()
    valid_metrics = [m for m in metrics if m in df_agg.columns]

    fig, axes = plt.subplots(1, len(conditions), figsize=(6 * len(conditions), max(8, len(valid_metrics) * 0.35)))

    for ax_idx, cond in enumerate(conditions):
        ax = axes[ax_idx]
        dc = df_agg[df_agg['condition'] == cond]

        corr_matrix = np.full((len(valid_metrics), len(targets)), np.nan)
        pval_matrix = np.full((len(valid_metrics), len(targets)), np.nan)

        for i, metric in enumerate(valid_metrics):
            for j, target in enumerate(targets):
                x = dc[metric].dropna()
                y = dc[target].dropna()
                common = x.index.intersection(y.index)
                if len(common) >= 3:
                    rho, p = stats.spearmanr(dc.loc[common, metric].astype(float),
                                             dc.loc[common, target].astype(float))
                    corr_matrix[i, j] = rho
                    pval_matrix[i, j] = p

        annot = []
        for i in range(len(valid_metrics)):
            row = []
            for j in range(len(targets)):
                if np.isnan(corr_matrix[i, j]):
                    row.append('')
                else:
                    star = ''
                    if pval_matrix[i, j] < 0.001:
                        star = '***'
                    elif pval_matrix[i, j] < 0.01:
                        star = '**'
                    elif pval_matrix[i, j] < 0.05:
                        star = '*'
                    row.append(f'{corr_matrix[i, j]:.2f}{star}')
            annot.append(row)

        sns.heatmap(corr_matrix, ax=ax, xticklabels=targets, yticklabels=valid_metrics,
                    annot=annot, fmt='', cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                    cbar_kws={'shrink': 0.5})
        ax.set_title(cond, fontsize=10, fontweight='bold')
        ax.tick_params(labelsize=7)

    fig.suptitle('Spearman Correlations: Cofolding Metrics vs Experimental (* p<0.05, ** p<0.01, *** p<0.001)',
                 fontsize=11, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {outpath}")


if __name__ == '__main__':
    main()
