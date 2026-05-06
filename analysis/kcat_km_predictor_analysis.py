#!/usr/bin/env python3
"""Systematic analysis of which Boltz2 cofolding metrics best predict kcat/Km.

Ranks all metrics by Pearson and Spearman correlation across 4 conditions,
generates focused top-predictor scatter plots, tests composite metrics
for titer prediction (Ki × kcat/Km), and outputs a summary CSV.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

ANALYSIS_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'

COFOLDING_METRICS = [
    'confidence_score', 'ptm', 'iptm', 'ligand_iptm',
    'complex_plddt', 'complex_iplddt', 'complex_pde', 'complex_ipde',
    'chain_ptm_A_protein', 'chain_ptm_L_substrate', 'chain_ptm_S_substrate2',
    'chain_ptm_T_atp',
    'pair_iptm_A_protein__L_substrate', 'pair_iptm_A_protein__S_substrate2',
    'pair_iptm_A_protein__T_atp',
    'pair_iptm_L_substrate__L_substrate',
    'pair_iptm_L_substrate__S_substrate2', 'pair_iptm_L_substrate__T_atp',
    'pair_iptm_S_substrate2__S_substrate2',
    'pair_iptm_S_substrate2__T_atp',
    'mean_dist', 'median_dist', 'frac_entrance', 'frac_active_site',
    'frac_close', 'mean_contacts', 'mean_entrance_contacts', 'mean_mg_s_dist',
]

CONDITIONS = ['constrained', 'no_constraint', 'no_atp', 'no_cofactor']


def rank_correlations(df_agg, target='kcat_km'):
    """Rank all metrics by correlation with target, per condition."""
    rows = []
    for condition in CONDITIONS:
        dc = df_agg[(df_agg['condition'] == condition) & df_agg[target].notna()].copy()
        available = [m for m in COFOLDING_METRICS
                     if m in dc.columns and dc[m].notna().sum() >= 4
                     and dc[m].dropna().std() > 1e-12]

        for metric in available:
            valid = dc[[metric, target]].dropna()
            if len(valid) < 4:
                continue
            x = valid[metric].astype(float).values
            y = valid[target].astype(float).values
            r, p_r = stats.pearsonr(x, y)
            rho, p_s = stats.spearmanr(x, y)
            rows.append({
                'condition': condition,
                'metric': metric,
                'n': len(valid),
                'pearson_r': r,
                'pearson_p': p_r,
                'spearman_rho': rho,
                'spearman_p': p_s,
                'abs_spearman': abs(rho),
            })

    df_rank = pd.DataFrame(rows)
    df_rank = df_rank.sort_values(['condition', 'abs_spearman'], ascending=[True, False])
    return df_rank


def plot_top_predictors(df_agg, df_rank, target='kcat_km', top_n=5):
    """Scatter plots: top N metrics per condition vs kcat/Km."""
    fig, axes = plt.subplots(len(CONDITIONS), top_n, figsize=(4 * top_n, 4 * len(CONDITIONS)))

    for row_idx, condition in enumerate(CONDITIONS):
        dc = df_agg[(df_agg['condition'] == condition) & df_agg[target].notna()].copy()
        cond_rank = df_rank[df_rank['condition'] == condition].head(top_n)

        for col_idx, (_, mrow) in enumerate(cond_rank.iterrows()):
            ax = axes[row_idx, col_idx]
            metric = mrow['metric']

            valid = dc[[metric, target, 'variant']].dropna()
            x = valid[metric].astype(float).values
            y = valid[target].astype(float).values

            ax.scatter(x, y, s=50, alpha=0.8, edgecolors='k', linewidth=0.5, zorder=5)

            for _, vrow in valid.iterrows():
                ax.annotate(vrow['variant'], (float(vrow[metric]), float(vrow[target])),
                            fontsize=5, alpha=0.7, ha='center', va='bottom')

            r, p_r = mrow['pearson_r'], mrow['pearson_p']
            rho, p_s = mrow['spearman_rho'], mrow['spearman_p']

            if len(x) >= 3 and np.std(x) > 1e-12:
                z = np.polyfit(x, y, 1)
                xline = np.linspace(x.min(), x.max(), 100)
                ax.plot(xline, np.polyval(z, xline), 'r--', alpha=0.5, linewidth=1.5)

            title = f'{metric}\nr={r:.2f} (p={p_r:.3f})\nρ={rho:.2f} (p={p_s:.3f})'
            ax.set_title(title, fontsize=7, fontweight='bold')
            ax.set_xlabel(metric, fontsize=6)
            if col_idx == 0:
                ax.set_ylabel(f'{target}', fontsize=8)
            ax.tick_params(labelsize=6)

        axes[row_idx, 0].text(-0.3, 0.5, condition,
                              transform=axes[row_idx, 0].transAxes,
                              fontsize=10, fontweight='bold', rotation=90,
                              ha='center', va='center')

    fig.suptitle(f'Top {top_n} Cofolding Metric Predictors of kcat/Km by Condition\n(Ranked by |Spearman ρ|)',
                 fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    outpath = os.path.join(ANALYSIS_DIR, 'kcat_km_top_predictors.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")


def plot_condition_comparison(df_rank, target='kcat_km'):
    """Heatmap: top metrics across conditions, showing Spearman rho."""
    top_metrics_per_cond = []
    for cond in CONDITIONS:
        top = df_rank[df_rank['condition'] == cond].head(8)['metric'].tolist()
        top_metrics_per_cond.extend(top)
    unique_metrics = list(dict.fromkeys(top_metrics_per_cond))

    matrix = np.full((len(unique_metrics), len(CONDITIONS)), np.nan)
    pval_matrix = np.full((len(unique_metrics), len(CONDITIONS)), np.nan)

    for j, cond in enumerate(CONDITIONS):
        cond_data = df_rank[df_rank['condition'] == cond]
        for i, metric in enumerate(unique_metrics):
            match = cond_data[cond_data['metric'] == metric]
            if not match.empty:
                matrix[i, j] = match.iloc[0]['spearman_rho']
                pval_matrix[i, j] = match.iloc[0]['spearman_p']

    fig, ax = plt.subplots(figsize=(8, max(6, len(unique_metrics) * 0.4)))

    annot = []
    for i in range(len(unique_metrics)):
        row = []
        for j in range(len(CONDITIONS)):
            if np.isnan(matrix[i, j]):
                row.append('')
            else:
                star = ''
                if pval_matrix[i, j] < 0.001:
                    star = '***'
                elif pval_matrix[i, j] < 0.01:
                    star = '**'
                elif pval_matrix[i, j] < 0.05:
                    star = '*'
                row.append(f'{matrix[i, j]:.2f}{star}')
        annot.append(row)

    import seaborn as sns
    sns.heatmap(matrix, ax=ax, xticklabels=CONDITIONS, yticklabels=unique_metrics,
                annot=annot, fmt='', cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                cbar_kws={'shrink': 0.7, 'label': 'Spearman ρ'})
    ax.set_title(f'Top Cofolding Metric Correlations with {target}\n(* p<0.05, ** p<0.01, *** p<0.001)',
                 fontsize=11, fontweight='bold')
    ax.tick_params(labelsize=8)

    outpath = os.path.join(ANALYSIS_DIR, 'kcat_km_predictor_heatmap.png')
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")


def test_composite_metrics(df_agg):
    """Test whether composite metrics predict titer better than individual ones.

    Since titer depends on both Ki and kcat/Km, combinations might be more
    predictive. Limited to variants with measured Ki (n=5).
    """
    results = []

    for condition in CONDITIONS:
        dc = df_agg[df_agg['condition'] == condition].copy()

        # Test individual metrics vs titer (all variants with titer)
        dc_titer = dc.dropna(subset=['titer'])
        available = [m for m in COFOLDING_METRICS
                     if m in dc_titer.columns and dc_titer[m].notna().sum() >= 4
                     and dc_titer[m].dropna().std() > 1e-12]

        for metric in available:
            valid = dc_titer[[metric, 'titer']].dropna()
            if len(valid) < 4:
                continue
            rho, p = stats.spearmanr(valid[metric].astype(float), valid['titer'].astype(float))
            results.append({
                'condition': condition,
                'metric': metric,
                'target': 'titer',
                'type': 'individual',
                'n': len(valid),
                'spearman_rho': rho,
                'spearman_p': p,
            })

        # For variants with Ki: test entrance% alone, kcat/km predictor alone,
        # and composite = -entrance% + best_kcat_km_metric
        dc_ki = dc.dropna(subset=['ki', 'titer'])
        if len(dc_ki) < 4:
            continue

        if 'frac_entrance' in dc_ki.columns and dc_ki['frac_entrance'].notna().sum() >= 4:
            valid = dc_ki[['frac_entrance', 'titer']].dropna()
            if len(valid) >= 4:
                rho, p = stats.spearmanr(valid['frac_entrance'].astype(float),
                                          valid['titer'].astype(float))
                results.append({
                    'condition': condition,
                    'metric': 'frac_entrance',
                    'target': 'titer (Ki-subset)',
                    'type': 'individual',
                    'n': len(valid),
                    'spearman_rho': rho,
                    'spearman_p': p,
                })

        # Composite: try combining entrance metrics with substrate-binding metrics
        for sub_metric in available:
            if sub_metric in ('frac_entrance', 'frac_close', 'frac_active_site'):
                continue
            valid = dc_ki[['frac_entrance', sub_metric, 'titer']].dropna()
            if len(valid) < 4:
                continue
            fe = valid['frac_entrance'].astype(float).values
            sm = valid[sub_metric].astype(float).values
            titer = valid['titer'].astype(float).values

            if np.std(fe) < 1e-12 or np.std(sm) < 1e-12:
                continue

            fe_norm = (fe - fe.mean()) / fe.std()
            sm_norm = (sm - sm.mean()) / sm.std()
            composite = -fe_norm + sm_norm

            rho, p = stats.spearmanr(composite, titer)
            results.append({
                'condition': condition,
                'metric': f'composite(-entrance+{sub_metric})',
                'target': 'titer (Ki-subset)',
                'type': 'composite',
                'n': len(valid),
                'spearman_rho': rho,
                'spearman_p': p,
            })

    df_comp = pd.DataFrame(results)
    df_comp['abs_rho'] = df_comp['spearman_rho'].abs()
    df_comp = df_comp.sort_values(['condition', 'target', 'abs_rho'], ascending=[True, True, False])
    return df_comp


def main():
    csv_path = os.path.join(ANALYSIS_DIR, 'boltz2_aggregated_metrics.csv')
    df_agg = pd.read_csv(csv_path)
    print(f"Loaded {len(df_agg)} rows from {csv_path}")
    print(f"Variants with kcat_km: {df_agg.dropna(subset=['kcat_km'])['variant'].nunique()}")

    # 1. Rank all metrics by correlation with kcat/Km
    print("\n=== Ranking metrics by correlation with kcat/Km ===")
    df_rank = rank_correlations(df_agg, target='kcat_km')

    print("\nTop 5 per condition:")
    for cond in CONDITIONS:
        print(f"\n  {cond}:")
        top = df_rank[df_rank['condition'] == cond].head(5)
        for _, row in top.iterrows():
            print(f"    {row['metric']:45s}  ρ={row['spearman_rho']:+.3f} (p={row['spearman_p']:.3f})  "
                  f"r={row['pearson_r']:+.3f} (p={row['pearson_p']:.3f})  n={row['n']}")

    # Save ranking CSV
    rank_csv = os.path.join(ANALYSIS_DIR, 'kcat_km_metric_rankings.csv')
    df_rank.to_csv(rank_csv, index=False)
    print(f"\nSaved metric rankings to {rank_csv}")

    # 2. Top predictors scatter plots
    print("\n=== Generating top predictor plots ===")
    plot_top_predictors(df_agg, df_rank, target='kcat_km', top_n=5)

    # 3. Condition comparison heatmap
    print("\n=== Generating condition comparison heatmap ===")
    plot_condition_comparison(df_rank, target='kcat_km')

    # 4. Composite metrics for titer prediction
    print("\n=== Testing composite metrics for titer prediction ===")
    df_comp = test_composite_metrics(df_agg)

    comp_csv = os.path.join(ANALYSIS_DIR, 'composite_titer_predictors.csv')
    df_comp.to_csv(comp_csv, index=False)
    print(f"Saved composite analysis to {comp_csv}")

    # Show best composites per condition
    print("\nBest composite predictors of titer (Ki-subset):")
    composites = df_comp[(df_comp['type'] == 'composite') & (df_comp['target'] == 'titer (Ki-subset)')]
    for cond in CONDITIONS:
        top = composites[composites['condition'] == cond].head(3)
        if not top.empty:
            print(f"\n  {cond}:")
            for _, row in top.iterrows():
                print(f"    {row['metric']:60s}  ρ={row['spearman_rho']:+.3f} (p={row['spearman_p']:.3f})")

    # Also rank metrics for other experimental targets
    print("\n\n=== Ranking metrics for all experimental targets ===")
    all_rankings = []
    for target in ['kcat_km', 'kcat', 'km', 'titer', 'growth_rate', 'ki']:
        df_t = rank_correlations(df_agg, target=target)
        df_t['target'] = target
        all_rankings.append(df_t)

    df_all = pd.concat(all_rankings, ignore_index=True)
    all_csv = os.path.join(ANALYSIS_DIR, 'all_metric_rankings.csv')
    df_all.to_csv(all_csv, index=False)
    print(f"Saved all-target rankings to {all_csv}")

    # Summary: best metric per target per condition
    print("\n=== SUMMARY: Best predictor per target per condition ===")
    print(f"{'Target':<15} {'Condition':<18} {'Best Metric':<45} {'Spearman':>10} {'p-value':>10}")
    print("-" * 100)
    for target in ['kcat_km', 'kcat', 'km', 'titer', 'growth_rate', 'ki']:
        dt = df_all[df_all['target'] == target]
        for cond in CONDITIONS:
            top = dt[dt['condition'] == cond]
            if not top.empty:
                best = top.iloc[0]
                print(f"{target:<15} {cond:<18} {best['metric']:<45} {best['spearman_rho']:>+10.3f} {best['spearman_p']:>10.4f}")

    print("\nDone!")


if __name__ == '__main__':
    main()
