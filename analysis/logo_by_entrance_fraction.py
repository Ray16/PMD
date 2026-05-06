#!/usr/bin/env python3
"""Generate sequence logo plots comparing low vs high entrance fraction homologs.

Splits the 93 cofolded homologs into low entrance fraction (bottom third)
and high entrance fraction (top third) groups, aligns their sequences to
PMDsc, and generates comparative logo plots at catalytic and gateway residues.
"""

import os
import csv
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logomaker
from Bio.Align import PairwiseAligner
from collections import Counter

ANALYSIS_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'
CONFIG_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/cofold_configs'

PMDSC_SEQ = (
    "MTVYTASVTAPVNIATLKYWGKRDTKLNLPTNSSISVTLSQDDLRTLTSAATAPEFERDTLWLNGEPHSIDNERTQNCLRDLR"
    "QLRKEMESKDASLPTLSQWKLHIVSENNFPTAAGLASSAAGFAALVSAIAKLYQLPQSTSEISRIARKGSGSACRSLFGGYVAWEM"
    "GKAEDGHDSMAVQIADSSDWPQMKACVLVVSDIKKDVSSTQGMQLTVATSELFKERIEHVVPKRFEVMRKAIVEKDFATFAKETMM"
    "DSNSFHATCLDSFPPIFYMNDTSKRIISWCHTINQFYGETIVAYTFDAGPNAVLYYLAENESKLFAFIYKLFGSVPGWDKKFTTEQL"
    "EAFNHQFESSNFTARELDLELQKDVARVILTQVGSGPQETNESLIDAKTGLPKE"
)

CATALYTIC_RESIDUES = {
    18: 'K18', 19: 'Y19', 20: 'W20', 22: 'K22',
    120: 'S120', 121: 'S121', 153: 'S153', 155: 'S155',
    158: 'R158', 208: 'S208', 302: 'D302',
}

GATEWAY_RESIDUES = {
    74: 'R74', 147: 'R147', 212: 'M212',
}

ENTRANCE_CONTACT_RESIDUES = {
    25: 'T25', 28: 'N28', 72: 'N72', 73: 'G73',
    74: 'R74', 154: 'G154', 209: 'T209',
}


def get_sequence_from_yaml(yaml_path):
    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    for entry in config.get('sequences', []):
        if 'protein' in entry:
            return entry['protein']['sequence']
    return None


def align_and_extract(pmdsc_seq, target_seq, positions):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(pmdsc_seq, target_seq)[0]
    pmdsc_blocks, target_blocks = alignment.aligned

    mapping = {}
    for (ps, pe), (ts, te) in zip(pmdsc_blocks, target_blocks):
        for offset in range(pe - ps):
            mapping[ps + offset + 1] = ts + offset

    result = {}
    for pos in positions:
        if pos in mapping and mapping[pos] < len(target_seq):
            result[pos] = target_seq[mapping[pos]]
        else:
            result[pos] = '-'
    return result


def make_comparison_logo(rows_low, rows_high, positions, labels, title, outpath,
                         highlight_positions=None):
    n_pos = len(positions)
    aa_order = list('ACDEFGHIKLMNPQRSTVWY')

    fig, axes = plt.subplots(2, 2, figsize=(max(14, n_pos * 0.9), 10),
                              gridspec_kw={'height_ratios': [3, 1]})

    for col, (rows, group_label, color) in enumerate([
        (rows_low, f'Low entrance fraction (n={len(rows_low)})', '#2ca02c'),
        (rows_high, f'High entrance fraction (n={len(rows_high)})', '#d62728'),
    ]):
        counts = np.zeros((n_pos, len(aa_order)))
        for row in rows:
            for i, pos in enumerate(positions):
                aa = row.get(pos, '-')
                if aa in aa_order:
                    counts[i, aa_order.index(aa)] += 1

        totals = counts.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1
        freq = counts / totals

        with np.errstate(divide='ignore', invalid='ignore'):
            info = np.log2(20) + np.sum(freq * np.where(freq > 0, np.log2(freq), 0), axis=1)
        info = np.nan_to_num(info, nan=0.0)
        df_info = pd.DataFrame(freq, columns=aa_order, index=range(n_pos)).multiply(info, axis=0)

        logo = logomaker.Logo(df_info, ax=axes[0, col], color_scheme='chemistry',
                              font_name='DejaVu Sans Mono')
        axes[0, col].set_ylabel('Information (bits)', fontsize=10)
        axes[0, col].set_title(group_label, fontsize=12, fontweight='bold', color=color)
        axes[0, col].set_ylim(0, np.log2(20) + 0.3)

        if highlight_positions:
            for i, pos in enumerate(positions):
                if pos in highlight_positions:
                    axes[0, col].axvspan(i - 0.5, i + 0.5, alpha=0.15, color='red', zorder=0)

        pmdsc_aa = [PMDSC_SEQ[p - 1] if p <= len(PMDSC_SEQ) else '?' for p in positions]
        axes[0, col].set_xticks(range(n_pos))
        axes[0, col].set_xticklabels([f'{l}\n({aa})' for l, aa in zip(labels, pmdsc_aa)],
                                      fontsize=7, ha='center')

        conservation = info / np.log2(20)
        axes[1, col].bar(range(n_pos), conservation, color=color, alpha=0.7,
                         edgecolor='k', linewidth=0.5)
        axes[1, col].set_ylabel('Conservation', fontsize=10)
        axes[1, col].set_ylim(0, 1.05)
        axes[1, col].set_xticks(range(n_pos))
        axes[1, col].set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
        axes[1, col].axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    fig.suptitle(title, fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")


def make_difference_bar(rows_low, rows_high, positions, labels, outpath):
    """Bar chart showing conservation difference (low - high) at each position."""
    n_pos = len(positions)
    aa_order = list('ACDEFGHIKLMNPQRSTVWY')

    cons_vals = {}
    pmdsc_match = {}
    for group_name, rows in [('low', rows_low), ('high', rows_high)]:
        counts = np.zeros((n_pos, len(aa_order)))
        for row in rows:
            for i, pos in enumerate(positions):
                aa = row.get(pos, '-')
                if aa in aa_order:
                    counts[i, aa_order.index(aa)] += 1
        totals = counts.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1
        freq = counts / totals
        with np.errstate(divide='ignore', invalid='ignore'):
            info = np.log2(20) + np.sum(freq * np.where(freq > 0, np.log2(freq), 0), axis=1)
        info = np.nan_to_num(info, nan=0.0)
        cons_vals[group_name] = info / np.log2(20)

        match_frac = []
        for i, pos in enumerate(positions):
            pmdsc_aa = PMDSC_SEQ[pos - 1]
            idx = aa_order.index(pmdsc_aa) if pmdsc_aa in aa_order else -1
            if idx >= 0 and totals[i, 0] > 0:
                match_frac.append(counts[i, idx] / totals[i, 0])
            else:
                match_frac.append(0)
        pmdsc_match[group_name] = np.array(match_frac)

    diff = cons_vals['low'] - cons_vals['high']
    match_diff = pmdsc_match['low'] - pmdsc_match['high']

    fig, axes = plt.subplots(2, 1, figsize=(max(14, n_pos * 0.7), 8))

    colors = ['#2ca02c' if d > 0 else '#d62728' for d in diff]
    axes[0].bar(range(n_pos), diff, color=colors, alpha=0.8, edgecolor='k', linewidth=0.5)
    axes[0].axhline(y=0, color='black', linewidth=0.8)
    axes[0].set_ylabel('Conservation difference\n(low entrance - high entrance)', fontsize=10)
    axes[0].set_title('Conservation Difference at Functional Residues\n'
                      '(green = more conserved in low-entrance group)', fontsize=12, fontweight='bold')
    axes[0].set_xticks(range(n_pos))
    axes[0].set_xticklabels(labels, fontsize=8, rotation=45, ha='right')

    colors2 = ['#2ca02c' if d > 0 else '#d62728' for d in match_diff]
    axes[1].bar(range(n_pos), match_diff, color=colors2, alpha=0.8, edgecolor='k', linewidth=0.5)
    axes[1].axhline(y=0, color='black', linewidth=0.8)
    axes[1].set_ylabel('PMDsc identity difference\n(low entrance - high entrance)', fontsize=10)
    axes[1].set_title('PMDsc Residue Identity Difference\n'
                      '(green = more PMDsc-like in low-entrance group)', fontsize=12, fontweight='bold')
    axes[1].set_xticks(range(n_pos))
    axes[1].set_xticklabels(labels, fontsize=8, rotation=45, ha='right')

    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")


def main():
    csv_path = os.path.join(ANALYSIS_DIR, 'homolog_entrance_analysis.csv')
    homologs = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            homologs.append(row)

    print(f"Loaded {len(homologs)} homologs from entrance analysis")

    n = len(homologs)
    third = n // 3
    low_group = homologs[:third]
    high_group = homologs[-third:]

    print(f"\nLow entrance fraction group (n={len(low_group)}):")
    print(f"  Entrance range: {low_group[0]['frac_entrance']} - {low_group[-1]['frac_entrance']}")
    print(f"  Close range: {low_group[0]['frac_close']} - {low_group[-1]['frac_close']}")

    print(f"\nHigh entrance fraction group (n={len(high_group)}):")
    print(f"  Entrance range: {high_group[0]['frac_entrance']} - {high_group[-1]['frac_entrance']}")
    print(f"  Close range: {high_group[0]['frac_close']} - {high_group[-1]['frac_close']}")

    all_positions = sorted(set(
        list(CATALYTIC_RESIDUES.keys()) +
        list(GATEWAY_RESIDUES.keys()) +
        list(ENTRANCE_CONTACT_RESIDUES.keys())
    ))
    all_labels = []
    for p in all_positions:
        if p in CATALYTIC_RESIDUES:
            all_labels.append(CATALYTIC_RESIDUES[p])
        elif p in GATEWAY_RESIDUES:
            all_labels.append(GATEWAY_RESIDUES[p])
        elif p in ENTRANCE_CONTACT_RESIDUES:
            all_labels.append(ENTRANCE_CONTACT_RESIDUES[p])

    def load_group_alignments(group):
        rows = []
        for h in group:
            name = h['name']
            yaml_path = os.path.join(CONFIG_DIR, f'{name}_2mvap.yaml')
            if not os.path.exists(yaml_path):
                continue
            seq = get_sequence_from_yaml(yaml_path)
            if not seq:
                continue
            try:
                extracted = align_and_extract(PMDSC_SEQ, seq, all_positions)
                extracted['name'] = name
                extracted['frac_entrance'] = float(h['frac_entrance'])
                extracted['frac_close'] = float(h['frac_close'])
                extracted['kingdom'] = h['kingdom']
                rows.append(extracted)
            except Exception as e:
                print(f"  Skipping {name}: {e}")
        return rows

    print("\nAligning low-entrance group sequences...")
    rows_low = load_group_alignments(low_group)
    print(f"  Successfully aligned: {len(rows_low)}")

    print("Aligning high-entrance group sequences...")
    rows_high = load_group_alignments(high_group)
    print(f"  Successfully aligned: {len(rows_high)}")

    cat_positions = sorted(CATALYTIC_RESIDUES.keys())
    cat_labels = [CATALYTIC_RESIDUES[p] for p in cat_positions]

    gw_positions = sorted(GATEWAY_RESIDUES.keys())
    cat_gw_positions = sorted(set(cat_positions + gw_positions))
    cat_gw_labels = []
    for p in cat_gw_positions:
        if p in CATALYTIC_RESIDUES:
            cat_gw_labels.append(CATALYTIC_RESIDUES[p])
        else:
            cat_gw_labels.append(GATEWAY_RESIDUES[p])

    entrance_positions = sorted(ENTRANCE_CONTACT_RESIDUES.keys())
    entrance_labels = [ENTRANCE_CONTACT_RESIDUES[p] for p in entrance_positions]

    print("\nGenerating comparative logo plots...")

    make_comparison_logo(
        rows_low, rows_high, cat_gw_positions, cat_gw_labels,
        'Catalytic & Gateway Residue Conservation:\nLow vs High Entrance Fraction Homologs',
        os.path.join(ANALYSIS_DIR, 'logo_low_vs_high_entrance_cat_gateway.png'),
        highlight_positions=set(gw_positions),
    )

    make_comparison_logo(
        rows_low, rows_high, entrance_positions, entrance_labels,
        'Entrance Contact Residue Conservation:\nLow vs High Entrance Fraction Homologs',
        os.path.join(ANALYSIS_DIR, 'logo_low_vs_high_entrance_contacts.png'),
    )

    make_comparison_logo(
        rows_low, rows_high, all_positions, all_labels,
        'All Functional Residue Conservation:\nLow vs High Entrance Fraction Homologs',
        os.path.join(ANALYSIS_DIR, 'logo_low_vs_high_entrance_all.png'),
        highlight_positions=set(entrance_positions),
    )

    make_difference_bar(
        rows_low, rows_high, all_positions, all_labels,
        os.path.join(ANALYSIS_DIR, 'logo_entrance_conservation_difference.png'),
    )

    # Print summary table
    aa_order = list('ACDEFGHIKLMNPQRSTVWY')
    print(f"\n{'Position':<8} {'PMDsc':>6} {'Low top AA':>10} {'Low freq':>9} {'High top AA':>11} {'High freq':>10} {'Delta cons':>10}")
    print("-" * 70)
    for pos, label in zip(all_positions, all_labels):
        pmdsc_aa = PMDSC_SEQ[pos - 1]
        for group_name, rows in [('low', rows_low), ('high', rows_high)]:
            counts = Counter()
            for row in rows:
                aa = row.get(pos, '-')
                if aa in aa_order:
                    counts[aa] += 1
            total = sum(counts.values())
            if total > 0:
                top_aa, top_count = counts.most_common(1)[0]
                freq = top_count / total
                freqs_arr = np.array([counts.get(a, 0) / total for a in aa_order])
                with np.errstate(divide='ignore', invalid='ignore'):
                    info = np.log2(20) + np.sum(freqs_arr * np.where(freqs_arr > 0, np.log2(freqs_arr), 0))
                cons = max(0, info) / np.log2(20)
            else:
                top_aa, freq, cons = '-', 0, 0
            if group_name == 'low':
                low_aa, low_freq, low_cons = top_aa, freq, cons
            else:
                high_aa, high_freq, high_cons = top_aa, freq, cons

        delta = low_cons - high_cons
        print(f"{label:<8} {pmdsc_aa:>6} {low_aa:>8} {low_freq:>8.0%} {high_aa:>10} {high_freq:>9.0%} {delta:>+10.3f}")

    # Kingdom breakdown
    print(f"\nKingdom breakdown:")
    low_kingdoms = Counter(r['kingdom'] for r in rows_low)
    high_kingdoms = Counter(r['kingdom'] for r in rows_high)
    all_kingdoms = sorted(set(list(low_kingdoms.keys()) + list(high_kingdoms.keys())))
    print(f"  {'Kingdom':<12} {'Low entrance':>13} {'High entrance':>14}")
    for k in all_kingdoms:
        print(f"  {k:<12} {low_kingdoms.get(k, 0):>13} {high_kingdoms.get(k, 0):>14}")

    print("\nDone!")


if __name__ == '__main__':
    main()
