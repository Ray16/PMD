#!/usr/bin/env python3
"""Generate sequence logo plots for catalytic site conservation across
PMDsc structural homologs from Foldseek results.

Aligns all homolog sequences to PMDsc, extracts columns at catalytic
residue positions, and generates logo plots showing conservation.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logomaker
from Bio.Align import PairwiseAligner
from collections import Counter

ANALYSIS_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'
FOLDSEEK_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/foldseek_results'

PMDSC_SEQ = (
    "MTVYTASVTAPVNIATLKYWGKRDTKLNLPTNSSISVTLSQDDLRTLTSAATAPEFERDTLWLNGEPHSIDNERTQNCLRDLR"
    "QLRKEMESKDASLPTLSQWKLHIVSENNFPTAAGLASSAAGFAALVSAIAKLYQLPQSTSEISRIARKGSGSACRSLFGGYVAWEM"
    "GKAEDGHDSMAVQIADSSDWPQMKACVLVVSDIKKDVSSTQGMQLTVATSELFKERIEHVVPKRFEVMRKAIVEKDFATFAKETMM"
    "DSNSFHATCLDSFPPIFYMNDTSKRIISWCHTINQFYGETIVAYTFDAGPNAVLYYLAENESKLFAFIYKLFGSVPGWDKKFTTEQL"
    "EAFNHQFESSNFTARELDLELQKDVARVILTQVGSGPQETNESLIDAKTGLPKE"
)

CATALYTIC_RESIDUES = {
    18: 'K18',   19: 'Y19',   20: 'W20',   22: 'K22',
    120: 'S120', 121: 'S121', 153: 'S153', 155: 'S155',
    158: 'R158', 208: 'S208', 302: 'D302',
}

# Gateway / entrance residues (not catalytic but functionally important)
GATEWAY_RESIDUES = {
    74: 'R74', 147: 'R147', 212: 'M212',
}

# Extend to include nearby context (±2 residues around each catalytic site)
CONTEXT_WINDOW = 2

KINGDOM_COLORS = {
    'Fungi': '#2ca02c',
    'Animal': '#1f77b4',
    'Plant': '#8c564b',
    'Protozoa': '#9467bd',
    'Bacteria': '#d62728',
    'Archaea': '#ff7f0e',
    'Unknown': '#7f7f7f',
}


def load_all_sequences():
    """Load all homolog sequences from Foldseek m8 results."""
    sequences = {}

    # PMDsc itself
    sequences['PMDsc_WT'] = {'seq': PMDSC_SEQ, 'organism': 'Saccharomyces cerevisiae',
                              'kingdom': 'Fungi', 'seq_identity': 100.0}

    # Parse Foldseek results
    for fname in ['alis_pdb100.m8', 'alis_afdb-swissprot.m8']:
        fpath = os.path.join(FOLDSEEK_DIR, fname)
        if not os.path.exists(fpath):
            continue
        with open(fpath) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 21:
                    continue

                target_raw = fields[1].split(' ')[0]
                seq_identity = float(fields[2])
                aln_len = int(fields[3])
                t_full_seq = fields[18].replace('-', '') if len(fields) > 18 else ''
                organism = fields[20] if len(fields) > 20 else ''

                if not t_full_seq or len(t_full_seq) < 150 or aln_len < 200:
                    continue
                if seq_identity > 95 or seq_identity < 15:
                    continue

                # Clean ID
                if target_raw.startswith('AF-'):
                    clean_id = target_raw.split('-')[1]
                else:
                    clean_id = target_raw.split('-')[0].upper()

                if clean_id in sequences:
                    continue

                # Classify kingdom
                kingdom = 'Unknown'
                for genus, k in [
                    ('Saccharomyces', 'Fungi'), ('Eremothecium', 'Fungi'),
                    ('Debaryomyces', 'Fungi'), ('Candida', 'Fungi'),
                    ('Fusarium', 'Fungi'), ('Schizosaccharomyces', 'Fungi'),
                    ('Aspergillus', 'Fungi'), ('Ganoderma', 'Fungi'),
                    ('Cryptococcus', 'Fungi'),
                    ('Mus', 'Animal'), ('Bos', 'Animal'), ('Rattus', 'Animal'),
                    ('Danio', 'Animal'), ('Homo', 'Animal'),
                    ('Trypanosoma', 'Protozoa'), ('Dictyostelium', 'Protozoa'),
                    ('Arabidopsis', 'Plant'), ('Panax', 'Plant'),
                    ('Saccharolobus', 'Archaea'), ('Sulfolobus', 'Archaea'),
                    ('Thermoplasma', 'Archaea'), ('Picrophilus', 'Archaea'),
                    ('Haloferax', 'Archaea'),
                    ('Staphylococcus', 'Bacteria'), ('Enterococcus', 'Bacteria'),
                    ('Streptococcus', 'Bacteria'), ('Bacillus', 'Bacteria'),
                    ('Listeria', 'Bacteria'), ('Anaerolinea', 'Bacteria'),
                    ('Legionella', 'Bacteria'), ('Synechococcus', 'Bacteria'),
                    ('Prochlorococcus', 'Bacteria'), ('Clostridium', 'Bacteria'),
                    ('Rhodococcus', 'Bacteria'), ('Corynebacterium', 'Bacteria'),
                ]:
                    if genus.lower() in organism.lower():
                        kingdom = k
                        break

                sequences[clean_id] = {
                    'seq': t_full_seq,
                    'organism': organism,
                    'kingdom': kingdom,
                    'seq_identity': seq_identity,
                }

    return sequences


def align_and_extract(pmdsc_seq, target_seq):
    """Align target to PMDsc and return column mapping."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(pmdsc_seq, target_seq)[0]

    # Walk through alignment to build position mapping
    aln_str = str(alignment).split('\n')
    aln_q = aln_str[0]
    aln_t = aln_str[2]

    mapping = {}
    q_pos, t_pos = 0, 0
    for i in range(min(len(aln_q), len(aln_t))):
        q_gap = aln_q[i] in ('-', '.')
        t_gap = aln_t[i] in ('-', '.')
        if not q_gap and not t_gap:
            mapping[q_pos + 1] = t_pos  # q_resnum (1-indexed) -> t_pos (0-indexed)
            q_pos += 1
            t_pos += 1
        elif q_gap:
            t_pos += 1
        elif t_gap:
            q_pos += 1

    return mapping


def build_alignment_matrix(sequences, residue_positions):
    """Build a matrix of amino acids at specified PMDsc positions across all sequences."""
    results = []

    for seq_id, info in sequences.items():
        if seq_id == 'PMDsc_WT':
            row = {'id': seq_id, 'kingdom': info['kingdom'],
                   'identity': info['seq_identity']}
            for pos in residue_positions:
                if pos <= len(PMDSC_SEQ):
                    row[pos] = PMDSC_SEQ[pos - 1]
                else:
                    row[pos] = '-'
            results.append(row)
            continue

        try:
            mapping = align_and_extract(PMDSC_SEQ, info['seq'])
        except Exception:
            continue

        row = {'id': seq_id, 'kingdom': info['kingdom'],
               'identity': info['seq_identity']}
        for pos in residue_positions:
            if pos in mapping:
                t_idx = mapping[pos]
                if t_idx < len(info['seq']):
                    row[pos] = info['seq'][t_idx]
                else:
                    row[pos] = '-'
            else:
                row[pos] = '-'
        results.append(row)

    return results


def make_logo(alignment_rows, positions, labels, title, outpath,
              highlight_positions=None):
    """Generate a sequence logo from alignment data at specified positions."""
    n_pos = len(positions)
    aa_order = list('ACDEFGHIKLMNPQRSTVWY')

    # Build count matrix
    counts = np.zeros((n_pos, len(aa_order)))
    for row in alignment_rows:
        for i, pos in enumerate(positions):
            aa = row.get(pos, '-')
            if aa in aa_order:
                counts[i, aa_order.index(aa)] += 1

    # Convert to frequency
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    freq = counts / totals

    # Create dataframe for logomaker (integer index required)
    df = pd.DataFrame(freq, columns=aa_order, index=range(n_pos))

    # Information content (bits)
    with np.errstate(divide='ignore', invalid='ignore'):
        info = np.log2(20) + np.sum(freq * np.where(freq > 0, np.log2(freq), 0), axis=1)
    info = np.nan_to_num(info, nan=0.0)
    df_info = df.multiply(info, axis=0)

    fig, axes = plt.subplots(2, 1, figsize=(max(12, n_pos * 0.8), 7),
                              gridspec_kw={'height_ratios': [3, 1]})

    # Logo plot
    logo = logomaker.Logo(df_info, ax=axes[0], color_scheme='chemistry',
                          font_name='DejaVu Sans Mono')
    axes[0].set_ylabel('Information (bits)', fontsize=11)
    axes[0].set_title(title, fontsize=13, fontweight='bold', pad=10)
    axes[0].set_ylim(0, np.log2(20) + 0.3)

    # Highlight catalytic vs gateway positions
    if highlight_positions:
        for i, pos in enumerate(positions):
            if pos in highlight_positions:
                axes[0].axvspan(i - 0.5, i + 0.5, alpha=0.15, color='red', zorder=0)

    # PMDsc reference annotation
    pmdsc_aa = [PMDSC_SEQ[p - 1] if p <= len(PMDSC_SEQ) else '?' for p in positions]
    axes[0].set_xticks(range(n_pos))
    axes[0].set_xticklabels([f'{l}\n({aa})' for l, aa in zip(labels, pmdsc_aa)],
                             fontsize=8, ha='center')

    # Conservation bar plot
    conservation = info / np.log2(20)
    colors = ['#d62728' if pos in (highlight_positions or {}) else '#1f77b4'
              for pos in positions]
    axes[1].bar(range(n_pos), conservation, color=colors, alpha=0.8, edgecolor='k',
                linewidth=0.5)
    axes[1].set_ylabel('Conservation', fontsize=10)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xticks(range(n_pos))
    axes[1].set_xticklabels(labels, fontsize=8, rotation=45, ha='right')
    axes[1].axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    n_seqs = len(alignment_rows)
    axes[1].text(0.98, 0.95, f'n = {n_seqs} sequences', transform=axes[1].transAxes,
                 fontsize=9, ha='right', va='top',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")


def make_conservation_by_kingdom(alignment_rows, positions, labels):
    """Conservation bar chart split by kingdom."""
    kingdoms = sorted(set(r['kingdom'] for r in alignment_rows))
    n_pos = len(positions)
    aa_order = list('ACDEFGHIKLMNPQRSTVWY')

    fig, axes = plt.subplots(len(kingdoms), 1,
                              figsize=(max(12, n_pos * 0.7), 2.5 * len(kingdoms)),
                              sharex=True)
    if len(kingdoms) == 1:
        axes = [axes]

    for ax, kingdom in zip(axes, kingdoms):
        rows_k = [r for r in alignment_rows if r['kingdom'] == kingdom]
        if not rows_k:
            continue

        counts = np.zeros((n_pos, len(aa_order)))
        for row in rows_k:
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
        conservation = info / np.log2(20)

        # Most common AA at each position
        top_aa = []
        for i in range(n_pos):
            if counts[i].sum() > 0:
                top_idx = np.argmax(counts[i])
                top_aa.append(f"{aa_order[top_idx]} ({counts[i][top_idx]:.0f}/{totals[i][0]:.0f})")
            else:
                top_aa.append('-')

        color = KINGDOM_COLORS.get(kingdom, '#7f7f7f')
        ax.bar(range(n_pos), conservation, color=color, alpha=0.8,
               edgecolor='k', linewidth=0.5)
        ax.set_ylabel('Cons.', fontsize=9)
        ax.set_ylim(0, 1.1)
        ax.set_title(f'{kingdom} (n={len(rows_k)})', fontsize=10,
                     fontweight='bold', color=color)
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.3)

        for i, (cons, aa_str) in enumerate(zip(conservation, top_aa)):
            if cons > 0.3:
                ax.text(i, cons + 0.02, aa_str.split(' ')[0], ha='center',
                        fontsize=7, fontweight='bold')

    axes[-1].set_xticks(range(n_pos))
    pmdsc_aa = [PMDSC_SEQ[p - 1] if p <= len(PMDSC_SEQ) else '?' for p in positions]
    axes[-1].set_xticklabels([f'{l}\n({aa})' for l, aa in zip(labels, pmdsc_aa)],
                              fontsize=8, ha='center')

    fig.suptitle('Catalytic & Gateway Residue Conservation by Kingdom',
                 fontsize=13, fontweight='bold', y=1.01)
    fig.tight_layout()
    outpath = os.path.join(ANALYSIS_DIR, 'catalytic_conservation_by_kingdom.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")


def main():
    print("Loading sequences from Foldseek results...")
    sequences = load_all_sequences()
    print(f"Loaded {len(sequences)} unique sequences")

    kingdoms = Counter(v['kingdom'] for v in sequences.values())
    for k, v in kingdoms.most_common():
        print(f"  {k}: {v}")

    # Define positions to analyze
    cat_positions = sorted(CATALYTIC_RESIDUES.keys())
    gw_positions = sorted(GATEWAY_RESIDUES.keys())
    all_positions = sorted(set(cat_positions + gw_positions))

    cat_labels = [CATALYTIC_RESIDUES[p] for p in cat_positions]
    gw_labels = [GATEWAY_RESIDUES[p] for p in gw_positions]
    all_labels = []
    for p in all_positions:
        if p in CATALYTIC_RESIDUES:
            all_labels.append(CATALYTIC_RESIDUES[p])
        else:
            all_labels.append(GATEWAY_RESIDUES[p])

    # Also include context window around catalytic sites
    context_positions = set()
    for p in cat_positions:
        for dp in range(-CONTEXT_WINDOW, CONTEXT_WINDOW + 1):
            if 1 <= p + dp <= len(PMDSC_SEQ):
                context_positions.add(p + dp)
    context_positions = sorted(context_positions)
    context_labels = []
    for p in context_positions:
        if p in CATALYTIC_RESIDUES:
            context_labels.append(f"*{CATALYTIC_RESIDUES[p]}")
        elif p in GATEWAY_RESIDUES:
            context_labels.append(f"#{GATEWAY_RESIDUES[p]}")
        else:
            aa = PMDSC_SEQ[p - 1]
            context_labels.append(f"{aa}{p}")

    # Build alignment matrix
    print(f"\nAligning {len(sequences)} sequences to PMDsc...")
    alignment_rows = build_alignment_matrix(sequences, sorted(set(all_positions + context_positions)))
    print(f"Successfully aligned: {len(alignment_rows)} sequences")

    # 1. Logo: catalytic residues only
    print("\nGenerating logo plots...")
    make_logo(alignment_rows, cat_positions, cat_labels,
              'Catalytic Residue Conservation Across MDD/PMD Homologs',
              os.path.join(ANALYSIS_DIR, 'logo_catalytic_residues.png'))

    # 2. Logo: catalytic + gateway residues
    make_logo(alignment_rows, all_positions, all_labels,
              'Catalytic & Gateway Residue Conservation Across MDD/PMD Homologs',
              os.path.join(ANALYSIS_DIR, 'logo_catalytic_and_gateway.png'),
              highlight_positions=set(gw_positions))

    # 3. Logo: catalytic sites with context
    make_logo(alignment_rows, context_positions, context_labels,
              'Catalytic Site Conservation with ±2 Residue Context',
              os.path.join(ANALYSIS_DIR, 'logo_catalytic_with_context.png'),
              highlight_positions=set(cat_positions))

    # 4. Conservation by kingdom
    make_conservation_by_kingdom(alignment_rows, all_positions, all_labels)

    # Print conservation summary
    print("\n=== Conservation Summary ===")
    aa_order = list('ACDEFGHIKLMNPQRSTVWY')
    print(f"{'Position':<10} {'PMDsc':>6} {'Top AA':>7} {'Freq':>6} {'Cons':>6}")
    print("-" * 40)
    for pos, label in zip(all_positions, all_labels):
        pmdsc_aa = PMDSC_SEQ[pos - 1]
        counts = Counter()
        total = 0
        for row in alignment_rows:
            aa = row.get(pos, '-')
            if aa in aa_order:
                counts[aa] += 1
                total += 1
        if total > 0:
            top_aa, top_count = counts.most_common(1)[0]
            freq = top_count / total
            # Info content
            freqs = np.array([counts.get(a, 0) / total for a in aa_order])
            info = np.log2(20) + np.sum(freqs * np.where(freqs > 0, np.log2(freqs), 0))
            cons = info / np.log2(20)
            match = '=' if top_aa == pmdsc_aa else '!'
            print(f"{label:<10} {pmdsc_aa:>6} {top_aa:>5}{match} {freq:>5.0%} {cons:>5.1%}")

    print("\nDone!")


if __name__ == '__main__':
    main()
