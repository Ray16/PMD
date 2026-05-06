#!/usr/bin/env python3
"""Generate alignment comparison of top homolog candidates vs PMDsc.

Shows residue identity at all functionally important positions for
candidates with low entrance fraction and retained substrate binding.
"""

import os
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner

ANALYSIS_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'
CONFIG_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/cofold_configs'

PMDSC_SEQ = (
    "MTVYTASVTAPVNIATLKYWGKRDTKLNLPTNSSISVTLSQDDLRTLTSAATAPEFERDTLWLNGEPHSIDNERTQNCLRDLR"
    "QLRKEMESKDASLPTLSQWKLHIVSENNFPTAAGLASSAAGFAALVSAIAKLYQLPQSTSEISRIARKGSGSACRSLFGGYVAWEM"
    "GKAEDGHDSMAVQIADSSDWPQMKACVLVVSDIKKDVSSTQGMQLTVATSELFKERIEHVVPKRFEVMRKAIVEKDFATFAKETMM"
    "DSNSFHATCLDSFPPIFYMNDTSKRIISWCHTINQFYGETIVAYTFDAGPNAVLYYLAENESKLFAFIYKLFGSVPGWDKKFTTEQL"
    "EAFNHQFESSNFTARELDLELQKDVARVILTQVGSGPQETNESLIDAKTGLPKE"
)

POSITIONS = {
    18: ('K18', 'substrate binding'),
    19: ('Y19', 'substrate binding'),
    20: ('W20', 'substrate binding'),
    22: ('K22', 'substrate binding'),
    25: ('T25', 'entrance contact'),
    28: ('N28', 'entrance contact'),
    72: ('N72', 'entrance contact'),
    73: ('G73', 'entrance contact'),
    74: ('R74', 'gateway'),
    120: ('S120', 'catalytic'),
    121: ('S121', 'catalytic'),
    145: ('I145', 'structural'),
    147: ('R147', 'gateway'),
    153: ('S153', 'catalytic'),
    155: ('S155', 'catalytic'),
    158: ('R158', 'catalytic'),
    208: ('S208', 'catalytic'),
    209: ('T209', 'entrance contact'),
    212: ('M212', 'gateway'),
    302: ('D302', 'catalytic'),
}

CANDIDATES = [
    {'name': 'Q751D8', 'config': 'homolog_Q751D8_2mvap.yaml',
     'organism': 'Eremothecium gossypii', 'kingdom': 'Fungi',
     'seq_id': 68.8, 'close': 0.06, 'sub_bind': '4/4'},
    {'name': '3QT6', 'config': 'homolog_3QT6_2mvap.yaml',
     'organism': 'Staphylococcus epidermidis', 'kingdom': 'Bacteria',
     'seq_id': 28.9, 'close': 0.02, 'sub_bind': '4/4'},
    {'name': 'MDD_Ddis', 'config': 'MDD_Ddis_2mvap.yaml',
     'organism': 'Dictyostelium discoideum', 'kingdom': 'Protozoa',
     'seq_id': 44.0, 'close': 0.00, 'sub_bind': '3/4'},
    {'name': 'Q4WNV9', 'config': 'homolog_Q4WNV9_2mvap.yaml',
     'organism': 'Aspergillus fumigatus', 'kingdom': 'Fungi',
     'seq_id': 55.4, 'close': 0.10, 'sub_bind': '1/4'},
    {'name': 'Q6BY07', 'config': 'homolog_Q6BY07_2mvap.yaml',
     'organism': 'Debaryomyces hansenii', 'kingdom': 'Fungi',
     'seq_id': 64.8, 'close': 0.40, 'sub_bind': '4/4'},
    {'name': 'G9BIY1', 'config': 'homolog_G9BIY1_2mvap.yaml',
     'organism': 'Ganoderma lucidum', 'kingdom': 'Fungi',
     'seq_id': 49.6, 'close': 0.66, 'sub_bind': '4/4'},
]

ROLE_COLORS = {
    'substrate binding': '#1f77b4',
    'entrance contact': '#d62728',
    'gateway': '#ff7f0e',
    'catalytic': '#2ca02c',
    'structural': '#9467bd',
}


def get_sequence_from_yaml(yaml_path):
    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    for entry in config.get('sequences', []):
        if 'protein' in entry:
            return entry['protein']['sequence']
    return None


def align_positions(target_seq, positions):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(PMDSC_SEQ, target_seq)[0]
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


def main():
    sorted_pos = sorted(POSITIONS.keys())
    labels = [POSITIONS[p][0] for p in sorted_pos]
    roles = [POSITIONS[p][1] for p in sorted_pos]
    pmdsc_aas = [PMDSC_SEQ[p - 1] for p in sorted_pos]

    all_data = []
    row_labels = ['PMDsc WT']
    all_data.append(pmdsc_aas)

    for cand in CANDIDATES:
        yaml_path = os.path.join(CONFIG_DIR, cand['config'])
        seq = get_sequence_from_yaml(yaml_path)
        if not seq:
            continue
        mapped = align_positions(seq, sorted_pos)
        aas = [mapped.get(p, '-') for p in sorted_pos]
        all_data.append(aas)
        row_labels.append(f"{cand['name']} ({cand['organism']}, {cand['seq_id']}%)\n"
                          f"close={cand['close']:.0%}, sub={cand['sub_bind']}")

    n_rows = len(all_data)
    n_cols = len(sorted_pos)

    fig, ax = plt.subplots(figsize=(max(16, n_cols * 0.75), 1.2 + n_rows * 0.9))

    for i, (row_aas, row_label) in enumerate(zip(all_data, row_labels)):
        y = n_rows - 1 - i
        for j, (aa, pmdsc_aa) in enumerate(zip(row_aas, pmdsc_aas)):
            if i == 0:
                color = ROLE_COLORS.get(roles[j], '#7f7f7f')
                ax.text(j, y, aa, ha='center', va='center', fontsize=14,
                        fontweight='bold', color='white',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor=color, alpha=0.9))
            else:
                if aa == pmdsc_aa:
                    ax.text(j, y, aa, ha='center', va='center', fontsize=14,
                            fontweight='bold', color='#2ca02c',
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='#d4edda', alpha=0.9))
                elif aa == '-':
                    ax.text(j, y, '-', ha='center', va='center', fontsize=12,
                            color='#999999')
                else:
                    ax.text(j, y, aa, ha='center', va='center', fontsize=14,
                            fontweight='bold', color='#d62728',
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8d7da', alpha=0.9))

    ax.set_xlim(-1.5, n_cols - 0.5)
    ax.set_ylim(-0.8, n_rows - 0.2)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(labels, fontsize=9, rotation=45, ha='right')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels[::-1], fontsize=9, ha='right')

    for j, role in enumerate(roles):
        color = ROLE_COLORS.get(role, '#7f7f7f')
        ax.axvspan(j - 0.45, j + 0.45, alpha=0.05, color=color, zorder=0)

    ax.axhline(y=n_rows - 1.5, color='black', linewidth=2, linestyle='-')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, bottom=False)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=r) for r, c in ROLE_COLORS.items()]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=2,
              title='Residue role', title_fontsize=9)

    ax.set_title('Top Homolog Candidates vs PMDsc: Functional Residue Comparison\n'
                 'Green = matches PMDsc, Red = differs, Top row = PMDsc reference',
                 fontsize=13, fontweight='bold', pad=15)

    fig.tight_layout()
    outpath = os.path.join(ANALYSIS_DIR, 'logo_top_candidates_vs_pmdsc.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {outpath}")

    # Print summary table
    print(f"\n{'':>45} " + ' '.join(f'{l:>5}' for l in labels))
    print(f"{'':>45} " + ' '.join(f'{PMDSC_SEQ[p-1]:>5}' for p in sorted_pos))
    print('-' * (50 + 6 * n_cols))
    for cand_data, cand_label in zip(all_data[1:], row_labels[1:]):
        short_label = cand_label.split('\n')[0][:43]
        matches = sum(1 for a, b in zip(cand_data, pmdsc_aas) if a == b and a != '-')
        total_mapped = sum(1 for a in cand_data if a != '-')
        line = ' '.join(f'{a:>5}' for a in cand_data)
        print(f"{short_label:<45} {line}  ({matches}/{total_mapped} match)")


if __name__ == '__main__':
    main()
