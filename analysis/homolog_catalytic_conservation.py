#!/usr/bin/env python3
"""Re-filter homologs by proper catalytic vs binding vs gateway conservation.

Classifies residues into functional groups and checks which homologs
conserve the catalytic machinery while varying at gateway/entrance positions.
"""

import os
import csv
import yaml
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner

ANALYSIS_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.join(os.path.dirname(ANALYSIS_DIR), 'structures', 'homologs', 'cofold_configs')

PMDSC_SEQ = (
    "MTVYTASVTAPVNIATLKYWGKRDTKLNLPTNSSISVTLSQDDLRTLTSAATAPEFERDTLWLNGEPHSIDNERTQNCLRDLR"
    "QLRKEMESKDASLPTLSQWKLHIVSENNFPTAAGLASSAAGFAALVSAIAKLYQLPQSTSEISRIARKGSGSACRSLFGGYVAWEM"
    "GKAEDGHDSMAVQIADSSDWPQMKACVLVVSDIKKDVSSTQGMQLTVATSELFKERIEHVVPKRFEVMRKAIVEKDFATFAKETMM"
    "DSNSFHATCLDSFPPIFYMNDTSKRIISWCHTINQFYGETIVAYTFDAGPNAVLYYLAENESKLFAFIYKLFGSVPGWDKKFTTEQL"
    "EAFNHQFESSNFTARELDLELQKDVARVILTQVGSGPQETNESLIDAKTGLPKE"
)

CATALYTIC_CORE = {302: 'D302', 158: 'R158'}
ATP_BINDING = {121: 'S121', 155: 'S155', 208: 'S208'}
SUBSTRATE_ANCHORING = {120: 'S120', 153: 'S153'}
SUBSTRATE_BINDING = {18: 'K18', 19: 'Y19', 20: 'W20', 22: 'K22'}
GATEWAY = {74: 'R74', 147: 'R147', 212: 'M212'}
ENTRANCE_CONTACTS = {25: 'T25', 28: 'N28', 72: 'N72', 73: 'G73', 209: 'T209'}

ALL_GROUPS = [
    ('Catalytic core', CATALYTIC_CORE),
    ('ATP binding', ATP_BINDING),
    ('Substrate anchoring', SUBSTRATE_ANCHORING),
    ('Substrate binding', SUBSTRATE_BINDING),
    ('Gateway', GATEWAY),
    ('Entrance contacts', ENTRANCE_CONTACTS),
]

ALL_POSITIONS = sorted(set(
    list(CATALYTIC_CORE) + list(ATP_BINDING) + list(SUBSTRATE_ANCHORING) +
    list(SUBSTRATE_BINDING) + list(GATEWAY) + list(ENTRANCE_CONTACTS)
))


def get_sequence_from_yaml(yaml_path):
    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    for entry in config.get('sequences', []):
        if 'protein' in entry:
            return entry['protein']['sequence']
    return None


def align_and_extract(target_seq, positions):
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


def count_conserved(mapped_aas, group):
    n = 0
    for pos, label in group.items():
        pmdsc_aa = PMDSC_SEQ[pos - 1]
        if mapped_aas.get(pos, '-') == pmdsc_aa:
            n += 1
    return n


def get_position_label(pos):
    for _, group in ALL_GROUPS:
        if pos in group:
            return group[pos]
    return f'?{pos}'


def get_position_role(pos):
    for name, group in ALL_GROUPS:
        if pos in group:
            return name
    return 'Unknown'


def main():
    csv_path = os.path.join(ANALYSIS_DIR, 'homolog_entrance_analysis.csv')
    homologs = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            homologs.append(row)
    print(f"Loaded {len(homologs)} homologs from entrance analysis\n")

    results = []
    for h in homologs:
        name = h['name']
        yaml_path = os.path.join(CONFIG_DIR, f'{name}_2mvap.yaml')
        if not os.path.exists(yaml_path):
            continue
        seq = get_sequence_from_yaml(yaml_path)
        if not seq:
            continue
        try:
            mapped = align_and_extract(seq, ALL_POSITIONS)
        except Exception as e:
            print(f"  Skipping {name}: {e}")
            continue

        row = dict(h)
        row['cat_core'] = count_conserved(mapped, CATALYTIC_CORE)
        row['atp_bind'] = count_conserved(mapped, ATP_BINDING)
        row['sub_anchor'] = count_conserved(mapped, SUBSTRATE_ANCHORING)
        row['sub_bind'] = count_conserved(mapped, SUBSTRATE_BINDING)
        row['gateway_diff'] = len(GATEWAY) - count_conserved(mapped, GATEWAY)
        row['entrance_diff'] = len(ENTRANCE_CONTACTS) - count_conserved(mapped, ENTRANCE_CONTACTS)
        for pos in ALL_POSITIONS:
            row[f'pos_{pos}'] = mapped.get(pos, '-')
        results.append(row)

    print(f"Successfully aligned: {len(results)}\n")

    # Save CSV
    pos_cols = [f'pos_{p}' for p in ALL_POSITIONS]
    fieldnames = [
        'rank', 'name', 'clean_id', 'organism', 'kingdom', 'seq_identity',
        'frac_entrance', 'frac_active_site', 'frac_close',
        'cat_core', 'atp_bind', 'sub_anchor', 'sub_bind',
        'gateway_diff', 'entrance_diff',
    ] + pos_cols

    out_csv = os.path.join(ANALYSIS_DIR, 'homolog_catalytic_conservation.csv')
    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    print(f"Saved CSV: {out_csv}\n")

    # --- Console output ---

    # Filter: catalytic core conserved (D302 + R158)
    core_pass = [r for r in results if r['cat_core'] == 2]
    print(f"=== HOMOLOGS WITH CATALYTIC CORE CONSERVED (D302 + R158): {len(core_pass)}/{len(results)} ===\n")

    header = (f"{'Name':<25} {'Kingdom':<10} {'SeqID':>5} {'Close%':>6} "
              f"{'Cat':>3} {'ATP':>3} {'Anc':>3} {'Sub':>3} {'GwDiff':>6} "
              f"  R74  R147 M212   T25  N28  N72  G73 T209")
    print(header)
    print("-" * len(header))

    for r in sorted(core_pass, key=lambda x: float(x['frac_close'])):
        gw = f"  {r['pos_74']:>3}  {r['pos_147']:>3} {r['pos_212']:>3}"
        ent = f"   {r['pos_25']:>3}  {r['pos_28']:>3}  {r['pos_72']:>3}  {r['pos_73']:>3} {r['pos_209']:>3}"
        print(f"{r['name']:<25} {r['kingdom']:<10} {r['seq_identity']:>5} {r['frac_close']:>6} "
              f"{r['cat_core']:>1}/2 {r['atp_bind']:>1}/3 {r['sub_anchor']:>1}/2 {r['sub_bind']:>1}/4 "
              f"{r['gateway_diff']:>4}/3{gw}{ent}")

    # Further filter: core + ATP + anchoring
    full_machinery = [r for r in core_pass
                      if r['atp_bind'] >= 2 and r['sub_anchor'] >= 1]
    print(f"\n=== FURTHER FILTER: core + >=2/3 ATP + >=1/2 anchoring: {len(full_machinery)} ===\n")

    if full_machinery:
        print(header)
        print("-" * len(header))
        for r in sorted(full_machinery, key=lambda x: float(x['frac_close'])):
            gw = f"  {r['pos_74']:>3}  {r['pos_147']:>3} {r['pos_212']:>3}"
            ent = f"   {r['pos_25']:>3}  {r['pos_28']:>3}  {r['pos_72']:>3}  {r['pos_73']:>3} {r['pos_209']:>3}"
            print(f"{r['name']:<25} {r['kingdom']:<10} {r['seq_identity']:>5} {r['frac_close']:>6} "
                  f"{r['cat_core']:>1}/2 {r['atp_bind']:>1}/3 {r['sub_anchor']:>1}/2 {r['sub_bind']:>1}/4 "
                  f"{r['gateway_diff']:>4}/3{gw}{ent}")

    # --- Figure: heatmap of core-passing homologs ---
    if not core_pass:
        print("\nNo homologs pass catalytic core filter - no figure generated.")
        return

    fig_data = sorted(core_pass, key=lambda x: float(x['frac_close']))
    n_rows = len(fig_data)
    n_cols = len(ALL_POSITIONS)

    fig_height = max(4, 0.45 * n_rows + 2)
    fig, ax = plt.subplots(figsize=(max(14, n_cols * 0.8), fig_height))

    role_colors = {
        'Catalytic core': '#e74c3c',
        'ATP binding': '#3498db',
        'Substrate anchoring': '#9b59b6',
        'Substrate binding': '#1f77b4',
        'Gateway': '#ff7f0e',
        'Entrance contacts': '#2ecc71',
    }

    col_labels = [get_position_label(p) for p in ALL_POSITIONS]
    pmdsc_aas = [PMDSC_SEQ[p - 1] for p in ALL_POSITIONS]

    # PMDsc reference row
    for j, (pos, pmdsc_aa) in enumerate(zip(ALL_POSITIONS, pmdsc_aas)):
        role = get_position_role(pos)
        color = role_colors.get(role, '#7f7f7f')
        ax.text(j, n_rows, pmdsc_aa, ha='center', va='center', fontsize=11,
                fontweight='bold', color='white',
                bbox=dict(boxstyle='round,pad=0.25', facecolor=color, alpha=0.9))

    for i, r in enumerate(fig_data):
        y = n_rows - 1 - i
        for j, (pos, pmdsc_aa) in enumerate(zip(ALL_POSITIONS, pmdsc_aas)):
            aa = r.get(f'pos_{pos}', '-')
            if aa == pmdsc_aa:
                ax.text(j, y, aa, ha='center', va='center', fontsize=10,
                        fontweight='bold', color='#2ca02c',
                        bbox=dict(boxstyle='round,pad=0.25', facecolor='#d4edda', alpha=0.9))
            elif aa == '-':
                ax.text(j, y, '-', ha='center', va='center', fontsize=9, color='#999999')
            else:
                ax.text(j, y, aa, ha='center', va='center', fontsize=10,
                        fontweight='bold', color='#d62728',
                        bbox=dict(boxstyle='round,pad=0.25', facecolor='#f8d7da', alpha=0.9))

    # Annotations on right side
    for i, r in enumerate(fig_data):
        y = n_rows - 1 - i
        ax.text(n_cols + 0.3, y, f"close={r['frac_close']}", fontsize=7,
                va='center', ha='left', color='#555555')

    row_labels = ['PMDsc WT'] + [
        f"{r['name'].replace('homolog_','').replace('MDD_','').replace('MPD_','')} "
        f"({r['kingdom']}, {r['seq_identity']}%)"
        for r in fig_data
    ]

    ax.set_xlim(-0.7, n_cols + 1.5)
    ax.set_ylim(-0.7, n_rows + 0.7)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, fontsize=8, rotation=45, ha='right')
    ax.set_yticks(range(n_rows + 1))
    ax.set_yticklabels(row_labels[::-1], fontsize=7, ha='right')

    for j, pos in enumerate(ALL_POSITIONS):
        role = get_position_role(pos)
        color = role_colors.get(role, '#7f7f7f')
        ax.axvspan(j - 0.45, j + 0.45, alpha=0.04, color=color, zorder=0)

    ax.axhline(y=n_rows - 0.5, color='black', linewidth=2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, bottom=False)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=r) for r, c in role_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=7, ncol=2,
              title='Residue group', title_fontsize=8)

    ax.set_title('Homologs with Catalytic Core Conserved (D302 + R158)\n'
                 'Green = matches PMDsc, Red = differs. Sorted by close binding fraction.',
                 fontsize=12, fontweight='bold', pad=15)

    fig.tight_layout()
    fig_path = os.path.join(ANALYSIS_DIR, 'homolog_catalytic_filter.png')
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved figure: {fig_path}")


if __name__ == '__main__':
    main()
