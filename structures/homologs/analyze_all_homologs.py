#!/usr/bin/env python3
"""Analyze entrance fractions for all cofolded homologs.

Computes L2-L1 distance, binding mode classification, and entrance
fraction for each homolog. Outputs ranked CSV and summary figures.
"""

import os
import csv
import glob
import json
import yaml
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'output_2mvap')
CONFIG_DIR = os.path.join(SCRIPT_DIR, 'cofold_configs')
ANALYSIS_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'

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

# Substrate binding residues in PMDsc
SUBSTRATE_BINDING = {18, 19, 20, 22}

CLOSE_THRESHOLD = 15.0
CATALYTIC_CONTACT_THRESHOLD = 5

TAXONOMY = {
    'Saccharomyces': 'Fungi', 'Eremothecium': 'Fungi', 'Debaryomyces': 'Fungi',
    'Candida': 'Fungi', 'Fusarium': 'Fungi', 'Schizosaccharomyces': 'Fungi',
    'Aspergillus': 'Fungi', 'Ganoderma': 'Fungi', 'Cryptococcus': 'Fungi',
    'Mus': 'Animal', 'Bos': 'Animal', 'Rattus': 'Animal', 'Danio': 'Animal',
    'Homo': 'Animal', 'Trypanosoma': 'Protozoa', 'Dictyostelium': 'Protozoa',
    'Arabidopsis': 'Plant', 'Panax': 'Plant',
    'Saccharolobus': 'Archaea', 'Sulfolobus': 'Archaea', 'Thermoplasma': 'Archaea',
    'Picrophilus': 'Archaea', 'Haloferax': 'Archaea',
    'Staphylococcus': 'Bacteria', 'Enterococcus': 'Bacteria', 'Streptococcus': 'Bacteria',
    'Bacillus': 'Bacteria', 'Listeria': 'Bacteria', 'Anaerolinea': 'Bacteria',
    'Legionella': 'Bacteria', 'Synechococcus': 'Bacteria', 'Prochlorococcus': 'Bacteria',
    'Clostridium': 'Bacteria', 'Rhodococcus': 'Bacteria', 'Corynebacterium': 'Bacteria',
    'Geobacillus': 'Bacteria', 'Halalkalibacterium': 'Bacteria',
    'Roseiflexus': 'Bacteria', 'Heliomicrobium': 'Bacteria',
    'Nocardioides': 'Bacteria', 'Parvibaculum': 'Bacteria',
    'Phocaeicola': 'Bacteria', 'Rhizobium': 'Bacteria',
    'Microcystis': 'Bacteria', 'Picosynechococcus': 'Bacteria',
    'Parasynechococcus': 'Bacteria', 'Limosilactobacillus': 'Bacteria',
    'Lachnoclostridium': 'Bacteria', 'Desulforamulus': 'Bacteria',
    'Thermoanaerobacter': 'Bacteria', 'Syntrophomonas': 'Bacteria',
}


def classify_organism(org):
    for genus, k in TAXONOMY.items():
        if genus.lower() in org.lower():
            return k
    return 'Unknown'


def get_catalytic_resnums(yaml_path):
    """Extract catalytic residue numbers from YAML config (L1 pocket contacts)."""
    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    resnums = set()
    for constraint in config.get('constraints', []):
        if 'pocket' in constraint and constraint['pocket'].get('binder') == 'L':
            for contact in constraint['pocket'].get('contacts', []):
                resnums.add(contact[1])
        if 'contact' in constraint:
            t2 = constraint['contact'].get('token2', [])
            if len(t2) >= 2 and t2[0] == 'A':
                resnums.add(t2[1])
    return resnums


def get_sequence_from_yaml(yaml_path):
    """Extract protein sequence from YAML config."""
    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    for entry in config.get('sequences', []):
        if 'protein' in entry:
            return entry['protein']['sequence']
    return None


def check_substrate_binding_conservation(sequence):
    """Check if K18/Y19/W20/K22 equivalents are conserved."""
    try:
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        alignment = aligner.align(PMDSC_SEQ, sequence)[0]
        aln_str = str(alignment).split('\n')
        aln_q, aln_t = aln_str[0], aln_str[2]

        conserved = 0
        mapping = {}
        q_pos, t_pos = 0, 0
        for i in range(min(len(aln_q), len(aln_t))):
            q_gap = aln_q[i] in ('-', '.')
            t_gap = aln_t[i] in ('-', '.')
            if not q_gap and not t_gap:
                q_resnum = q_pos + 1
                if q_resnum in SUBSTRATE_BINDING:
                    pmdsc_aa = PMDSC_SEQ[q_pos]
                    target_aa = sequence[t_pos]
                    if pmdsc_aa == target_aa:
                        conserved += 1
                    mapping[q_resnum] = (pmdsc_aa, target_aa)
                q_pos += 1
                t_pos += 1
            elif q_gap:
                t_pos += 1
            elif t_gap:
                q_pos += 1

        return conserved, len(SUBSTRATE_BINDING), mapping
    except Exception:
        return 0, 4, {}


def analyze_homolog(name):
    """Compute entrance fraction for a single homolog."""
    # Find prediction directory
    pred_dir = None
    for subdir in os.listdir(OUTPUT_DIR):
        if subdir.startswith(name) and os.path.isdir(os.path.join(OUTPUT_DIR, subdir)):
            inner = os.path.join(OUTPUT_DIR, subdir, f'boltz_results_{subdir}', 'predictions', subdir)
            if os.path.isdir(inner):
                pred_dir = inner
                break

    if not pred_dir:
        return None

    # Get catalytic residues from YAML
    yaml_path = os.path.join(CONFIG_DIR, f'{name}_2mvap.yaml')
    if not os.path.exists(yaml_path):
        # Try alternative naming
        for f in os.listdir(CONFIG_DIR):
            if name.replace('homolog_', '') in f:
                yaml_path = os.path.join(CONFIG_DIR, f)
                break

    if os.path.exists(yaml_path):
        cat_res = get_catalytic_resnums(yaml_path)
    else:
        cat_res = set()

    pdb_files = sorted(glob.glob(os.path.join(pred_dir, '*.pdb')))
    if not pdb_files:
        return None

    modes, dists = [], []
    for pdb_path in pdb_files:
        chains = {}
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith(('ATOM', 'HETATM')):
                    continue
                if line[76:78].strip() in ('H', 'D'):
                    continue
                c = line[21]
                chains.setdefault(c, []).append({
                    'resnum': int(line[22:26]),
                    'xyz': np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                })

        if 'A' not in chains or 'L' not in chains or 'S' not in chains:
            continue

        l1_com = np.mean([a['xyz'] for a in chains['L']], axis=0)
        s_com = np.mean([a['xyz'] for a in chains['S']], axis=0)
        dist = np.linalg.norm(s_com - l1_com)
        dists.append(dist)

        s_coords = np.array([a['xyz'] for a in chains['S']])
        n_cat = 0
        if cat_res:
            n_cat = sum(1 for a in chains['A']
                        if a['resnum'] in cat_res
                        and np.min(np.linalg.norm(s_coords - a['xyz'], axis=1)) < 4.0)

        if dist > 25:
            modes.append('distant')
        elif dist > CLOSE_THRESHOLD:
            modes.append('intermediate')
        elif n_cat >= CATALYTIC_CONTACT_THRESHOLD:
            modes.append('active_site')
        else:
            modes.append('entrance')

    if not modes:
        return None

    n = len(modes)
    return {
        'n_models': n,
        'mean_dist': np.mean(dists),
        'median_dist': np.median(dists),
        'frac_entrance': sum(1 for m in modes if m == 'entrance') / n,
        'frac_active_site': sum(1 for m in modes if m == 'active_site') / n,
        'frac_close': sum(1 for m in modes if m in ('entrance', 'active_site')) / n,
    }


def main():
    print("=== Analyzing All Homolog Cofolding Results ===\n")

    # Find all completed jobs
    completed = []
    if not os.path.isdir(OUTPUT_DIR):
        print(f"No output directory found: {OUTPUT_DIR}")
        return

    for subdir in sorted(os.listdir(OUTPUT_DIR)):
        if subdir.endswith('.log') or not os.path.isdir(os.path.join(OUTPUT_DIR, subdir)):
            continue
        inner = os.path.join(OUTPUT_DIR, subdir, f'boltz_results_{subdir}', 'predictions', subdir)
        if os.path.isdir(inner):
            pdb_count = len(glob.glob(os.path.join(inner, '*.pdb')))
            if pdb_count >= 10:
                completed.append(subdir)

    print(f"Completed jobs: {len(completed)}")

    # Load manifest for metadata
    manifest = {}
    manifest_path = os.path.join(SCRIPT_DIR, 'batch_manifest.csv')
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            for row in csv.DictReader(f):
                manifest[row['safe_name']] = row

    # Analyze each
    results = []
    for subdir in completed:
        name = subdir.replace('_2mvap', '')
        r = analyze_homolog(name)
        if not r:
            continue

        # Get metadata
        meta = manifest.get(name, {})
        organism = meta.get('organism', '')
        kingdom = meta.get('kingdom', '')
        seq_identity = meta.get('seq_identity', '')
        clean_id = meta.get('clean_id', name.replace('homolog_', '').replace('MDD_', '').replace('MPD_', ''))

        # If not in manifest, try to get from YAML
        if not organism:
            yaml_path = os.path.join(CONFIG_DIR, f'{name}_2mvap.yaml')
            if os.path.exists(yaml_path):
                seq = get_sequence_from_yaml(yaml_path)
                if seq:
                    kingdom = 'Unknown'

        # For original 10, fill in manually
        original_meta = {
            'MDD_Spombe': ('O13963', 'Schizosaccharomyces pombe', 'Fungi', '54.9'),
            'MDD_Cneo': ('6XR5', 'Cryptococcus neoformans', 'Fungi', '50.7'),
            'MDD_Hsap': ('3D4J', 'Homo sapiens', 'Animal', '44.4'),
            'MDD_Ddis': ('Q54YQ9', 'Dictyostelium discoideum', 'Protozoa', '44.0'),
            'MDD_Sepi': ('4DPU', 'Staphylococcus epidermidis', 'Bacteria', '27.6'),
            'MDD_Saur': ('2HK2', 'Staphylococcus aureus', 'Bacteria', '27.5'),
            'MDD_Efae': ('6E2U', 'Enterococcus faecalis', 'Bacteria', '28.0'),
            'MDD_Spyo': ('2GS8', 'Streptococcus pyogenes', 'Bacteria', '23.4'),
            'MPD_Ather': ('6N0X', 'Anaerolinea thermophila', 'Bacteria', '31.7'),
            'MDD_Ssol': ('5GMD', 'Saccharolobus solfataricus', 'Archaea', '32.8'),
        }
        if name in original_meta:
            clean_id, organism, kingdom, seq_identity = original_meta[name]

        if not kingdom and organism:
            kingdom = classify_organism(organism)

        # Check substrate binding conservation
        yaml_path = os.path.join(CONFIG_DIR, f'{name}_2mvap.yaml')
        sub_conserved, sub_total, sub_mapping = 0, 4, {}
        if os.path.exists(yaml_path):
            seq = get_sequence_from_yaml(yaml_path)
            if seq:
                sub_conserved, sub_total, sub_mapping = check_substrate_binding_conservation(seq)

        results.append({
            'name': name,
            'clean_id': clean_id,
            'organism': organism,
            'kingdom': kingdom,
            'seq_identity': seq_identity,
            'n_models': r['n_models'],
            'mean_dist': r['mean_dist'],
            'median_dist': r['median_dist'],
            'frac_entrance': r['frac_entrance'],
            'frac_active_site': r['frac_active_site'],
            'frac_close': r['frac_close'],
            'substrate_binding_conserved': sub_conserved,
            'substrate_binding_total': sub_total,
        })

    # Sort by frac_close (ascending = best Ki predicted)
    results.sort(key=lambda x: x['frac_close'])

    # Add PMDsc references
    pmdsc_refs = [
        {'name': 'PMDsc_HKQ', 'clean_id': 'HKQ', 'organism': 'S. cerevisiae (engineered)',
         'kingdom': 'Fungi', 'seq_identity': '100', 'n_models': 50,
         'mean_dist': 33.5, 'median_dist': 36.1, 'frac_entrance': 0.02,
         'frac_active_site': 0.0, 'frac_close': 0.02, 'substrate_binding_conserved': 4, 'substrate_binding_total': 4},
        {'name': 'PMDsc_WT', 'clean_id': 'WT', 'organism': 'S. cerevisiae',
         'kingdom': 'Fungi', 'seq_identity': '100', 'n_models': 50,
         'mean_dist': 26.3, 'median_dist': 32.4, 'frac_entrance': 0.28,
         'frac_active_site': 0.04, 'frac_close': 0.32, 'substrate_binding_conserved': 4, 'substrate_binding_total': 4},
    ]

    # Save CSV
    csv_path = os.path.join(ANALYSIS_DIR, 'homolog_entrance_analysis.csv')
    fieldnames = ['rank', 'name', 'clean_id', 'organism', 'kingdom', 'seq_identity',
                  'n_models', 'mean_dist', 'median_dist', 'frac_entrance',
                  'frac_active_site', 'frac_close', 'substrate_binding_conserved']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, r in enumerate(results):
            writer.writerow({
                'rank': i + 1,
                'name': r['name'],
                'clean_id': r['clean_id'],
                'organism': r['organism'],
                'kingdom': r['kingdom'],
                'seq_identity': r['seq_identity'],
                'n_models': r['n_models'],
                'mean_dist': f"{r['mean_dist']:.1f}",
                'median_dist': f"{r['median_dist']:.1f}",
                'frac_entrance': f"{r['frac_entrance']:.2f}",
                'frac_active_site': f"{r['frac_active_site']:.2f}",
                'frac_close': f"{r['frac_close']:.2f}",
                'substrate_binding_conserved': f"{r['substrate_binding_conserved']}/{r['substrate_binding_total']}",
            })
    print(f"\nSaved results to {csv_path}")

    # Print top and bottom
    print(f"\n=== TOP 20 (lowest substrate inhibition predicted) ===")
    print(f"{'Rank':<5} {'Name':<25} {'Kingdom':<10} {'Organism':<30} {'SeqID':>6} {'Ent%':>5} {'Act%':>5} {'Close%':>6} {'SubBind':>8}")
    print("-" * 110)
    for i, r in enumerate(results[:20]):
        print(f"{i+1:<5} {r['name']:<25} {r['kingdom']:<10} {r['organism'][:28]:<30} "
              f"{r['seq_identity']:>5}% {r['frac_entrance']:>4.0%} {r['frac_active_site']:>4.0%} "
              f"{r['frac_close']:>5.0%} {r['substrate_binding_conserved']}/{r['substrate_binding_total']}")

    # Reference lines
    print(f"\n--- PMDsc References ---")
    for ref in pmdsc_refs:
        print(f"{'ref':<5} {ref['name']:<25} {ref['kingdom']:<10} {ref['organism'][:28]:<30} "
              f"{ref['seq_identity']:>5}% {ref['frac_entrance']:>4.0%} {ref['frac_active_site']:>4.0%} "
              f"{ref['frac_close']:>5.0%} {ref['substrate_binding_conserved']}/{ref['substrate_binding_total']}")

    print(f"\n=== BOTTOM 10 (highest substrate inhibition predicted) ===")
    for i, r in enumerate(results[-10:]):
        idx = len(results) - 10 + i + 1
        print(f"{idx:<5} {r['name']:<25} {r['kingdom']:<10} {r['organism'][:28]:<30} "
              f"{r['seq_identity']:>5}% {r['frac_entrance']:>4.0%} {r['frac_active_site']:>4.0%} "
              f"{r['frac_close']:>5.0%} {r['substrate_binding_conserved']}/{r['substrate_binding_total']}")

    # Generate figure
    print("\nGenerating summary figure...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # Left: bar chart of close% ranked
    names = [r['name'].replace('homolog_', '').replace('MDD_', '').replace('MPD_', '') for r in results]
    close_fracs = [r['frac_close'] for r in results]
    kingdoms_list = [r['kingdom'] for r in results]
    kingdom_colors = {
        'Fungi': '#2ca02c', 'Animal': '#1f77b4', 'Plant': '#8c564b',
        'Protozoa': '#9467bd', 'Bacteria': '#d62728', 'Archaea': '#ff7f0e', 'Unknown': '#7f7f7f',
    }
    colors = [kingdom_colors.get(k, '#7f7f7f') for k in kingdoms_list]

    axes[0].barh(range(len(results)), close_fracs, color=colors, alpha=0.8, edgecolor='k', linewidth=0.3)
    axes[0].axvline(x=0.02, color='green', linestyle='--', alpha=0.7, label='HKQ (Ki=80)')
    axes[0].axvline(x=0.32, color='orange', linestyle='--', alpha=0.7, label='WT (Ki=18)')
    axes[0].set_xlabel('Close Binding Fraction (entrance + active_site)', fontsize=11)
    axes[0].set_ylabel('Homolog (ranked)', fontsize=11)
    axes[0].set_title('Predicted Substrate Inhibition\n(lower = less inhibition = higher Ki)', fontsize=12, fontweight='bold')
    axes[0].legend(fontsize=9)
    axes[0].invert_yaxis()
    if len(results) > 30:
        axes[0].set_yticks([])
    else:
        axes[0].set_yticks(range(len(results)))
        axes[0].set_yticklabels(names, fontsize=6)

    # Right: kingdom distribution of top vs bottom
    n_half = len(results) // 2
    top_kingdoms = Counter(r['kingdom'] for r in results[:n_half])
    bot_kingdoms = Counter(r['kingdom'] for r in results[n_half:])
    all_kingdoms = sorted(set(list(top_kingdoms.keys()) + list(bot_kingdoms.keys())))

    x = np.arange(len(all_kingdoms))
    w = 0.35
    axes[1].bar(x - w/2, [top_kingdoms.get(k, 0) for k in all_kingdoms], w,
                label='Low inhibition (top half)', color='#2ca02c', alpha=0.8)
    axes[1].bar(x + w/2, [bot_kingdoms.get(k, 0) for k in all_kingdoms], w,
                label='High inhibition (bottom half)', color='#d62728', alpha=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(all_kingdoms, rotation=45, ha='right')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Kingdom Distribution:\nLow vs High Substrate Inhibition', fontsize=12, fontweight='bold')
    axes[1].legend()

    fig.tight_layout()
    fig_path = os.path.join(ANALYSIS_DIR, 'homolog_substrate_inhibition_ranking.png')
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved figure: {fig_path}")

    print(f"\n=== SUMMARY ===")
    print(f"Total homologs analyzed: {len(results)}")
    better_than_hkq = sum(1 for r in results if r['frac_close'] <= 0.02)
    better_than_wt = sum(1 for r in results if r['frac_close'] < 0.32)
    print(f"Better than HKQ (close ≤ 2%): {better_than_hkq}")
    print(f"Better than WT (close < 32%): {better_than_wt}")

    # Kingdom breakdown of top candidates
    top_10 = results[:10]
    top_kingdoms = Counter(r['kingdom'] for r in top_10)
    print(f"\nTop 10 kingdom breakdown: {dict(top_kingdoms.most_common())}")


if __name__ == '__main__':
    main()
