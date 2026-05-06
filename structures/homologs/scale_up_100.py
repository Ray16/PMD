#!/usr/bin/env python3
"""Scale up homolog cofolding to 100 total.

Selects 90 additional diverse homologs beyond the original 10,
generates YAML configs, and creates a batch run script.
"""

import os
import csv
import yaml
import numpy as np
from Bio.Align import PairwiseAligner
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FOLDSEEK_DIR = os.path.join(SCRIPT_DIR, 'foldseek_results')
CONFIG_DIR = os.path.join(SCRIPT_DIR, 'cofold_configs')
SEQ_DIR = os.path.join(SCRIPT_DIR, 'sequences')

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
L1_POCKET_RESIDUES = [18, 19, 20, 22, 120, 153, 155, 158, 208, 302]
ATP_POCKET_RESIDUES = [121, 153, 155, 208]
MG_CONTACT_RESIDUE = 302
MVAP_SMILES = 'C[C@](O)(CCOP(=O)(O)O)CC(=O)O'

ALREADY_DONE = {
    'O13963', '6XR5', '3D4J', 'Q54YQ9', '4DPU', '2HK2', '6E2U', '2GS8', '6N0X', '5GMD',
}

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


def load_all_candidates():
    """Load all sequences from Foldseek m8 results."""
    candidates = []
    seen_ids = set()

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
                target_desc = fields[1].split(' ', 1)[1] if ' ' in fields[1] else ''
                seq_identity = float(fields[2])
                aln_len = int(fields[3])
                prob = float(fields[10])
                evalue = float(fields[11])
                bit_score = float(fields[12])
                t_len = int(fields[14])
                t_full_seq = fields[18].replace('-', '') if len(fields) > 18 else ''
                organism = fields[20] if len(fields) > 20 else ''

                if not t_full_seq or len(t_full_seq) < 150:
                    continue
                if aln_len < 200:
                    continue
                if seq_identity > 95:
                    continue

                # Clean ID
                if target_raw.startswith('AF-'):
                    clean_id = target_raw.split('-')[1]
                    db = 'alphafold'
                else:
                    clean_id = target_raw.split('-')[0].upper()
                    db = 'pdb'

                if clean_id in seen_ids:
                    continue
                seen_ids.add(clean_id)

                kingdom = classify_organism(organism)

                candidates.append({
                    'clean_id': clean_id,
                    'db': db,
                    'desc': target_desc,
                    'organism': organism,
                    'kingdom': kingdom,
                    'seq_identity': seq_identity,
                    'aln_len': aln_len,
                    'prob': prob,
                    'evalue': evalue,
                    'bit_score': bit_score,
                    't_len': t_len,
                    'sequence': t_full_seq,
                })

    # Sort by bit_score descending
    candidates.sort(key=lambda x: x['bit_score'], reverse=True)
    return candidates


def select_90_diverse(candidates):
    """Select 90 diverse homologs, excluding already-done ones."""
    selected = []
    # Track organism-level diversity (avoid too many strains of same species)
    species_count = Counter()

    for c in candidates:
        if c['clean_id'] in ALREADY_DONE:
            continue
        if len(selected) >= 90:
            break

        # Limit per-species redundancy (max 2 per species for PDB, 1 for AF)
        species = c['organism'].split(' ')[0] + ' ' + (c['organism'].split(' ')[1] if len(c['organism'].split(' ')) > 1 else '')
        max_per_species = 2 if c['db'] == 'pdb' else 1
        if species_count[species] >= max_per_species:
            continue

        species_count[species] += 1
        selected.append(c)

    return selected


def align_and_map(target_seq):
    """Align target to PMDsc and map catalytic residues."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(PMDSC_SEQ, target_seq)[0]

    # Use indices array: shape (2, aln_len)
    # Row 0 = target (PMDsc) indices, Row 1 = query (target_seq) indices
    # -1 means gap
    indices = alignment.indices
    mapping = {}
    for col in range(indices.shape[1]):
        q_idx = indices[0, col]  # PMDsc position (0-indexed)
        t_idx = indices[1, col]  # target position (0-indexed)
        if q_idx >= 0 and t_idx >= 0:
            q_resnum = q_idx + 1  # 1-indexed
            t_resnum = t_idx + 1
            if q_resnum in CATALYTIC_RESIDUES:
                mapping[q_resnum] = t_resnum

    return mapping


def generate_yaml(name, sequence, residue_mapping):
    """Generate Boltz2 YAML config."""
    l1_contacts = [['A', int(residue_mapping[r])] for r in L1_POCKET_RESIDUES if r in residue_mapping]
    atp_contacts = [['A', int(residue_mapping[r])] for r in ATP_POCKET_RESIDUES if r in residue_mapping]
    mg_res = int(residue_mapping[MG_CONTACT_RESIDUE]) if MG_CONTACT_RESIDUE in residue_mapping else None

    config = {
        'version': 1,
        'sequences': [
            {'protein': {'id': 'A', 'sequence': sequence}},
            {'ligand': {'id': 'L', 'smiles': MVAP_SMILES}},
            {'ligand': {'id': 'S', 'smiles': MVAP_SMILES}},
            {'ligand': {'id': 'T', 'ccd': 'ATP'}},
            {'ligand': {'id': 'M', 'ccd': 'MG'}},
        ],
        'constraints': [],
    }
    if l1_contacts:
        config['constraints'].append({'pocket': {'binder': 'L', 'contacts': l1_contacts, 'max_distance': 6.0}})
    if atp_contacts:
        config['constraints'].append({'pocket': {'binder': 'T', 'contacts': atp_contacts, 'max_distance': 6.0}})
    if mg_res:
        config['constraints'].append({'contact': {'token1': ['M', 'MG'], 'token2': ['A', mg_res], 'max_distance': 4.0}})

    return config


def main():
    os.makedirs(CONFIG_DIR, exist_ok=True)
    os.makedirs(SEQ_DIR, exist_ok=True)

    print("Loading all candidates from Foldseek results...")
    candidates = load_all_candidates()
    print(f"Total unique candidates: {len(candidates)}")

    print(f"\nSelecting 90 diverse homologs (excluding {len(ALREADY_DONE)} already done)...")
    selected = select_90_diverse(candidates)
    print(f"Selected: {len(selected)}")

    kingdoms = Counter(c['kingdom'] for c in selected)
    print(f"Kingdom distribution: {dict(kingdoms.most_common())}")

    # Process each
    successful = []
    failed = []
    for i, c in enumerate(selected):
        cid = c['clean_id']
        safe_name = f"homolog_{cid}"

        try:
            mapping = align_and_map(c['sequence'])
        except Exception as e:
            failed.append((cid, f"alignment error: {e}"))
            continue

        if len(mapping) < 5:
            failed.append((cid, f"only {len(mapping)} residues mapped"))
            continue

        config = generate_yaml(safe_name, c['sequence'], mapping)
        yaml_path = os.path.join(CONFIG_DIR, f'{safe_name}_2mvap.yaml')
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)

        # Save FASTA
        fasta_path = os.path.join(SEQ_DIR, f'{cid}.fasta')
        if not os.path.exists(fasta_path):
            with open(fasta_path, 'w') as f:
                f.write(f">{cid} | {c['desc'][:80]} | {c['organism']}\n")
                for j in range(0, len(c['sequence']), 80):
                    f.write(c['sequence'][j:j+80] + '\n')

        successful.append({
            **c,
            'safe_name': safe_name,
            'mapped': len(mapping),
            'catalytic_residues': mapping,
            'yaml': yaml_path,
        })

        if (i + 1) % 20 == 0:
            print(f"  Processed {i+1}/{len(selected)}...")

    print(f"\nSuccessful: {len(successful)}, Failed: {len(failed)}")
    if failed:
        print(f"Failed IDs: {[f[0] for f in failed[:10]]}")

    # Save manifest
    manifest_path = os.path.join(SCRIPT_DIR, 'batch_manifest.csv')
    with open(manifest_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'safe_name', 'clean_id', 'db', 'organism', 'kingdom',
            'seq_identity', 'bit_score', 'seq_len', 'mapped',
        ])
        writer.writeheader()
        for s in successful:
            writer.writerow({
                'safe_name': s['safe_name'],
                'clean_id': s['clean_id'],
                'db': s['db'],
                'organism': s['organism'],
                'kingdom': s['kingdom'],
                'seq_identity': f"{s['seq_identity']:.1f}",
                'bit_score': f"{s['bit_score']:.0f}",
                'seq_len': len(s['sequence']),
                'mapped': s['mapped'],
            })
    print(f"Saved manifest: {manifest_path}")

    # Generate batch run script
    # 10 GPUs, process 10 at a time, wait, then next batch
    yaml_names = [s['safe_name'] for s in successful]
    batch_size = 10
    batches = [yaml_names[i:i+batch_size] for i in range(0, len(yaml_names), batch_size)]

    script = f"""#!/bin/bash
# Boltz-2 cofolding: 90 additional homologs in batches of 10
# Each batch uses GPUs 10-19, waits for completion, then starts next batch

CONFIG_DIR={CONFIG_DIR}
OUTPUT_DIR={os.path.join(SCRIPT_DIR, 'output_2mvap')}
LOG=$OUTPUT_DIR/batch_run.log

mkdir -p $OUTPUT_DIR

echo "$(date): Starting batch cofolding — {len(yaml_names)} homologs in {len(batches)} batches" | tee $LOG

"""
    for batch_idx, batch in enumerate(batches):
        script += f"""
# === Batch {batch_idx + 1}/{len(batches)} ({len(batch)} jobs) ===
echo "$(date): Starting batch {batch_idx + 1}/{len(batches)}" | tee -a $LOG
GPU=10
"""
        for name in batch:
            script += f"""
NAME={name}_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \\
            "$CONFIG_DIR/${{NAME}}.yaml" \\
            --out_dir "$OUTPUT_DIR/$NAME" \\
            --use_msa_server \\
            --accelerator gpu \\
            --recycling_steps 3 \\
            --sampling_steps 200 \\
            --diffusion_samples 50 \\
            --output_format pdb \\
            --no_kernels \\
            > "$OUTPUT_DIR/${{NAME}}.log" 2>&1
        echo "$(date): $NAME done (exit code $?)" | tee -a $LOG
    ) &
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))
"""
        script += f"""
echo "$(date): Batch {batch_idx + 1} launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch {batch_idx + 1} complete" | tee -a $LOG
"""

    script += """
echo "$(date): All batches finished!" | tee -a $LOG

# Run analysis automatically
echo "$(date): Running entrance fraction analysis..." | tee -a $LOG
python3 """ + os.path.join(SCRIPT_DIR, 'analyze_all_homologs.py') + """ 2>&1 | tee -a $LOG

echo "$(date): All done!" | tee -a $LOG
"""

    script_path = os.path.join(SCRIPT_DIR, 'run_batch_cofold.sh')
    with open(script_path, 'w') as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    print(f"\nSaved batch script: {script_path}")
    print(f"  {len(batches)} batches of {batch_size}, using GPUs 10-19")
    print(f"  Estimated time: {len(batches) * 30} min ({len(batches) * 0.5:.1f} hours)")

    # Summary
    print(f"\n=== Summary ===")
    print(f"New configs generated: {len(successful)}")
    print(f"Total (with original 10): {len(successful) + 10}")
    print(f"Batch script: {script_path}")


if __name__ == '__main__':
    main()
