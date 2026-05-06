#!/usr/bin/env python3
"""Generate Boltz2 cofolding YAML configs for PMDsc structural homologs.

Selects a diverse set of homologs from Foldseek results, aligns them to PMDsc
to map catalytic residue positions, and generates YAML configs for 2-MVAP
cofolding (same pipeline as the PMDsc variant analysis).
"""

import os
import re
import csv
import subprocess
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
FOLDSEEK_DIR = os.path.join(SCRIPT_DIR, 'foldseek_results')
SEQ_DIR = os.path.join(SCRIPT_DIR, 'sequences')
CONFIG_DIR = os.path.join(SCRIPT_DIR, 'cofold_configs')
BOLTZ_DIR = os.path.join(PROJECT_DIR, 'structures', 'boltz2')

PMDSC_SEQ = (
    "MTVYTASVTAPVNIATLKYWGKRDTKLNLPTNSSISVTLSQDDLRTLTSAATAPEFERDTLWLNGEPHSIDNERTQNCLRDLR"
    "QLRKEMESKDASLPTLSQWKLHIVSENNFPTAAGLASSAAGFAALVSAIAKLYQLPQSTSEISRIARKGSGSACRSLFGGYVAWEM"
    "GKAEDGHDSMAVQIADSSDWPQMKACVLVVSDIKKDVSSTQGMQLTVATSELFKERIEHVVPKRFEVMRKAIVEKDFATFAKETMM"
    "DSNSFHATCLDSFPPIFYMNDTSKRIISWCHTINQFYGETIVAYTFDAGPNAVLYYLAENESKLFAFIYKLFGSVPGWDKKFTTEQL"
    "EAFNHQFESSNFTARELDLELQKDVARVILTQVGSGPQETNESLIDAKTGLPKE"
)

# PMDsc catalytic residues (1-indexed in the full sequence including M1)
CATALYTIC_RESIDUES = {
    18: 'K18',   # substrate binding
    19: 'Y19',   # substrate binding
    20: 'W20',   # substrate binding
    22: 'K22',   # substrate binding
    120: 'S120', # substrate anchoring
    121: 'S121', # ATP binding
    153: 'S153', # substrate anchoring
    155: 'S155', # ATP binding
    158: 'R158', # decarboxylation
    208: 'S208', # ATP binding
    302: 'D302', # catalytic base
}

# L1 pocket constraint residues (same as cofold_2mvap/PMDsc_WT_2mvap.yaml)
L1_POCKET_RESIDUES = [18, 19, 20, 22, 120, 153, 155, 158, 208, 302]
# ATP pocket constraint residues
ATP_POCKET_RESIDUES = [121, 153, 155, 208]
# Mg2+ contact residue
MG_CONTACT_RESIDUE = 302

MVAP_SMILES = 'C[C@](O)(CCOP(=O)(O)O)CC(=O)O'

# Curated diverse homolog set — covering fungi, animals, plants, protozoa, bacteria, archaea
SELECTED_HOMOLOGS = [
    # Close eukaryotic homologs (MDD, same reaction on MVAPP)
    {'id': 'O13963', 'name': 'MDD_Spombe',    'desc': 'S. pombe MDD'},
    {'id': '6XR5',   'name': 'MDD_Cneo',      'desc': 'C. neoformans MDD (PDB)'},
    {'id': '3D4J',   'name': 'MDD_Hsap',      'desc': 'Human MDD (PDB)'},
    {'id': 'Q54YQ9', 'name': 'MDD_Ddis',      'desc': 'D. discoideum MDD'},
    # Bacterial MDD (lower identity, different entrance architecture)
    {'id': '4DPU',   'name': 'MDD_Sepi',      'desc': 'S. epidermidis MDD (PDB)'},
    {'id': '2HK2',   'name': 'MDD_Saur',      'desc': 'S. aureus MDD (PDB)'},
    {'id': '6E2U',   'name': 'MDD_Efae',      'desc': 'E. faecalis MDD (PDB)'},
    {'id': '2GS8',   'name': 'MDD_Spyo',      'desc': 'S. pyogenes MDD (PDB)'},
    # Archaeal/Chloroflexi MVAP decarboxylase (acts on MVAP like PMDsc!)
    {'id': '6N0X',   'name': 'MPD_Ather',     'desc': 'A. thermophila MVAP decarboxylase (PDB)'},
    {'id': '5GMD',   'name': 'MDD_Ssol',      'desc': 'S. solfataricus MDD (PDB, archaeal)'},
]


def pairwise_align(seq1, seq2):
    """Simple Needleman-Wunsch-like alignment using subprocess if available,
    otherwise a basic gap-aware alignment."""
    try:
        import Bio.pairwise2 as pw2
        from Bio.pairwise2 import format_alignment
        alignments = pw2.align.globalms(seq1, seq2, 2, -1, -2, -0.5,
                                         one_alignment_only=True)
        if alignments:
            return alignments[0][0], alignments[0][1]
    except ImportError:
        pass

    try:
        from Bio.Align import PairwiseAligner
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        alignment = aligner.align(seq1, seq2)[0]
        aln_strs = str(alignment).split('\n')
        return aln_strs[0], aln_strs[2]
    except (ImportError, Exception):
        pass

    # Fallback: use EMBOSS needle via subprocess if available
    try:
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f1:
            f1.write(f">query\n{seq1}\n")
            f1_path = f1.name
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f2:
            f2.write(f">target\n{seq2}\n")
            f2_path = f2.name
        out_path = f1_path + '.aln'
        subprocess.run(['needle', '-asequence', f1_path, '-bsequence', f2_path,
                        '-outfile', out_path, '-gapopen', '10', '-gapextend', '0.5',
                        '-auto'], capture_output=True, timeout=30)
        # Parse needle output
        aln1, aln2 = '', ''
        with open(out_path) as f:
            for line in f:
                if line.startswith('query'):
                    aln1 += line.split()[2] if len(line.split()) > 2 else ''
                elif line.startswith('target'):
                    aln2 += line.split()[2] if len(line.split()) > 2 else ''
        os.unlink(f1_path)
        os.unlink(f2_path)
        os.unlink(out_path)
        if aln1 and aln2:
            return aln1, aln2
    except Exception:
        pass

    # Last resort: simple identity-based alignment (no gaps)
    # This is very crude but ensures we always return something
    print("  WARNING: No alignment tool available, using crude positional mapping")
    return seq1, seq2


def map_residues_via_alignment(pmdsc_seq, target_seq):
    """Align target to PMDsc and map catalytic residue positions.

    Returns dict: pmdsc_resnum -> target_resnum (1-indexed).
    """
    aln_q, aln_t = pairwise_align(pmdsc_seq, target_seq)

    # Map alignment positions to sequence positions
    mapping = {}
    q_pos = 0  # 0-indexed position in query (PMDsc)
    t_pos = 0  # 0-indexed position in target

    for aln_idx in range(min(len(aln_q), len(aln_t))):
        q_char = aln_q[aln_idx]
        t_char = aln_t[aln_idx]

        q_is_gap = q_char in ('-', '.')
        t_is_gap = t_char in ('-', '.')

        if not q_is_gap and not t_is_gap:
            q_resnum = q_pos + 1  # 1-indexed
            t_resnum = t_pos + 1
            if q_resnum in CATALYTIC_RESIDUES:
                mapping[q_resnum] = t_resnum
            q_pos += 1
            t_pos += 1
        elif q_is_gap and not t_is_gap:
            t_pos += 1
        elif not q_is_gap and t_is_gap:
            q_pos += 1
        # both gaps: skip

    return mapping


def generate_yaml_config(name, sequence, residue_mapping):
    """Generate Boltz2 cofolding YAML config for a homolog."""
    # Map L1 pocket residues
    l1_contacts = []
    for pmdsc_res in L1_POCKET_RESIDUES:
        if pmdsc_res in residue_mapping:
            t_res = residue_mapping[pmdsc_res]
            l1_contacts.append(['A', t_res])

    # Map ATP pocket residues
    atp_contacts = []
    for pmdsc_res in ATP_POCKET_RESIDUES:
        if pmdsc_res in residue_mapping:
            t_res = residue_mapping[pmdsc_res]
            atp_contacts.append(['A', t_res])

    # Map Mg2+ contact
    mg_contact_res = residue_mapping.get(MG_CONTACT_RESIDUE)

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
        config['constraints'].append({
            'pocket': {
                'binder': 'L',
                'contacts': l1_contacts,
                'max_distance': 6.0,
            }
        })

    if atp_contacts:
        config['constraints'].append({
            'pocket': {
                'binder': 'T',
                'contacts': atp_contacts,
                'max_distance': 6.0,
            }
        })

    if mg_contact_res:
        config['constraints'].append({
            'contact': {
                'token1': ['M', 'MG'],
                'token2': ['A', mg_contact_res],
                'max_distance': 4.0,
            }
        })

    return config


def load_homolog_sequence(homolog_id):
    """Load sequence from foldseek results or local FASTA files."""
    # Check local sequences directory
    fasta_path = os.path.join(SEQ_DIR, f'{homolog_id}.fasta')
    if os.path.exists(fasta_path):
        with open(fasta_path) as f:
            lines = f.readlines()
            seq = ''.join(l.strip() for l in lines if not l.startswith('>'))
            return seq

    # Check in foldseek m8 results (field 19 = target full sequence)
    for fname in ['alis_pdb100.m8', 'alis_afdb-swissprot.m8']:
        fpath = os.path.join(FOLDSEEK_DIR, fname)
        if not os.path.exists(fpath):
            continue
        with open(fpath) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 19:
                    continue
                target_raw = fields[1].split(' ')[0]
                # Match PDB IDs
                pdb_id = target_raw.split('-')[0].upper()
                # Match AF IDs
                af_parts = target_raw.split('-')
                af_id = af_parts[1] if len(af_parts) > 1 and target_raw.startswith('AF-') else ''

                if pdb_id == homolog_id.upper() or af_id == homolog_id:
                    seq = fields[18].replace('-', '')
                    if len(seq) > 50:
                        return seq

    return None


def generate_run_script(homologs_with_configs):
    """Generate a batch run script for Boltz2 cofolding."""
    script = f"""#!/bin/bash
# Boltz-2 cofolding: PMDsc structural homologs with 2x MVAP + ATP + Mg2+
# Generated by prepare_homolog_cofold.py

CONFIG_DIR={CONFIG_DIR}
OUTPUT_DIR={os.path.join(SCRIPT_DIR, 'output_2mvap')}
LOG=$OUTPUT_DIR/run.log

mkdir -p $OUTPUT_DIR

echo "$(date): Starting homolog cofolding" | tee $LOG

GPU=${{1:-0}}  # Use GPU passed as argument, default 0
for YAML in $CONFIG_DIR/*.yaml; do
    NAME=$(basename $YAML .yaml)
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \\
            "$YAML" \\
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
    )
done

echo "$(date): All homolog cofolding jobs finished." | tee -a $LOG
"""
    script_path = os.path.join(SCRIPT_DIR, 'run_homolog_cofold.sh')
    with open(script_path, 'w') as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    print(f"\nSaved run script: {script_path}")
    print(f"Usage: bash {script_path} <GPU_ID>")
    return script_path


def main():
    os.makedirs(CONFIG_DIR, exist_ok=True)

    print("=== Preparing Boltz2 Cofolding Configs for Homologs ===\n")
    print(f"PMDsc sequence length: {len(PMDSC_SEQ)} aa")
    print(f"Catalytic residues: {list(CATALYTIC_RESIDUES.values())}")
    print(f"Selected homologs: {len(SELECTED_HOMOLOGS)}\n")

    successful = []
    for h in SELECTED_HOMOLOGS:
        hid = h['id']
        hname = h['name']
        print(f"--- {hname} ({hid}): {h['desc']} ---")

        # Load sequence
        seq = load_homolog_sequence(hid)
        if not seq:
            print(f"  ERROR: Could not find sequence for {hid}")
            continue
        print(f"  Sequence: {len(seq)} aa")

        # Align and map catalytic residues
        mapping = map_residues_via_alignment(PMDSC_SEQ, seq)
        print(f"  Mapped {len(mapping)}/{len(CATALYTIC_RESIDUES)} catalytic residues:")
        for pmdsc_res, target_res in sorted(mapping.items()):
            pmdsc_aa = PMDSC_SEQ[pmdsc_res - 1] if pmdsc_res <= len(PMDSC_SEQ) else '?'
            target_aa = seq[target_res - 1] if target_res <= len(seq) else '?'
            conserved = '=' if pmdsc_aa == target_aa else '!'
            print(f"    {CATALYTIC_RESIDUES[pmdsc_res]} (PMDsc {pmdsc_aa}{pmdsc_res}) "
                  f"-> {target_aa}{target_res} {conserved}")

        if len(mapping) < 5:
            print(f"  WARNING: Only {len(mapping)} residues mapped, skipping")
            continue

        # Generate YAML config
        config = generate_yaml_config(hname, seq, mapping)
        yaml_path = os.path.join(CONFIG_DIR, f'{hname}_2mvap.yaml')
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        print(f"  Saved config: {yaml_path}")

        # Save FASTA if not already saved
        fasta_path = os.path.join(SEQ_DIR, f'{hid}.fasta')
        if not os.path.exists(fasta_path):
            with open(fasta_path, 'w') as f:
                f.write(f">{hid} | {h['desc']}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

        successful.append({**h, 'seq_len': len(seq), 'mapped': len(mapping),
                          'yaml': yaml_path})

    # Generate run script
    if successful:
        generate_run_script(successful)

    # Summary
    print(f"\n=== Summary ===")
    print(f"Successfully prepared: {len(successful)}/{len(SELECTED_HOMOLOGS)} homologs")
    print(f"Configs saved to: {CONFIG_DIR}")
    print(f"\nHomologs ready for cofolding:")
    print(f"{'Name':<20} {'ID':<12} {'Length':>6} {'Mapped':>7} {'Description'}")
    print("-" * 80)
    for h in successful:
        print(f"{h['name']:<20} {h['id']:<12} {h['seq_len']:>6} "
              f"{h['mapped']:>5}/{len(CATALYTIC_RESIDUES):<2} {h['desc']}")


if __name__ == '__main__':
    main()
