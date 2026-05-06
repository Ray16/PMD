#!/usr/bin/env python3
"""Query the Foldseek web API with PMDsc (1FI4) to find structural homologs.

Searches pdb100 and afdb-swissprot databases, parses results, and outputs
a curated list of homologs for Boltz2 cofolding analysis.
"""

import os
import sys
import time
import json
import tarfile
import io
import requests
import csv

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_PATH = os.path.join(SCRIPT_DIR, '..', 'PDB', '1FI4.pdb')
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'foldseek_results')
SEQUENCES_DIR = os.path.join(SCRIPT_DIR, 'sequences')

FOLDSEEK_API = 'https://search.foldseek.com/api'
DATABASES = ['pdb100', 'afdb-swissprot']
MODE = '3diaa'

POLL_INTERVAL = 10
MAX_POLLS = 60


def submit_search(pdb_path, databases, mode='3diaa'):
    """Submit a structure search to the Foldseek web API."""
    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    data = {
        'mode': mode,
    }
    files = {
        'q': ('query.pdb', pdb_content, 'application/octet-stream'),
    }
    for db in databases:
        data.setdefault('database[]', [])
    form_data = [('mode', mode)]
    for db in databases:
        form_data.append(('database[]', db))

    resp = requests.post(
        f'{FOLDSEEK_API}/ticket',
        data=form_data,
        files=files,
    )
    resp.raise_for_status()
    result = resp.json()
    ticket_id = result['id']
    print(f"Submitted job: {ticket_id}")
    return ticket_id


def poll_status(ticket_id):
    """Poll until job completes."""
    for i in range(MAX_POLLS):
        resp = requests.get(f'{FOLDSEEK_API}/ticket/{ticket_id}')
        resp.raise_for_status()
        status_data = resp.json()
        status = status_data.get('status', 'UNKNOWN')

        if status == 'COMPLETE':
            print(f"Job completed after {i * POLL_INTERVAL}s")
            return True
        elif status == 'ERROR':
            print(f"Job failed: {status_data}")
            return False
        else:
            print(f"  Poll {i+1}/{MAX_POLLS}: status={status}")
            time.sleep(POLL_INTERVAL)

    print("Timed out waiting for results")
    return False


def download_results(ticket_id):
    """Download and extract result files."""
    resp = requests.get(f'{FOLDSEEK_API}/result/download/{ticket_id}',
                        stream=True)
    resp.raise_for_status()

    content = b''
    for chunk in resp.iter_content(chunk_size=8192):
        content += chunk

    results = {}
    try:
        tar = tarfile.open(fileobj=io.BytesIO(content), mode='r:gz')
        for member in tar.getmembers():
            f = tar.extractfile(member)
            if f:
                file_content = f.read().decode('utf-8', errors='replace')
                results[member.name] = file_content
                outpath = os.path.join(RESULTS_DIR, os.path.basename(member.name))
                with open(outpath, 'w') as out:
                    out.write(file_content)
                print(f"  Extracted: {member.name} ({len(file_content)} bytes)")
    except tarfile.ReadError:
        outpath = os.path.join(RESULTS_DIR, 'raw_results.txt')
        with open(outpath, 'wb') as out:
            out.write(content)
        print(f"  Saved raw results ({len(content)} bytes)")
        results['raw'] = content.decode('utf-8', errors='replace')

    return results


def parse_m8_results(results):
    """Parse Foldseek extended m8 format (21 fields).

    Fields: query, target+desc, identity, aln_len, mismatches, gap_opens,
            q_start, q_end, t_start, t_end, prob, evalue, bit_score,
            q_len, t_len, q_aln_seq, t_aln_seq, coords, t_full_seq, taxid, organism
    """
    all_hits = []

    for filename, content in results.items():
        if not content.strip() or '_report' in filename:
            continue

        lines = content.strip().split('\n')
        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 15:
                continue

            target_full = fields[1]
            # Parse target: "pdbid-assembly.cif.gz_chain description" or "AF-XXXX-F1-model_v6 description"
            target_parts = target_full.split(' ', 1)
            target_id = target_parts[0]
            target_desc = target_parts[1] if len(target_parts) > 1 else ''

            # Extract clean PDB ID or UniProt ID
            if target_id.startswith('AF-'):
                # AlphaFold: AF-Q5HII1-F1-model_v6 → Q5HII1
                parts = target_id.split('-')
                clean_id = parts[1] if len(parts) > 1 else target_id
                db_source = 'alphafold'
            else:
                # PDB: "1fi4-assembly1.cif.gz_A-2" → "1fi4", chain "A"
                pdb_part = target_id.split('-')[0]
                chain_part = target_id.split('_')[-1].split('-')[0] if '_' in target_id else 'A'
                clean_id = pdb_part.upper()
                db_source = 'pdb'

            hit = {
                'query': fields[0],
                'target_raw': target_id,
                'target_desc': target_desc,
                'clean_id': clean_id,
                'db_source': db_source,
                'seq_identity': float(fields[2]),
                'alignment_length': int(fields[3]),
                'mismatches': int(fields[4]),
                'gap_opens': int(fields[5]),
                'q_start': int(fields[6]),
                'q_end': int(fields[7]),
                't_start': int(fields[8]),
                't_end': int(fields[9]),
                'prob': float(fields[10]),
                'evalue': float(fields[11]),
                'bit_score': float(fields[12]),
                'q_len': int(fields[13]),
                't_len': int(fields[14]),
                'source_file': filename,
            }

            if len(fields) > 15:
                hit['q_aln_seq'] = fields[15]
            if len(fields) > 16:
                hit['t_aln_seq'] = fields[16]
            if len(fields) > 18:
                hit['t_full_seq'] = fields[18]
            if len(fields) > 19:
                hit['taxid'] = fields[19]
            if len(fields) > 20:
                hit['organism'] = fields[20]

            # Extract chain for PDB hits
            if db_source == 'pdb' and '_' in target_id:
                chain = target_id.split('_')[-1].split('-')[0]
                hit['chain'] = chain

            all_hits.append(hit)

    return all_hits


def curate_homologs(hits, min_prob=0.5, min_identity=15, max_identity=95,
                    min_aln_length=150):
    """Filter and deduplicate hits to select diverse homologs."""
    filtered = []
    for h in hits:
        identity = h['seq_identity']
        aln_len = h['alignment_length']
        prob = h.get('prob', 0)

        if identity < min_identity or identity > max_identity:
            continue
        if aln_len < min_aln_length:
            continue
        if prob < min_prob:
            continue

        filtered.append(h)

    # Sort by prob descending, then by bit_score
    filtered.sort(key=lambda x: (x.get('prob', 0), x['bit_score']), reverse=True)

    # Deduplicate by clean_id (take best hit per protein)
    seen = set()
    unique = []
    for h in filtered:
        cid = h['clean_id']
        if cid not in seen:
            seen.add(cid)
            unique.append(h)

    return unique


def extract_sequences(hits, max_n=15):
    """Extract sequences from Foldseek results (already included in output).

    Falls back to fetching from PDB/UniProt if the sequence field is missing.
    """
    sequences = {}

    for h in hits[:max_n]:
        clean_id = h['clean_id']
        organism = h.get('organism', 'Unknown')

        # Foldseek includes the target sequence in field 19 (t_full_seq)
        if h.get('t_full_seq') and len(h['t_full_seq']) > 50:
            seq = h['t_full_seq'].replace('-', '')
            header = f">{clean_id} | {h.get('target_desc', '')} | {organism}"
            sequences[clean_id] = {
                'header': header,
                'sequence': seq,
                'source': 'foldseek',
                'organism': organism,
                'hit': h,
            }
            print(f"  Extracted {clean_id}: {len(seq)} aa ({organism})")
            continue

        # Fallback: fetch from PDB or UniProt
        if h['db_source'] == 'pdb':
            try:
                url = f'https://www.rcsb.org/fasta/entry/{clean_id}/download'
                resp = requests.get(url, timeout=15)
                if resp.status_code == 200 and resp.text.strip().startswith('>'):
                    lines = resp.text.strip().split('\n')
                    header = lines[0]
                    seq = ''.join(l.strip() for l in lines[1:] if not l.startswith('>'))
                    sequences[clean_id] = {
                        'header': header,
                        'sequence': seq,
                        'source': 'PDB',
                        'organism': organism,
                        'hit': h,
                    }
                    print(f"  Fetched {clean_id}: {len(seq)} aa from PDB")
                    continue
            except Exception:
                pass
        else:
            try:
                url = f'https://rest.uniprot.org/uniprotkb/{clean_id}.fasta'
                resp = requests.get(url, timeout=15)
                if resp.status_code == 200 and resp.text.strip().startswith('>'):
                    lines = resp.text.strip().split('\n')
                    header = lines[0]
                    seq = ''.join(l.strip() for l in lines[1:])
                    sequences[clean_id] = {
                        'header': header,
                        'sequence': seq,
                        'source': 'UniProt',
                        'organism': organism,
                        'hit': h,
                    }
                    print(f"  Fetched {clean_id}: {len(seq)} aa from UniProt")
                    continue
            except Exception:
                pass

        print(f"  Could not get sequence for {clean_id}")

    return sequences


def save_results(hits, sequences):
    """Save curated results as CSV and FASTA files."""
    csv_path = os.path.join(RESULTS_DIR, 'curated_homologs.csv')
    fieldnames = ['rank', 'clean_id', 'db_source', 'description', 'organism',
                  'seq_identity', 'alignment_length', 'prob', 'evalue',
                  'bit_score', 'q_len', 't_len', 'has_sequence', 'sequence_length']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, h in enumerate(hits):
            seq_info = sequences.get(h['clean_id'], {})
            writer.writerow({
                'rank': i + 1,
                'clean_id': h['clean_id'],
                'db_source': h['db_source'],
                'description': h.get('target_desc', ''),
                'organism': h.get('organism', ''),
                'seq_identity': f"{h['seq_identity']:.1f}",
                'alignment_length': h['alignment_length'],
                'prob': f"{h['prob']:.4f}",
                'evalue': f"{h['evalue']:.2e}",
                'bit_score': f"{h['bit_score']:.1f}",
                'q_len': h.get('q_len', ''),
                't_len': h.get('t_len', ''),
                'has_sequence': h['clean_id'] in sequences,
                'sequence_length': len(seq_info.get('sequence', '')),
            })
    print(f"\nSaved curated homologs to {csv_path}")

    # Save FASTA files
    for target_id, info in sequences.items():
        safe_name = target_id.replace('/', '_').replace(' ', '_')
        fasta_path = os.path.join(SEQUENCES_DIR, f'{safe_name}.fasta')
        with open(fasta_path, 'w') as f:
            f.write(f"{info['header']}\n")
            seq = info['sequence']
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
        print(f"  Saved {fasta_path}")


def load_existing_results():
    """Load previously downloaded Foldseek results from disk."""
    results = {}
    for fname in os.listdir(RESULTS_DIR):
        if fname.endswith('.m8'):
            fpath = os.path.join(RESULTS_DIR, fname)
            with open(fpath) as f:
                results[fname] = f.read()
    return results


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(SEQUENCES_DIR, exist_ok=True)

    # Check if we already have results on disk
    existing_m8 = [f for f in os.listdir(RESULTS_DIR) if f.endswith('.m8')]
    if existing_m8:
        print(f"Found existing results: {existing_m8}")
        print("Parsing existing results (skip API query)...")
        results = load_existing_results()
    else:
        pdb_path = os.path.abspath(PDB_PATH)
        if not os.path.exists(pdb_path):
            print(f"ERROR: PDB file not found: {pdb_path}")
            sys.exit(1)

        print(f"Querying Foldseek with: {pdb_path}")
        print(f"Databases: {DATABASES}")
        print(f"Mode: {MODE}")

        # Step 1: Submit search
        ticket_id = submit_search(pdb_path, DATABASES, MODE)

        # Step 2: Poll for completion
        success = poll_status(ticket_id)
        if not success:
            print("Search failed or timed out")
            sys.exit(1)

        # Step 3: Download results
        print("\nDownloading results...")
        results = download_results(ticket_id)

    # Step 4: Parse results
    print("\nParsing results...")
    hits = parse_m8_results(results)
    print(f"Total hits: {len(hits)}")

    if not hits:
        print("No hits found. Check raw results in foldseek_results/")
        # Save raw for debugging
        raw_path = os.path.join(RESULTS_DIR, 'all_results_raw.json')
        with open(raw_path, 'w') as f:
            json.dump({k: v[:500] for k, v in results.items()}, f, indent=2)
        sys.exit(1)

    # Step 5: Curate
    print("\nCurating homologs...")
    curated = curate_homologs(hits)
    print(f"Curated homologs (after filtering): {len(curated)}")

    for i, h in enumerate(curated[:20]):
        print(f"  {i+1:3d}. {h['clean_id']:15s}  identity={h['seq_identity']:.1f}%  "
              f"aln_len={h['alignment_length']}  prob={h['prob']:.3f}  evalue={h['evalue']:.1e}  "
              f"bit_score={h['bit_score']:.0f}  {h.get('organism', '')}")

    # Step 6: Extract/fetch sequences for top hits
    print(f"\nExtracting sequences for top {min(15, len(curated))} hits...")
    sequences = extract_sequences(curated, max_n=15)
    print(f"Successfully obtained {len(sequences)} sequences")

    # Step 7: Save everything
    save_results(curated, sequences)

    print(f"\n=== Summary ===")
    print(f"Total hits: {len(hits)}")
    print(f"Curated homologs: {len(curated)}")
    print(f"Sequences fetched: {len(sequences)}")
    print(f"Results saved to: {RESULTS_DIR}")
    print(f"Sequences saved to: {SEQUENCES_DIR}")


if __name__ == '__main__':
    main()
