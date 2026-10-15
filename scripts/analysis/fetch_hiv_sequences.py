#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Fetch HIV-1 Group M gp120/Env sequences from NCBI and compute Shannon entropy.

Since LANL downloads require browser auth, we use NCBI Entrez + Clustal alignment.
Alternative: uses pre-aligned sequences if available.
"""

import json, math, os, sys, time
from collections import Counter

import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))

HXB2_RELAY = {465: "R", 466: "P", 467: "G", 468: "G", 469: "G", 470: "D", 471: "M", 472: "R"}
PAPER_VALUES = {
    465: 0.15, 466: 0.08, 467: 0.05, 468: 0.05,
    469: 0.12, 470: 0.22, 471: 0.35, 472: 0.18,
}

# HXB2 Env polyprotein (gp160) sequence - the C5/relay region
# HXB2 GenBank K03455, Env positions 465-472 = RPGGGDMR
# Full gp120 C5 region for context (HXB2 ~461-489)
HXB2_C5_CONTEXT = "NNTSVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSENFTNNAKTIIVQLTTPVQINCTRPGNNTRRSIHIGPGRAFYATGDIIGDIRQAHCNLSRAKWNNTLKQIVIKLREQFG"
# The above is HXB2 Env positions ~420-530 area
# Relay = RPGGGDMR at 465-472


def shannon_entropy(column):
    valid = [aa for aa in column if aa not in ("-", ".", "X", "?", "*", "B", "Z", "J")]
    if not valid:
        return 0.0
    total = len(valid)
    counts = Counter(valid)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def fetch_hiv1_env_from_ncbi():
    """Fetch HIV-1 subtype reference Env sequences from NCBI."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Strategy: fetch the HIV-1 reference genomes for each major subtype
    # These are curated reference sequences with known clade assignments
    
    # Known reference accessions for HIV-1 Group M subtypes (Env/gp160)
    # From LANL reference set / NCBI RefSeq
    refs = {
        "B_HXB2": "AAB50262.1",      # Clade B reference (HXB2)
        "B_NL43": "AAK08484.1",       # Clade B (NL4-3)
        "C_ref": "AY772699",          # Clade C  
        "A1_ref": "AF004885",         # Clade A1
        "D_ref": "K03454",            # Clade D
        "AE_ref": "U54771",           # CRF01_AE
        "AG_ref": "L39106",           # CRF02_AG
        "F_ref": "AF005494",          # Clade F
        "G_ref": "U88826",            # Clade G
    }

    # Better approach: search for a large set of HIV-1 Env sequences
    query = (
        '"Human immunodeficiency virus 1"[Organism] AND '
        'envelope[Title] AND '
        '800:900[Sequence Length]'
    )
    
    print(f"Searching NCBI: {query}")
    r = requests.get(f"{base}/esearch.fcgi", params={
        "db": "protein",
        "term": query,
        "retmax": 1000,
        "retmode": "json",
        "usehistory": "y",
    }, timeout=30)
    
    data = r.json()["esearchresult"]
    count = data["count"]
    webenv = data.get("webenv", "")
    query_key = data.get("querykey", "")
    ids = data["idlist"]
    print(f"  Found {count} sequences, retrieving {len(ids)}")
    
    if not ids:
        print("  No sequences found!")
        return None
    
    # Fetch sequences in FASTA using webenv/query_key (avoids URI length limits)
    print("  Downloading sequences (via history)...")
    r2 = requests.get(f"{base}/efetch.fcgi", params={
        "db": "protein",
        "query_key": query_key,
        "WebEnv": webenv,
        "retstart": 0,
        "retmax": 500,
        "rettype": "fasta",
        "retmode": "text",
    }, timeout=120)
    
    if r2.status_code != 200:
        print(f"  Fetch failed: {r2.status_code}")
        return None
    
    text = r2.text
    sequences = []
    header = None
    seq_parts = []
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                sequences.append((header, "".join(seq_parts)))
            header = line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line.upper())
    if header is not None:
        sequences.append((header, "".join(seq_parts)))
    
    print(f"  Got {len(sequences)} sequences")
    return sequences


def find_relay_in_sequence(seq, motif="RPGGGDMR"):
    """Find the relay motif or close variants in a sequence."""
    # Exact match
    idx = seq.find(motif)
    if idx >= 0:
        return idx, motif, True
    
    # 7-of-8 match
    for i in range(len(seq) - 7):
        window = seq[i:i+8]
        matches = sum(1 for a, b in zip(window, motif) if a == b)
        if matches >= 7:
            return i, window, False
    
    # 6-of-8 match (broader)
    for i in range(len(seq) - 7):
        window = seq[i:i+8]
        matches = sum(1 for a, b in zip(window, motif) if a == b)
        if matches >= 6:
            return i, window, False
    
    return -1, "", False


def compute_entropy_from_unaligned(sequences):
    """Compute per-position entropy by aligning relay regions."""
    print("\nAligning relay regions across sequences...")
    
    relay_sequences = []
    for header, seq in sequences:
        idx, motif, exact = find_relay_in_sequence(seq)
        if idx >= 0:
            # Extract the 8-residue window
            window = seq[idx:idx+8]
            if len(window) == 8:
                relay_sequences.append((header, window, exact))
    
    print(f"  Found relay motif in {len(relay_sequences)} / {len(sequences)} sequences")
    exact_count = sum(1 for _, _, e in relay_sequences if e)
    print(f"  Exact RPGGGDMR: {exact_count}")
    print(f"  Variants (6-7/8 match): {len(relay_sequences) - exact_count}")
    
    if len(relay_sequences) < 20:
        print("  Too few sequences for reliable entropy calculation")
        return None
    
    # Compute per-position entropy
    print(f"\n{'Pos':>4}  {'HXB2':>4}  {'Entropy':>8}  {'Paper':>8}  {'Diff':>8}  {'N':>6}  Top AAs")
    print("-" * 80)
    
    results = {}
    for pos_offset in range(8):
        hxb2_pos = 465 + pos_offset
        column = [seq[pos_offset] for _, seq, _ in relay_sequences]
        entropy = shannon_entropy(column)
        paper_val = PAPER_VALUES.get(hxb2_pos, None)
        diff = entropy - paper_val if paper_val is not None else None
        
        counts = Counter(column)
        total = len(column)
        top = counts.most_common(5)
        top_str = ", ".join(f"{aa}:{c}({100*c/total:.0f}%)" for aa, c in top)
        
        hxb2_aa = HXB2_RELAY[hxb2_pos]
        diff_str = f"{diff:+.4f}" if diff is not None else "N/A"
        paper_str = f"{paper_val:.2f}" if paper_val is not None else "N/A"
        
        print(f"  {hxb2_pos:>3}  {hxb2_aa:>4}  {entropy:>8.4f}  {paper_str:>8}  {diff_str:>8}  {total:>6}  {top_str}")
        
        results[hxb2_pos] = {
            "hxb2_aa": hxb2_aa,
            "entropy": round(entropy, 4),
            "paper_value": paper_val,
            "difference": round(diff, 4) if diff is not None else None,
            "n_sequences": total,
            "aa_distribution": {aa: c for aa, c in counts.most_common()},
            "conservation_pct": round(100 * counts.most_common(1)[0][1] / total, 1) if counts else 0,
        }
    
    entropies = [r["entropy"] for r in results.values()]
    computed_mean = sum(entropies) / len(entropies)
    
    print(f"\n  Computed mean: {computed_mean:.4f}")
    print(f"  Paper mean:    {0.15:.4f}")
    print(f"  Difference:    {computed_mean - 0.15:+.4f}")
    
    # Show variant sequences
    variants = Counter()
    for _, seq, exact in relay_sequences:
        if not exact:
            variants[seq] += 1
    
    if variants:
        print(f"\n  Variant sequences ({len(variants)} unique):")
        for v, c in variants.most_common(10):
            diffs = "".join("^" if a != b else " " for a, b in zip(v, "RPGGGDMR"))
            print(f"    {v}  (n={c})  {diffs}")
    
    return results, computed_mean


def main():
    print("=" * 80)
    print("SHANNON ENTROPY VALIDATION -- from NCBI HIV-1 Env sequences")
    print("=" * 80)
    
    sequences = fetch_hiv1_env_from_ncbi()
    if sequences is None:
        return
    
    result = compute_entropy_from_unaligned(sequences)
    if result is None:
        return
    
    per_residue, computed_mean = result
    
    output = {
        "analysis": "Shannon entropy validation from NCBI HIV-1 Env sequences",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_sequences_searched": len(sequences),
        "relay_pocket": per_residue,
        "computed_mean": round(computed_mean, 4),
        "paper_mean": 0.15,
        "verdict": (
            "MATCH" if abs(computed_mean - 0.15) < 0.05
            else "CLOSE" if abs(computed_mean - 0.15) < 0.15
            else "MISMATCH"
        ),
    }
    
    out_path = os.path.join(SUPP_DIR, "results", "03-conservation", "entropy_validation.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
