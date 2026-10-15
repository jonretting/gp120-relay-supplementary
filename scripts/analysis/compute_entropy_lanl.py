#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Compute Shannon entropy from LANL HIV-1 Group M gp120 alignment.

Downloads the HIV-1 Group M Env protein alignment from the Los Alamos
HIV Sequence Database, extracts positions 465-472 (HXB2 numbering),
and computes per-residue Shannon entropy.

This validates (or invalidates) the entropy values currently in the paper.
"""

import json, math, os, sys, time
from collections import Counter

try:
    import requests
except ImportError:
    print("ERROR: requests library required. pip install requests")
    sys.exit(1)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
OUTPUT_FILE = os.path.join(SUPP_DIR, "results", "03-conservation", "entropy_validation.json")

# HXB2 reference for relay pocket
HXB2_RELAY = {465: "R", 466: "P", 467: "G", 468: "G", 469: "G", 470: "D", 471: "M", 472: "R"}

# Values currently in the paper (to compare against)
PAPER_VALUES = {
    465: 0.15, 466: 0.08, 467: 0.05, 468: 0.05,
    469: 0.12, 470: 0.22, 471: 0.35, 472: 0.18,
}
PAPER_MEAN = 0.15


def shannon_entropy(column):
    """Compute Shannon entropy for a column of aligned residues."""
    valid = [aa for aa in column if aa not in ("-", ".", "X", "?", "*")]
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


def parse_fasta(text):
    """Parse FASTA format text into list of (header, sequence) tuples."""
    sequences = []
    header = None
    seq_parts = []
    for line in text.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                sequences.append((header, "".join(seq_parts)))
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line.upper())
    if header is not None:
        sequences.append((header, "".join(seq_parts)))
    return sequences


def find_hxb2_index(sequences):
    """Find the HXB2 reference sequence in the alignment."""
    for i, (header, seq) in enumerate(sequences):
        h = header.upper()
        if "HXB2" in h or "K03455" in h or "H.FR.1983.HXB2" in h:
            return i, header, seq
    for i, (header, seq) in enumerate(sequences):
        h = header.upper()
        if "B.FR.1983" in h and "HXB" in h:
            return i, header, seq
    return None, None, None


def alignment_to_hxb2_map(hxb2_seq):
    """Map alignment columns to HXB2 residue numbers.
    
    In gp120, HXB2 numbering starts at residue 1 for the signal peptide.
    Env residue 1 = start of signal peptide.
    gp120 starts after signal peptide cleavage (~residue 31).
    The C5 domain relay pocket is at HXB2 positions 465-472.
    """
    hxb2_pos = 0
    col_to_hxb2 = {}
    for col_idx, aa in enumerate(hxb2_seq):
        if aa not in ("-", "."):
            hxb2_pos += 1
            col_to_hxb2[col_idx] = hxb2_pos
    return col_to_hxb2


def download_lanl_alignment():
    """Download HIV-1 Group M Env alignment from LANL.
    
    Uses the LANL pre-made alignment download interface.
    """
    print("Attempting to download LANL HIV-1 Group M alignment...")
    print()

    # LANL alignment download endpoint
    # The pre-made alignments are at specific URLs
    urls_to_try = [
        # Pre-made filtered web alignment - Env protein, Group M
        {
            "url": "https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi",
            "method": "POST",
            "data": {
                "ORGANISM": "HIV",
                "ALIGN_TYPE": "PRO",
                "SUBORGANISM": "HIV1",
                "PRE_USER": "predefined",
                "REGION": "ENV",
                "GESSION": "LANL",
                "action": "Get Alignment",
                "GROUP": "M",
            },
            "desc": "LANL web alignment tool (Env, Group M, protein)",
        },
        # Direct pre-computed alignment files
        {
            "url": "https://www.hiv.lanl.gov/content/sequence/NEWALIGN/ALIGN_DOWNLOADS/ENV_PRO/ENV_PRO.fasta",
            "method": "GET",
            "desc": "Pre-computed Env protein alignment (all)",
        },
        {
            "url": "https://www.hiv.lanl.gov/content/sequence/NEWALIGN/ALIGN_DOWNLOADS/GP120_PRO/GP120_PRO.fasta",
            "method": "GET",
            "desc": "Pre-computed gp120 protein alignment (all)",
        },
    ]

    for attempt in urls_to_try:
        print(f"  Trying: {attempt['desc']}")
        print(f"  URL: {attempt['url']}")
        try:
            if attempt["method"] == "POST":
                resp = requests.post(
                    attempt["url"],
                    data=attempt.get("data", {}),
                    timeout=60,
                    headers={"User-Agent": "Mozilla/5.0 (research)"},
                )
            else:
                resp = requests.get(
                    attempt["url"],
                    timeout=60,
                    headers={"User-Agent": "Mozilla/5.0 (research)"},
                )

            print(f"  Status: {resp.status_code}")
            print(f"  Content-Type: {resp.headers.get('Content-Type', 'unknown')}")
            print(f"  Size: {len(resp.text)} chars")

            if resp.status_code == 200:
                # Check if it's FASTA
                text = resp.text
                if text.strip().startswith(">"):
                    seqs = parse_fasta(text)
                    print(f"  Parsed {len(seqs)} sequences")
                    if len(seqs) > 10:
                        # Save raw alignment
                        raw_path = os.path.join(SUPP_DIR, "lanl_alignment_raw.fasta")
                        with open(raw_path, "w") as f:
                            f.write(text)
                        print(f"  Saved raw alignment: {raw_path}")
                        return seqs, attempt["desc"]
                else:
                    # Might be HTML - check for FASTA content embedded
                    if ">B.FR.1983" in text or ">HXB2" in text:
                        # Extract FASTA from HTML
                        start = text.find(">")
                        if start >= 0:
                            fasta_text = text[start:]
                            seqs = parse_fasta(fasta_text)
                            if len(seqs) > 10:
                                print(f"  Extracted {len(seqs)} sequences from HTML")
                                return seqs, attempt["desc"]
                    print(f"  Not FASTA format, first 200 chars: {text[:200]}")
        except Exception as e:
            print(f"  Error: {e}")
        print()
        time.sleep(1)

    return None, None


def compute_all_entropy(sequences, relay_positions=range(465, 473)):
    """Compute entropy at relay pocket positions."""
    hxb2_idx, hxb2_header, hxb2_seq = find_hxb2_index(sequences)

    if hxb2_idx is None:
        print("WARNING: HXB2 reference not found in alignment!")
        print("  Checking first 10 headers:")
        for i, (h, s) in enumerate(sequences[:10]):
            print(f"    {i}: {h[:80]}")
        return None

    print(f"\nHXB2 reference found: {hxb2_header[:80]}")
    print(f"  Alignment length: {len(hxb2_seq)} columns")
    print(f"  Total sequences: {len(sequences)}")

    col_to_hxb2 = alignment_to_hxb2_map(hxb2_seq)
    hxb2_to_col = {v: k for k, v in col_to_hxb2.items()}

    print(f"  HXB2 residues mapped: {len(col_to_hxb2)}")
    print(f"  HXB2 range: {min(col_to_hxb2.values())} to {max(col_to_hxb2.values())}")

    # Check if relay positions are in range
    for pos in relay_positions:
        if pos not in hxb2_to_col:
            print(f"  WARNING: HXB2 position {pos} not found in alignment map!")

    results = {}
    print(f"\n{'Pos':>4}  {'HXB2':>4}  {'Entropy':>8}  {'Paper':>8}  {'Diff':>8}  {'N_seq':>6}  {'Top AAs'}")
    print("-" * 75)

    for pos in relay_positions:
        col = hxb2_to_col.get(pos)
        if col is None:
            print(f"  {pos:>3}  MISSING")
            continue

        column = [seq[col] for _, seq in sequences if col < len(seq)]
        valid = [aa for aa in column if aa not in ("-", ".", "X", "?", "*")]
        entropy = shannon_entropy(column)
        paper_val = PAPER_VALUES.get(pos, None)
        diff = entropy - paper_val if paper_val is not None else None

        # Top amino acids
        counts = Counter(valid)
        total = len(valid)
        top = counts.most_common(5)
        top_str = ", ".join(f"{aa}:{c}/{total}({100*c/total:.0f}%)" for aa, c in top)

        hxb2_aa = HXB2_RELAY.get(pos, "?")
        diff_str = f"{diff:+.3f}" if diff is not None else "N/A"
        paper_str = f"{paper_val:.2f}" if paper_val is not None else "N/A"

        print(f"  {pos:>3}  {hxb2_aa:>4}  {entropy:>8.4f}  {paper_str:>8}  {diff_str:>8}  {total:>6}  {top_str}")

        results[pos] = {
            "hxb2_aa": hxb2_aa,
            "entropy": round(entropy, 4),
            "paper_value": paper_val,
            "difference": round(diff, 4) if diff is not None else None,
            "n_sequences": total,
            "aa_distribution": {aa: c for aa, c in counts.most_common()},
        }

    # Compute mean
    entropies = [r["entropy"] for r in results.values()]
    computed_mean = sum(entropies) / len(entropies) if entropies else 0
    print(f"\n  Computed mean entropy: {computed_mean:.4f}")
    print(f"  Paper mean entropy:   {PAPER_MEAN}")
    print(f"  Difference:           {computed_mean - PAPER_MEAN:+.4f}")

    # Also compute entropy for broader gp120 domains for comparison
    print("\n\nDOMAIN COMPARISON:")
    print("-" * 60)
    domains = {
        "V1 (131-157)": range(131, 158),
        "V2 (158-196)": range(158, 197),
        "C2 (196-296)": range(196, 297),
        "V3 (296-331)": range(296, 332),
        "C3 (331-385)": range(331, 386),
        "V4 (385-418)": range(385, 419),
        "C4 (418-461)": range(418, 462),
        "C5 (461-485)": range(461, 486),
        "Relay (465-472)": range(465, 473),
    }

    paper_domain_values = {
        "V1 (131-157)": 1.82, "V2 (158-196)": 1.45, "C2 (196-296)": 0.42,
        "V3 (296-331)": 0.95, "C3 (331-385)": 0.38, "V4 (385-418)": 1.12,
        "C4 (418-461)": 0.45, "C5 (461-485)": 0.25, "Relay (465-472)": 0.15,
    }

    domain_results = {}
    print(f"{'Domain':<20}  {'Computed':>10}  {'Paper':>10}  {'Diff':>10}  {'Positions':>10}")
    print("-" * 65)
    for domain_name, pos_range in domains.items():
        domain_entropies = []
        for pos in pos_range:
            col = hxb2_to_col.get(pos)
            if col is not None:
                column = [seq[col] for _, seq in sequences if col < len(seq)]
                domain_entropies.append(shannon_entropy(column))

        if domain_entropies:
            domain_mean = sum(domain_entropies) / len(domain_entropies)
            paper_val = paper_domain_values.get(domain_name)
            diff = domain_mean - paper_val if paper_val else None
            diff_str = f"{diff:+.3f}" if diff is not None else "N/A"
            paper_str = f"{paper_val:.2f}" if paper_val else "N/A"
            print(f"  {domain_name:<18}  {domain_mean:>10.4f}  {paper_str:>10}  {diff_str:>10}  {len(domain_entropies):>10}")
            domain_results[domain_name] = {
                "computed_mean": round(domain_mean, 4),
                "paper_value": paper_val,
                "n_positions": len(domain_entropies),
            }

    return results, domain_results, computed_mean


def main():
    print("=" * 75)
    print("SHANNON ENTROPY VALIDATION")
    print("Comparing paper values against LANL HIV-1 Group M alignment")
    print("=" * 75)
    print()

    sequences, source_desc = download_lanl_alignment()

    if sequences is None:
        print("\nFailed to download alignment from LANL.")
        print("Checking for local alignment file...")

        local_path = os.path.join(SUPP_DIR, "lanl_alignment_raw.fasta")
        if os.path.exists(local_path):
            with open(local_path) as f:
                sequences = parse_fasta(f.read())
            source_desc = f"Local file: {local_path}"
            print(f"  Loaded {len(sequences)} sequences from local file")
        else:
            print("  No local alignment found.")
            print("\nTo run manually:")
            print("  1. Go to https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html")
            print("  2. Select: Organism=HIV-1, Region=ENV or GP120, Alignment type=Protein")
            print("  3. Select: Subtype=M (Group M)")
            print("  4. Download the alignment as FASTA")
            print(f"  5. Save as: {local_path}")
            print("  6. Re-run this script")
            return

    print(f"\nSource: {source_desc}")
    print(f"Sequences: {len(sequences)}")

    result = compute_all_entropy(sequences)
    if result is None:
        print("\nFailed to compute entropy. Check alignment format.")
        return

    per_residue, domain_results, computed_mean = result

    # Save results
    output = {
        "analysis": "Shannon entropy validation against LANL alignment",
        "source": source_desc,
        "n_sequences": len(sequences),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "relay_pocket_per_residue": per_residue,
        "domain_comparison": domain_results,
        "computed_relay_mean": round(computed_mean, 4),
        "paper_relay_mean": PAPER_MEAN,
        "verdict": (
            "MATCH" if abs(computed_mean - PAPER_MEAN) < 0.05
            else "CLOSE" if abs(computed_mean - PAPER_MEAN) < 0.15
            else "MISMATCH"
        ),
    }

    with open(OUTPUT_FILE, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved results: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
