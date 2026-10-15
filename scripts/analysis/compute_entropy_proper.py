#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Proper Shannon entropy computation for HIV-1 gp120 relay pocket.

Strategy:
  1. Fetch HIV-1 gp120 sequences from NCBI (full-length gp120, ~480-520 residues)
  2. Pairwise-align each against HXB2 gp120 to establish numbering
  3. Extract positions 465-472 (HXB2 numbering) from each sequence
  4. Compute Shannon entropy per position
  5. Also compute domain-level entropy for comparison

Uses Biopython for sequence retrieval and pairwise alignment.
"""

import json, math, os, sys, time
from collections import Counter

import requests
from Bio import Align

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
OUTPUT_FILE = os.path.join(SUPP_DIR, "results", "03-conservation", "entropy_validation_proper.json")

# HXB2 gp120 sequence (GenBank K03455, Env processed to gp120 only)
# gp120 = Env residues ~31-511 (after signal peptide cleavage)
# HXB2 Env numbering used in the paper: 465-472 = RPGGGDMR in the C5 domain
# This is the full gp120 amino acid sequence from HXB2
HXB2_GP120 = (
    "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKA"  # 1-60
    "YDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCV"  # 61-120
    "KLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDI"  # 121-180
    "IPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNV"  # 181-240
    "STVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNT"  # 241-300
    "RKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQS"  # 301-360
    "SGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINM"  # 361-420
    "WQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVV"  # 421-480
    "KIEPLGVAPTKAKRRVVQREKR"  # 481-502
)

# Clean it (remove any whitespace)
HXB2_GP120 = HXB2_GP120.replace("\n", "").replace(" ", "")

# Find the actual position of RPGGGDMR in our HXB2 sequence (1-indexed)
_rpgggdmr_idx = HXB2_GP120.find("RPGGGDMR")
assert _rpgggdmr_idx >= 0, "RPGGGDMR not found in HXB2 gp120 sequence!"
RELAY_START = _rpgggdmr_idx + 1  # 1-indexed
RELAY_POSITIONS_HXB2 = list(range(RELAY_START, RELAY_START + 8))
HXB2_RELAY = {RELAY_START + i: aa for i, aa in enumerate("RPGGGDMR")}

# Paper uses positions 465-472 (different numbering convention);
# we map to our HXB2 sequence positions and label with paper convention for reporting
PAPER_POS_OFFSET = 465 - RELAY_START  # offset to convert our numbering to paper numbering

PAPER_VALUES = {
    465: 0.15, 466: 0.08, 467: 0.05, 468: 0.05,
    469: 0.12, 470: 0.22, 471: 0.35, 472: 0.18,
}

PAPER_DOMAIN_VALUES = {
    "V1 (131-157)": 1.82, "V2 (158-196)": 1.45, "C2 (196-296)": 0.42,
    "V3 (296-331)": 0.95, "C3 (331-385)": 0.38, "V4 (385-418)": 1.12,
    "C4 (418-461)": 0.45, "C5 (461-485)": 0.25, "Relay (465-472)": 0.15,
}


def shannon_entropy(column):
    valid = [aa for aa in column if aa not in ("-", ".", "X", "?", "*", "B", "Z", "J", "")]
    if len(valid) < 2:
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


def fetch_gp120_sequences():
    """Fetch HIV-1 gp120 sequences from NCBI."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # gp120 is ~480-520 aa; search for that length range to avoid gp160
    query = (
        '"Human immunodeficiency virus 1"[Organism] AND '
        '(gp120[Title] OR "envelope glycoprotein gp120"[Title]) AND '
        '450:550[Sequence Length]'
    )

    print(f"Searching NCBI: {query[:80]}...")
    r = requests.get(f"{base}/esearch.fcgi", params={
        "db": "protein",
        "term": query,
        "retmax": 2000,
        "retmode": "json",
        "usehistory": "y",
    }, timeout=30)

    data = r.json()["esearchresult"]
    count = data["count"]
    webenv = data.get("webenv", "")
    query_key = data.get("querykey", "")
    print(f"  Found {count} sequences")

    # Fetch in batches
    all_text = ""
    batch_size = 200
    fetched = 0
    max_fetch = min(int(count), 1000)

    while fetched < max_fetch:
        n = min(batch_size, max_fetch - fetched)
        print(f"  Fetching {fetched+1}-{fetched+n}...")
        r2 = requests.get(f"{base}/efetch.fcgi", params={
            "db": "protein",
            "query_key": query_key,
            "WebEnv": webenv,
            "retstart": fetched,
            "retmax": n,
            "rettype": "fasta",
            "retmode": "text",
        }, timeout=120)

        if r2.status_code != 200:
            print(f"  Fetch failed at batch {fetched}: {r2.status_code}")
            break

        all_text += r2.text
        fetched += n
        time.sleep(0.5)

    sequences = parse_fasta(all_text)
    print(f"  Total sequences parsed: {len(sequences)}")
    return sequences


def align_to_hxb2(query_seq, hxb2_seq):
    """Pairwise align query to HXB2 and return position mapping."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(hxb2_seq, query_seq)
    if not alignments:
        return None

    aln = alignments[0]

    # Build mapping: HXB2 position -> query amino acid
    hxb2_pos = 0
    query_pos = 0
    mapping = {}

    # Parse the alignment
    aligned = aln.aligned
    # aligned is a tuple of two arrays of (start, end) pairs
    hxb2_segments = aligned[0]
    query_segments = aligned[1]

    # Reconstruct full alignment strings
    hxb2_aligned = ""
    query_aligned = ""

    h_idx = 0
    q_idx = 0

    for (hs, he), (qs, qe) in zip(hxb2_segments, query_segments):
        # Gaps in query before this segment
        while h_idx < hs:
            hxb2_aligned += hxb2_seq[h_idx]
            query_aligned += "-"
            h_idx += 1
        while q_idx < qs:
            hxb2_aligned += "-"
            query_aligned += query_seq[q_idx]
            q_idx += 1

        # Aligned segment
        seg_len = he - hs
        hxb2_aligned += hxb2_seq[hs:he]
        query_aligned += query_seq[qs:qe]
        h_idx = he
        q_idx = qe

    # Trailing
    while h_idx < len(hxb2_seq):
        hxb2_aligned += hxb2_seq[h_idx]
        query_aligned += "-"
        h_idx += 1
    while q_idx < len(query_seq):
        hxb2_aligned += "-"
        query_aligned += query_seq[q_idx]
        q_idx += 1

    # Now map HXB2 positions to query residues
    hxb2_resnum = 0
    result = {}
    for h_aa, q_aa in zip(hxb2_aligned, query_aligned):
        if h_aa != "-":
            hxb2_resnum += 1
            if q_aa != "-":
                result[hxb2_resnum] = q_aa

    return result


def main():
    print("=" * 80)
    print("PROPER SHANNON ENTROPY COMPUTATION")
    print("HIV-1 gp120 relay pocket -- from NCBI + pairwise alignment to HXB2")
    print("=" * 80)

    print(f"\nHXB2 gp120 length: {len(HXB2_GP120)}")
    print(f"RPGGGDMR at seq position {RELAY_START} (paper convention: 465)")
    relay_from_hxb2 = HXB2_GP120[RELAY_START-1:RELAY_START+7]
    print(f"Verified: {relay_from_hxb2}")
    print(f"Offset to paper numbering: {PAPER_POS_OFFSET}")
    print()

    sequences = fetch_gp120_sequences()
    if not sequences:
        print("No sequences retrieved.")
        return

    # Align each sequence to HXB2 and extract positions
    print(f"\nAligning {len(sequences)} sequences to HXB2...")

    # Storage: position -> list of amino acids
    all_positions = {pos: [] for pos in range(1, len(HXB2_GP120) + 1)}
    aligned_count = 0
    failed_count = 0
    relay_exact = 0
    relay_variant = 0

    for i, (header, seq) in enumerate(sequences):
        if i > 0 and i % 100 == 0:
            print(f"  [{i}/{len(sequences)}] aligned={aligned_count}, failed={failed_count}")

        if len(seq) < 200 or len(seq) > 700:
            failed_count += 1
            continue

        try:
            mapping = align_to_hxb2(seq, HXB2_GP120)
        except Exception:
            failed_count += 1
            continue

        if mapping is None:
            failed_count += 1
            continue

        # Check if relay positions are mapped
        relay_mapped = sum(1 for p in RELAY_POSITIONS_HXB2 if p in mapping)
        if relay_mapped < 6:
            failed_count += 1
            continue

        aligned_count += 1

        # Record all positions
        for pos, aa in mapping.items():
            if pos in all_positions:
                all_positions[pos].append(aa)

        # Check relay sequence
        relay_seq = "".join(mapping.get(p, "?") for p in RELAY_POSITIONS_HXB2)
        if relay_seq == "RPGGGDMR":
            relay_exact += 1
        else:
            relay_variant += 1

    print(f"\n  Total aligned: {aligned_count}")
    print(f"  Failed/skipped: {failed_count}")
    print(f"  Exact RPGGGDMR: {relay_exact}")
    print(f"  Variants: {relay_variant}")

    if aligned_count < 20:
        print("  Too few sequences for reliable entropy computation!")
        return

    # Compute per-position entropy for relay pocket
    print(f"\n{'='*80}")
    print("RELAY POCKET ENTROPY (positions 465-472)")
    print(f"{'='*80}")
    print(f"{'Pos':>4}  {'HXB2':>4}  {'Entropy':>8}  {'Paper':>8}  {'Diff':>8}  {'N':>6}  {'Cons%':>6}  Top AAs")
    print("-" * 90)

    relay_results = {}
    for pos in RELAY_POSITIONS_HXB2:
        column = all_positions[pos]
        entropy = shannon_entropy(column)
        paper_pos = pos + PAPER_POS_OFFSET
        paper_val = PAPER_VALUES.get(paper_pos)
        diff = entropy - paper_val if paper_val is not None else None

        counts = Counter(column)
        total = len(column)
        top = counts.most_common(5)
        top_str = ", ".join(f"{aa}:{c}({100*c/total:.0f}%)" for aa, c in top)
        top1_pct = 100 * top[0][1] / total if top else 0

        hxb2_aa = HXB2_RELAY[pos]
        diff_str = f"{diff:+.4f}" if diff is not None else "N/A"
        paper_str = f"{paper_val:.2f}" if paper_val is not None else "N/A"

        print(f"  {paper_pos:>3}  {hxb2_aa:>4}  {entropy:>8.4f}  {paper_str:>8}  {diff_str:>8}  {total:>6}  {top1_pct:>5.1f}%  {top_str}")

        relay_results[str(paper_pos)] = {
            "hxb2_aa": hxb2_aa,
            "computed_entropy": round(entropy, 4),
            "paper_value": paper_val,
            "difference": round(diff, 4) if diff is not None else None,
            "n_sequences": total,
            "conservation_pct": round(top1_pct, 1),
            "aa_distribution": {aa: c for aa, c in counts.most_common()},
        }

    relay_entropies = [v["computed_entropy"] for v in relay_results.values()]
    computed_mean = sum(relay_entropies) / len(relay_entropies)
    print(f"\n  Computed relay mean: {computed_mean:.4f}")
    print(f"  Paper relay mean:   {PAPER_VALUES[465]}")

    # Compute domain-level entropy
    print(f"\n{'='*80}")
    print("DOMAIN COMPARISON")
    print(f"{'='*80}")

    # Domain boundaries use paper convention; adjust to our seq numbering
    off = -PAPER_POS_OFFSET
    domains = {
        "V1 (131-157)": range(131 + off, 158 + off),
        "V2 (158-196)": range(158 + off, 197 + off),
        "C2 (196-296)": range(196 + off, 297 + off),
        "V3 (296-331)": range(296 + off, 332 + off),
        "C3 (331-385)": range(331 + off, 386 + off),
        "V4 (385-418)": range(385 + off, 419 + off),
        "C4 (418-461)": range(418 + off, 462 + off),
        "C5 (461-485)": range(461 + off, 486 + off),
        "Relay (465-472)": range(RELAY_START, RELAY_START + 8),
    }

    domain_results = {}
    print(f"{'Domain':<20}  {'Computed':>10}  {'Paper':>10}  {'Diff':>10}  {'N_pos':>8}")
    print("-" * 65)

    for domain_name, pos_range in domains.items():
        domain_entropies = []
        for pos in pos_range:
            col = all_positions.get(pos, [])
            if len(col) >= 10:
                domain_entropies.append(shannon_entropy(col))

        if domain_entropies:
            domain_mean = sum(domain_entropies) / len(domain_entropies)
            paper_val = PAPER_DOMAIN_VALUES.get(domain_name)
            diff = domain_mean - paper_val if paper_val else None
            diff_str = f"{diff:+.4f}" if diff is not None else "N/A"
            paper_str = f"{paper_val:.2f}" if paper_val else "N/A"
            print(f"  {domain_name:<18}  {domain_mean:>10.4f}  {paper_str:>10}  {diff_str:>10}  {len(domain_entropies):>8}")
            domain_results[domain_name] = {
                "computed_mean": round(domain_mean, 4),
                "paper_value": paper_val,
                "n_positions": len(domain_entropies),
            }

    # Verdict
    max_diff = max(abs(v["difference"]) for v in relay_results.values() if v["difference"] is not None)
    if max_diff < 0.05:
        verdict = "MATCH -- paper values are accurate"
    elif max_diff < 0.15:
        verdict = "CLOSE -- paper values are approximate but directionally correct"
    else:
        verdict = "MISMATCH -- paper values need replacement"

    print(f"\n{'='*80}")
    print(f"VERDICT: {verdict}")
    print(f"Max per-residue difference: {max_diff:.4f}")
    print(f"{'='*80}")

    output = {
        "analysis": "Proper Shannon entropy computation -- NCBI HIV-1 gp120 + HXB2 pairwise alignment",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "method": "Biopython PairwiseAligner global alignment to HXB2 gp120",
        "n_sequences_fetched": len(sequences),
        "n_sequences_aligned": aligned_count,
        "n_exact_RPGGGDMR": relay_exact,
        "n_variants": relay_variant,
        "relay_pocket": relay_results,
        "domain_comparison": domain_results,
        "computed_relay_mean": round(computed_mean, 4),
        "paper_relay_mean": 0.15,
        "max_per_residue_diff": round(max_diff, 4),
        "verdict": verdict,
    }

    with open(OUTPUT_FILE, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
