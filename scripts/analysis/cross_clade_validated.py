#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Cross-Clade Structural Survey -- VALIDATED
Authoritative RCSB PDB survey for RPGGGDMR relay pocket motif in HIV-1 gp120.

Query strategy:
  - RCSB Search API v2: entity description contains "gp120" AND organism
    is "Human immunodeficiency virus"  (same query that originally produced
    the ~1,982-hit result set)
  - For each PDB entry, iterate polymer entities and pull the canonical
    one-letter sequence; search for exact RPGGGDMR and 7-of-8 partial matches
  - Record hit/miss, variant sequences, clade, resolution, method, date
  - Output: validated JSON with full audit trail

Runtime: ~30-60 min depending on RCSB latency (0.05s rate limit per request)
"""

import json, os, requests, time, sys, hashlib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
OUTPUT_FILE = os.path.join(SUPP_DIR, "results", "03-conservation", "cross_clade_survey_validated.json")

RELAY_SEQUENCE = "RPGGGDMR"
RELAY_RESIDUES_HXB2 = list(range(465, 473))

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent_id}"

RATE_LIMIT_S = 0.05


def build_search_query():
    """HIV-1 gp120 structures only -- no gp41, no gp160, no generic 'envelope'."""
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity.pdbx_description",
                        "operator": "contains_words",
                        "value": "gp120",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                        "operator": "contains_words",
                        "value": "Human immunodeficiency virus",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "return_all_hits": True,
        },
    }


def search_rcsb():
    """Return list of PDB IDs from RCSB search."""
    query = build_search_query()
    resp = requests.post(SEARCH_URL, json=query, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    pdb_ids = [hit["identifier"] for hit in data.get("result_set", [])]
    total_hits = data.get("total_count", len(pdb_ids))
    return pdb_ids, total_hits, query


def find_motif_in_sequence(seq, motif=RELAY_SEQUENCE):
    """Search for exact 8-mer or partial (7-of-8) matches."""
    idx = seq.find(motif)
    if idx >= 0:
        return motif, True, idx

    for start in range(len(motif) - 6):
        sub7 = motif[start : start + 7]
        idx = seq.find(sub7)
        if idx >= 0:
            window_start = max(0, idx - start)
            window = seq[window_start : window_start + 8]
            if len(window) == 8:
                return window, False, window_start
    return "", False, -1


def analyze_entry(pdb_id, session):
    """Fetch metadata + sequences for one PDB entry."""
    info_resp = session.get(ENTRY_URL.format(pdb_id=pdb_id), timeout=10)
    info_resp.raise_for_status()
    info = info_resp.json()

    title = info.get("struct", {}).get("title", "N/A")
    resolution = info.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
    method = info.get("exptl", [{}])[0].get("method", "N/A")
    dep_raw = info.get("rcsb_accession_info", {}).get("deposit_date", "N/A")
    deposition_date = dep_raw[:10] if dep_raw != "N/A" else "N/A"

    best_hit = None
    organism = ""

    for ent_id in range(1, 20):
        ent_resp = session.get(
            ENTITY_URL.format(pdb_id=pdb_id, ent_id=ent_id), timeout=10
        )
        if ent_resp.status_code != 200:
            break

        ent_data = ent_resp.json()
        desc = ent_data.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        seq = ent_data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")

        org_list = ent_data.get("rcsb_entity_source_organism", [])
        if org_list and not organism:
            organism = org_list[0].get("ncbi_scientific_name", "")

        if not seq:
            continue

        motif_found, is_exact, pos = find_motif_in_sequence(seq)
        if motif_found:
            if best_hit is None or (is_exact and not best_hit["exact_match"]):
                best_hit = {
                    "relay_sequence": motif_found,
                    "exact_match": is_exact,
                    "entity_id": ent_id,
                    "entity_description": desc[:100],
                    "seq_position": pos,
                }
            if is_exact:
                break

        time.sleep(RATE_LIMIT_S)

    entry = {
        "pdb_id": pdb_id,
        "title": title[:120],
        "resolution": resolution,
        "method": method,
        "deposition_date": deposition_date,
        "organism": organism,
        "has_relay_pocket": best_hit is not None,
        "relay_sequence": best_hit["relay_sequence"] if best_hit else "",
        "exact_match": best_hit["exact_match"] if best_hit else False,
        "entity_description": best_hit["entity_description"] if best_hit else "",
    }
    return entry


def main():
    start_wall = time.time()
    ts = time.strftime("%Y-%m-%d %H:%M:%S")

    print("=" * 90)
    print("CROSS-CLADE STRUCTURAL SURVEY -- VALIDATED")
    print("=" * 90)
    print(f"Target motif: {RELAY_SEQUENCE} (HXB2 {RELAY_RESIDUES_HXB2[0]}-{RELAY_RESIDUES_HXB2[-1]})")
    print(f"Started: {ts}")
    print(f"Query: entity description contains 'gp120' AND organism 'Human immunodeficiency virus'")
    print()

    pdb_ids, total_hits, query = search_rcsb()
    print(f"RCSB total_count: {total_hits}")
    print(f"PDB IDs retrieved: {len(pdb_ids)}")
    print()

    session = requests.Session()
    results = []
    errors = []

    for i, pdb_id in enumerate(pdb_ids):
        if i > 0 and i % 50 == 0:
            elapsed = time.time() - start_wall
            rate = i / elapsed if elapsed > 0 else 1
            eta = (len(pdb_ids) - i) / rate if rate > 0 else 0
            n_ok = sum(1 for r in results if r["has_relay_pocket"])
            print(
                f"  [{i}/{len(pdb_ids)}] {i/len(pdb_ids)*100:.0f}% -- "
                f"pocket hits so far: {n_ok} -- "
                f"elapsed {elapsed:.0f}s, ETA ~{eta:.0f}s"
            )

        try:
            entry = analyze_entry(pdb_id, session)
            results.append(entry)
        except requests.exceptions.Timeout:
            errors.append({"pdb_id": pdb_id, "error": "timeout"})
        except Exception as e:
            errors.append({"pdb_id": pdb_id, "error": str(e)[:200]})

        time.sleep(RATE_LIMIT_S)

    elapsed_total = time.time() - start_wall

    with_pocket = [r for r in results if r["has_relay_pocket"]]
    exact_matches = [r for r in results if r.get("exact_match")]
    variant_matches = [r for r in with_pocket if not r.get("exact_match")]

    variants = {}
    for r in variant_matches:
        seq = r["relay_sequence"]
        if seq:
            variants[seq] = variants.get(seq, 0) + 1

    methods = {}
    for r in with_pocket:
        m = r["method"]
        methods[m] = methods.get(m, 0) + 1

    res_vals = [r["resolution"] for r in with_pocket if r["resolution"]]
    dates = sorted(
        [r["deposition_date"] for r in with_pocket if r["deposition_date"] != "N/A"]
    )

    print()
    print("=" * 90)
    print("RESULTS")
    print("=" * 90)
    print(f"Total structures analyzed:     {len(results)}")
    print(f"Query errors/timeouts:         {len(errors)}")
    print(f"Structures with relay motif:   {len(with_pocket)}")
    print(f"  Exact RPGGGDMR:              {len(exact_matches)}")
    print(f"  Variant (7-of-8):            {len(variant_matches)}")
    print(f"Without motif (non-gp120/trimmed): {len(results) - len(with_pocket)}")
    print()
    if with_pocket:
        print(
            f"Conservation among gp120:      "
            f"{len(exact_matches)}/{len(with_pocket)} "
            f"({len(exact_matches)/len(with_pocket)*100:.1f}%)"
        )
    if res_vals:
        print(f"Resolution range:              {min(res_vals):.1f} - {max(res_vals):.1f} A")
    if dates:
        print(f"Deposition range:              {dates[0]} to {dates[-1]}")
    print(f"Methods:                       {dict(sorted(methods.items(), key=lambda x: -x[1]))}")
    if variants:
        print(f"Sequence variants:")
        for v, c in sorted(variants.items(), key=lambda x: -x[1]):
            print(f"  {v}: {c} structures")
    print(f"\nRuntime: {elapsed_total:.0f}s ({elapsed_total/60:.1f} min)")

    if errors:
        print(f"\nErrors ({len(errors)}):")
        for e in errors[:10]:
            print(f"  {e['pdb_id']}: {e['error'][:80]}")
        if len(errors) > 10:
            print(f"  ... and {len(errors)-10} more")

    query_hash = hashlib.sha256(json.dumps(query, sort_keys=True).encode()).hexdigest()[:16]

    output = {
        "analysis": "Cross-Clade Structural Survey -- VALIDATED",
        "version": "3.0",
        "timestamp": ts,
        "script": "cross_clade_validated.py",
        "query_description": "entity description contains 'gp120' AND organism 'Human immunodeficiency virus'",
        "query_hash_sha256_16": query_hash,
        "reference_sequence": RELAY_SEQUENCE,
        "reference_residues_HXB2": RELAY_RESIDUES_HXB2,
        "rcsb_query_total_hits": total_hits,
        "total_structures_analyzed": len(results),
        "with_relay_pocket": len(with_pocket),
        "exact_RPGGGDMR_match": len(exact_matches),
        "variant_matches": len(variant_matches),
        "without_pocket": len(results) - len(with_pocket),
        "query_errors": len(errors),
        "conservation_among_gp120": (
            f"{len(exact_matches)}/{len(with_pocket)} "
            f"({len(exact_matches)/len(with_pocket)*100:.1f}%)"
            if with_pocket
            else "N/A"
        ),
        "resolution_range_angstrom": (
            [round(min(res_vals), 1), round(max(res_vals), 1)] if res_vals else None
        ),
        "deposition_date_range": [dates[0], dates[-1]] if dates else None,
        "method_distribution": methods,
        "sequence_variants": variants,
        "runtime_seconds": round(elapsed_total),
        "structures": results,
        "errors": errors,
    }

    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
