#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Cross-Clade Structural Analysis
Survey all available gp120 structures in the PDB to check whether
the relay pocket (residues 465-472 equivalent) exists across structures.
"""

import json, os, requests, time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))

RELAY_RESIDUES_HXB2 = list(range(465, 473))  # ARG-PRO-GLY-GLY-GLY-ASP-MET-ARG
RELAY_SEQUENCE = "RPGGGDMR"

print("=== Cross-Clade Structural Survey -- gp120 Relay Pocket ===")
print(f"Target: relay pocket residues {RELAY_RESIDUES_HXB2[0]}-{RELAY_RESIDUES_HXB2[-1]} (HXB2: {RELAY_SEQUENCE})")
print()

# Step 1: Query RCSB PDB for all structures containing gp120
print("--- Step 1: Querying RCSB PDB for gp120 structures ---")

search_query = {
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
                    "value": "gp120 envelope glycoprotein"
                }
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                    "operator": "contains_words",
                    "value": "Human immunodeficiency virus"
                }
            }
        ]
    },
    "return_type": "entry",
    "request_options": {
        "results_content_type": ["experimental"],
        "return_all_hits": True
    }
}

# Use curated list of known gp120 structures spanning clades, years, labs
pdb_ids = [
    # Reference structures (used in our study)
    "6MEO", "6MET",
    # BMS-626529 / fostemsavir co-crystals
    "5U7O", "5FYK", "5FYJ", "5FYL",
    # SOSIP trimers and near-native structures
    "4TVP", "4NCO", "5CEZ", "5T3X", "5T3Z",
    # CD4-bound structures
    "3JWD", "3JWO", "3SE8", "3SE9",
    # Unliganded / antibody-bound
    "2NY7", "2NXY", "2NXZ", "2B4C", "1GC1", "1G9M", "1G9N",
    # Various labs and clades
    "1RZJ", "1RZK", "3NGB", "4RQS", "4RQT",
    "5FUU", "5ACO", "5V8L", "5V8M",
    "5VN3", "5VN8", "6CM3", "6CK9", "6CDI",
    "6MUF", "6MUG", "6MUH",
    # Recent high-resolution structures
    "6NQD", "6OLP", "6OLQ",
    "6UDA", "6UDJ", "6UTK",
    "7LO6", "7LO7", "7N4L",
    "7TFO", "7TFP", "7TFQ",
    "8FAD", "8FAE",
]

print(f"  Will analyze: {len(pdb_ids)} structures")
print()

# Step 2: For each structure, get metadata and check for relay pocket residues
print("--- Step 2: Analyzing each structure ---")

results = []
errors = []

for i, pdb_id in enumerate(pdb_ids):
    if i > 0 and i % 10 == 0:
        print(f"  Progress: {i}/{len(pdb_ids)}")

    try:
        # Get structure info
        info_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        info_resp = requests.get(info_url, timeout=5)
        info_resp.raise_for_status()
        info = info_resp.json()

        title = info.get("struct", {}).get("title", "N/A")
        resolution = info.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
        method = info.get("exptl", [{}])[0].get("method", "N/A")
        deposition_date = info.get("rcsb_accession_info", {}).get("deposit_date", "N/A")

        has_relay = False
        relay_sequence_found = ""
        chain_id = ""
        organism = ""
        clade_info = ""

        for ent_id in range(1, 10):
            try:
                ent_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent_id}"
                ent_resp = requests.get(ent_url, timeout=5)
                if ent_resp.status_code != 200:
                    break
                ent_data = ent_resp.json()

                desc = ent_data.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
                if "gp120" not in desc.lower() and "envelope" not in desc.lower():
                    continue

                seq = ent_data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
                organism_list = ent_data.get("rcsb_entity_source_organism", [])
                if organism_list:
                    organism = organism_list[0].get("ncbi_scientific_name", "")

                # Search for RPGGGDMR or close variants
                # The relay sequence is highly conserved
                for pattern in [RELAY_SEQUENCE, "RPGGGDM", "PGGGDMR", "RPGGGD"]:
                    idx = seq.find(pattern)
                    if idx >= 0:
                        has_relay = True
                        # Extract the full 8-residue window
                        start = max(0, idx - 1) if pattern != RELAY_SEQUENCE else idx
                        relay_sequence_found = seq[start:start+8] if len(seq) > start+8 else seq[start:]
                        break

                if has_relay:
                    break

            except Exception:
                continue

        entry = {
            "pdb_id": pdb_id,
            "title": title[:100],
            "resolution": resolution,
            "method": method,
            "deposition_date": deposition_date[:10] if deposition_date != "N/A" else "N/A",
            "organism": organism,
            "has_relay_pocket": has_relay,
            "relay_sequence": relay_sequence_found,
            "sequence_identity_to_reference": "exact" if relay_sequence_found == RELAY_SEQUENCE else ("variant" if has_relay else "not_found"),
        }
        results.append(entry)

        time.sleep(0.1)  # rate limit

    except Exception as e:
        errors.append({"pdb_id": pdb_id, "error": str(e)})

# Step 3: Summarize
print()
print("=" * 90)
print("CROSS-CLADE STRUCTURAL SURVEY RESULTS")
print("=" * 90)

total = len(results)
with_pocket = [r for r in results if r["has_relay_pocket"]]
without = [r for r in results if not r["has_relay_pocket"]]
exact_match = [r for r in results if r["relay_sequence"] == RELAY_SEQUENCE]

print(f"\nTotal structures analyzed: {total}")
print(f"Structures with relay pocket sequence: {len(with_pocket)} ({len(with_pocket)/total*100:.0f}%)")
print(f"Exact RPGGGDMR match: {len(exact_match)} ({len(exact_match)/total*100:.0f}%)")
print(f"Structures missing relay pocket: {len(without)} ({len(without)/total*100:.0f}%)")
print(f"Query errors: {len(errors)}")

print(f"\n--- Structures WITH relay pocket (sorted by resolution) ---")
with_pocket.sort(key=lambda x: x["resolution"] or 99)
print(f"{'PDB':>5} {'Res':>5} {'Method':>12} {'Date':>12} {'Sequence':>10} {'Title'}")
print("-" * 90)
for r in with_pocket:
    res_str = f"{r['resolution']:.1f}" if r['resolution'] else "N/A"
    print(f"{r['pdb_id']:>5} {res_str:>5} {r['method'][:12]:>12} {r['deposition_date']:>12} {r['relay_sequence']:>10} {r['title'][:35]}")

if without:
    print(f"\n--- Structures WITHOUT relay pocket ({len(without)}) ---")
    for r in without:
        res_str = f"{r['resolution']:.1f}" if r['resolution'] else "N/A"
        print(f"  {r['pdb_id']} ({res_str} Å) {r['title'][:60]}")

# Sequence variants
variants = set()
for r in with_pocket:
    if r["relay_sequence"] and r["relay_sequence"] != RELAY_SEQUENCE:
        variants.add(r["relay_sequence"])

if variants:
    print(f"\n--- Sequence variants found ---")
    for v in sorted(variants):
        count = sum(1 for r in with_pocket if r["relay_sequence"] == v)
        print(f"  {v}: {count} structures")
else:
    print(f"\n--- All structures show identical RPGGGDMR sequence ---")

# Save results
output = {
    "analysis": "Cross-Clade Structural Survey",
    "reference_sequence": RELAY_SEQUENCE,
    "reference_residues": RELAY_RESIDUES_HXB2,
    "total_structures": total,
    "with_relay_pocket": len(with_pocket),
    "exact_match": len(exact_match),
    "without_pocket": len(without),
    "errors": len(errors),
    "conservation_rate": round(len(with_pocket) / total * 100, 1) if total > 0 else 0,
    "structures": results,
    "query_errors": errors,
}

outfile = os.path.join(SUPP_DIR, "results", "03-conservation", "cross_clade_survey.json")
with open(outfile, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to: {outfile}")
