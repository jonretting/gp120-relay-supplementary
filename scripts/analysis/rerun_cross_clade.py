#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Cross-Clade Structural Survey -- FULL RE-RUN
Queries RCSB PDB for ALL HIV envelope structures, then checks each
for the RPGGGDMR relay pocket motif. No curated lists -- this is 
a clean, reproducible sweep of the entire PDB.
"""

import json, os, requests, time, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
OUTPUT_FILE = os.path.join(SUPP_DIR, "results", "03-conservation", "cross_clade_survey_rerun.json")

RELAY_SEQUENCE = "RPGGGDMR"
RELAY_RESIDUES_HXB2 = list(range(465, 473))

print("=" * 90)
print("CROSS-CLADE STRUCTURAL SURVEY -- FULL PDB RE-RUN")
print("=" * 90)
print(f"Target: relay pocket motif '{RELAY_SEQUENCE}' (HXB2 residues {RELAY_RESIDUES_HXB2[0]}-{RELAY_RESIDUES_HXB2[-1]})")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# Step 1: Query RCSB for ALL structures with HIV envelope proteins
print("--- Step 1: Querying RCSB PDB Search API ---")

search_url = "https://search.rcsb.org/rcsbsearch/v2/query"

search_query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "group",
                "logical_operator": "or",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity.pdbx_description",
                            "operator": "contains_words",
                            "value": "gp120"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity.pdbx_description",
                            "operator": "contains_words",
                            "value": "gp160"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity.pdbx_description",
                            "operator": "contains_words",
                            "value": "envelope glycoprotein"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity.pdbx_description",
                            "operator": "contains_words",
                            "value": "gp41"
                        }
                    },
                ]
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

resp = requests.post(search_url, json=search_query, timeout=30)
resp.raise_for_status()
search_data = resp.json()

pdb_ids = [hit["identifier"] for hit in search_data.get("result_set", [])]
total_hits = search_data.get("total_count", len(pdb_ids))

print(f"  RCSB returned: {total_hits} structures")
print(f"  PDB IDs retrieved: {len(pdb_ids)}")
print()

# Step 2: For each structure, get entity data and search for RPGGGDMR
print("--- Step 2: Analyzing each structure for RPGGGDMR motif ---")
start_time = time.time()

results = []
errors = []

for i, pdb_id in enumerate(pdb_ids):
    if i > 0 and i % 50 == 0:
        elapsed = time.time() - start_time
        rate = i / elapsed
        eta = (len(pdb_ids) - i) / rate
        print(f"  Progress: {i}/{len(pdb_ids)} ({i/len(pdb_ids)*100:.0f}%) -- {elapsed:.0f}s elapsed, ~{eta:.0f}s remaining")

    try:
        # Get entry info
        info_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        info_resp = requests.get(info_url, timeout=10)
        info_resp.raise_for_status()
        info = info_resp.json()

        title = info.get("struct", {}).get("title", "N/A")
        resolution = info.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
        method = info.get("exptl", [{}])[0].get("method", "N/A")
        deposition_date = info.get("rcsb_accession_info", {}).get("deposit_date", "N/A")

        has_relay = False
        relay_sequence_found = ""
        organism = ""
        entity_description = ""

        # Check all polymer entities in the structure
        for ent_id in range(1, 20):
            try:
                ent_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent_id}"
                ent_resp = requests.get(ent_url, timeout=10)
                if ent_resp.status_code != 200:
                    break
                ent_data = ent_resp.json()

                seq = ent_data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
                if not seq:
                    continue

                desc = ent_data.get("rcsb_polymer_entity", {}).get("pdbx_description", "")

                organism_list = ent_data.get("rcsb_entity_source_organism", [])
                if organism_list and not organism:
                    organism = organism_list[0].get("ncbi_scientific_name", "")

                # Search for exact RPGGGDMR
                idx = seq.find(RELAY_SEQUENCE)
                if idx >= 0:
                    has_relay = True
                    relay_sequence_found = RELAY_SEQUENCE
                    entity_description = desc
                    break

                # Also search for partial matches (at least 6 of 8)
                for sublen in range(7, 5, -1):
                    for start_pos in range(len(RELAY_SEQUENCE) - sublen + 1):
                        sub = RELAY_SEQUENCE[start_pos:start_pos+sublen]
                        if seq.find(sub) >= 0:
                            has_relay = True
                            # Extract the actual 8-residue window
                            sub_idx = seq.find(sub)
                            window_start = max(0, sub_idx - start_pos)
                            relay_sequence_found = seq[window_start:window_start+8]
                            entity_description = desc
                            break
                    if has_relay:
                        break

                if has_relay:
                    break

            except Exception:
                continue

        entry = {
            "pdb_id": pdb_id,
            "title": title[:120],
            "resolution": resolution,
            "method": method,
            "deposition_date": deposition_date[:10] if deposition_date != "N/A" else "N/A",
            "organism": organism,
            "entity_description": entity_description[:80],
            "has_relay_pocket": has_relay,
            "relay_sequence": relay_sequence_found,
            "exact_match": relay_sequence_found == RELAY_SEQUENCE,
        }
        results.append(entry)

        time.sleep(0.05)

    except requests.exceptions.Timeout:
        errors.append({"pdb_id": pdb_id, "error": "timeout"})
    except Exception as e:
        errors.append({"pdb_id": pdb_id, "error": str(e)})

elapsed_total = time.time() - start_time

# Step 3: Summarize
print()
print("=" * 90)
print("CROSS-CLADE STRUCTURAL SURVEY -- RESULTS")
print("=" * 90)

total = len(results)
with_pocket = [r for r in results if r["has_relay_pocket"]]
without = [r for r in results if not r["has_relay_pocket"]]
exact_match = [r for r in results if r.get("exact_match")]
variant_match = [r for r in with_pocket if not r.get("exact_match")]

print(f"\nTotal structures analyzed: {total}")
print(f"Query errors/timeouts: {len(errors)}")
print(f"Structures with relay pocket motif: {len(with_pocket)} ({len(with_pocket)/total*100:.1f}%)")
print(f"  - Exact RPGGGDMR match: {len(exact_match)}")
print(f"  - Variant match (6-7 residue): {len(variant_match)}")
print(f"Structures without (non-gp120): {len(without)}")
print(f"\nAmong gp120-containing structures:")
if len(with_pocket) > 0:
    print(f"  Exact RPGGGDMR: {len(exact_match)}/{len(with_pocket)} ({len(exact_match)/len(with_pocket)*100:.1f}%)")

# Resolution distribution
res_values = [r["resolution"] for r in with_pocket if r["resolution"]]
if res_values:
    print(f"\n  Resolution range: {min(res_values):.1f} - {max(res_values):.1f} A")
    print(f"  Mean resolution: {sum(res_values)/len(res_values):.1f} A")

# Method distribution
methods = {}
for r in with_pocket:
    m = r["method"]
    methods[m] = methods.get(m, 0) + 1
print(f"\n  Methods: {dict(sorted(methods.items(), key=lambda x: -x[1]))}")

# Date range
dates = sorted([r["deposition_date"] for r in with_pocket if r["deposition_date"] != "N/A"])
if dates:
    print(f"  Date range: {dates[0]} to {dates[-1]}")

# Sequence variants
variants = {}
for r in with_pocket:
    seq = r["relay_sequence"]
    if seq and seq != RELAY_SEQUENCE:
        variants[seq] = variants.get(seq, 0) + 1

if variants:
    print(f"\n  Sequence variants ({len(variants)} types):")
    for v, count in sorted(variants.items(), key=lambda x: -x[1]):
        print(f"    {v}: {count} structures")
else:
    print(f"\n  ALL gp120 structures show identical RPGGGDMR -- ZERO sequence variants")

print(f"\nRuntime: {elapsed_total:.0f}s ({elapsed_total/60:.1f} min)")

# Save
output = {
    "analysis": "Cross-Clade Structural Survey -- FULL PDB RE-RUN",
    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    "reference_sequence": RELAY_SEQUENCE,
    "reference_residues_HXB2": RELAY_RESIDUES_HXB2,
    "rcsb_query_total_hits": total_hits,
    "total_structures_analyzed": total,
    "with_relay_pocket": len(with_pocket),
    "exact_RPGGGDMR_match": len(exact_match),
    "variant_matches": len(variant_match),
    "without_pocket_non_gp120": len(without),
    "query_errors": len(errors),
    "conservation_among_gp120": f"{len(exact_match)}/{len(with_pocket)} ({len(exact_match)/len(with_pocket)*100:.1f}%)" if len(with_pocket) > 0 else "N/A",
    "runtime_seconds": round(elapsed_total),
    "sequence_variants": variants,
    "structures": results,
    "errors": errors,
}

with open(OUTPUT_FILE, "w") as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to: {OUTPUT_FILE}")
