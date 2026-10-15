# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
import json, os

SUPP_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))

d = json.load(open(os.path.join(SUPP_DIR, "results", "03-conservation", "cross_clade_survey.json")))
print(f"Total structures: {d['total_structures']}")
print(f"With pocket: {d['with_relay_pocket']}")
print(f"Exact match: {d['exact_match']}")
print(f"Without: {d['without_pocket']}")
print(f"Errors: {d['errors']}")
print(f"Conservation: {d['conservation_rate']}%")

# Check what structures were without
for s in d.get("structures", []):
    if not s["has_relay_pocket"]:
        print(f"  NO POCKET: {s['pdb_id']} - {s['title'][:60]}")
