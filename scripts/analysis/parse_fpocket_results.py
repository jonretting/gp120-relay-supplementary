#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse fpocket results and identify which pockets contain relay pocket residues."""

import os, re, glob

relay_residues = {465, 466, 467, 468, 469, 470, 471, 472}
relay_names = {465:'ARG', 466:'PRO', 467:'GLY', 468:'GLY', 469:'GLY', 470:'ASP', 471:'MET', 472:'ARG'}

pocket_dir = os.path.expanduser("~/hiv_6meo_out/pockets")
pocket_files = sorted(glob.glob(os.path.join(pocket_dir, "pocket*_atm.pdb")))

print(f"=== fpocket 4.2 -- Relay Pocket Detection ===")
print(f"Total pockets detected: {len(pocket_files)}")
print()

relay_found_in = {}

for pf in pocket_files:
    pocket_num = int(re.search(r'pocket(\d+)', pf).group(1))
    residues_in_pocket = set()
    chains_in_pocket = set()
    
    with open(pf) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                residues_in_pocket.add((chain, resnum, resname))
                chains_in_pocket.add(chain)
    
    relay_in_this = set()
    for chain, resnum, resname in residues_in_pocket:
        if chain == 'G' and resnum in relay_residues:
            relay_in_this.add(resnum)
    
    if relay_in_this:
        relay_found_in[pocket_num] = {
            'relay_residues': relay_in_this,
            'total_residues': len(residues_in_pocket),
            'chains': chains_in_pocket,
        }
        print(f"  Pocket {pocket_num}: {len(relay_in_this)}/8 relay residues")
        for rn in sorted(relay_in_this):
            print(f"    {relay_names.get(rn, '???')}{rn}")
        print(f"    Total residues in pocket: {len(residues_in_pocket)}, chains: {','.join(sorted(chains_in_pocket))}")
        print()

# Summary
all_relay_found = set()
for info in relay_found_in.values():
    all_relay_found.update(info['relay_residues'])

print(f"=== Summary ===")
print(f"  Relay residues found in fpocket predictions: {len(all_relay_found)}/8")
missing = relay_residues - all_relay_found
if missing:
    print(f"  Missing: {', '.join(relay_names.get(r, f'?{r}') + str(r) for r in sorted(missing))}")
else:
    print(f"  ALL relay pocket residues independently detected by fpocket")
print(f"  Pockets containing relay residues: {sorted(relay_found_in.keys())}")

# Parse the info file for druggability of relevant pockets
print()
print(f"=== Druggability Scores for Relay-Containing Pockets ===")
info_file = os.path.expanduser("~/hiv_6meo_out/hiv_6meo_info.txt")
with open(info_file) as f:
    content = f.read()

for pocket_num in sorted(relay_found_in.keys()):
    pattern = rf"Pocket {pocket_num} :.*?(?=Pocket \d+ :|$)"
    match = re.search(pattern, content, re.DOTALL)
    if match:
        block = match.group()
        score_m = re.search(r'Score\s*:\s*([\d.]+)', block)
        drug_m = re.search(r'Druggability Score\s*:\s*([\d.]+)', block)
        vol_m = re.search(r'Volume\s*:\s*([\d.]+)', block)
        if score_m:
            print(f"  Pocket {pocket_num}: score={score_m.group(1)}, druggability={drug_m.group(1) if drug_m else 'N/A'}, volume={vol_m.group(1) if vol_m else 'N/A'} A^3")
