#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse ConCavity results for relay pocket residues 465-472."""

relay_residues = {465, 466, 467, 468, 469, 470, 471, 472}
relay_names = {465:'ARG', 466:'PRO', 467:'GLY', 468:'GLY', 469:'GLY', 470:'ASP', 471:'MET', 472:'ARG'}

scores = {}
with open(os.path.expanduser('~/concavity_results/hiv_6meo_relay_pocket_residue.pdb')) as f:
    for line in f:
        if line.startswith('ATOM'):
            chain = line[21]
            if chain == 'G':
                resnum = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if atomname == 'CA':
                    scores[resnum] = (resname, bfactor)

print('=== ConCavity Binding Site Prediction -- Relay Pocket (465-472) ===')
print()
all_vals = sorted(scores.items(), key=lambda x: -x[1][1])
total = len(all_vals)

relay_scores = []
for resnum in sorted(relay_residues):
    if resnum in scores:
        resname, sc = scores[resnum]
        rank = next(i+1 for i, (rn, (nm, s)) in enumerate(all_vals) if rn == resnum)
        percentile = rank / total * 100
        print(f'  {resname}{resnum}: score={sc:.4f}  rank={rank}/{total}  top {percentile:.1f}%')
        relay_scores.append(sc)
    else:
        print(f'  {resnum}: NOT FOUND in chain G')

print()
print(f'=== Top 30 Highest-Scoring Residues (Chain G, {total} total) ===')
for i, (resnum, (resname, sc)) in enumerate(all_vals[:30]):
    marker = ' <<< RELAY POCKET' if resnum in relay_residues else ''
    print(f'  #{i+1:2d}  {resname}{resnum}: {sc:.4f}{marker}')

vals = [s for _, s in scores.values()]
mean_all = sum(vals)/len(vals)
median_all = sorted(vals)[len(vals)//2]

print(f'\n=== Statistics ===')
print(f'  All residues: mean={mean_all:.4f}, median={median_all:.4f}, max={max(vals):.4f}')
if relay_scores:
    mean_relay = sum(relay_scores)/len(relay_scores)
    print(f'  Relay pocket: mean={mean_relay:.4f}, max={max(relay_scores):.4f}, min={min(relay_scores):.4f}')
    print(f'  Relay vs all mean: {mean_relay/mean_all:.1f}x enrichment')
    above_median = sum(1 for s in relay_scores if s > median_all)
    above_75th = sum(1 for s in relay_scores if s > sorted(vals)[int(len(vals)*0.75)])
    above_90th = sum(1 for s in relay_scores if s > sorted(vals)[int(len(vals)*0.90)])
    print(f'  Relay above median: {above_median}/{len(relay_scores)}')
    print(f'  Relay above 75th percentile: {above_75th}/{len(relay_scores)}')
    print(f'  Relay above 90th percentile: {above_90th}/{len(relay_scores)}')

# Also check the contact shell / wall residues (from pharmacophore analysis)
wall_residues = [109, 112, 113, 200, 202, 252, 253, 256, 257, 375, 376, 377, 378, 382, 421, 422, 423, 425, 427, 429, 431, 433, 434, 435, 436, 473, 474, 475, 476]
print(f'\n=== Contact Shell / Wall Residues ({len(wall_residues)} residues) ===')
wall_scores = []
for resnum in sorted(wall_residues):
    if resnum in scores:
        resname, sc = scores[resnum]
        rank = next(i+1 for i, (rn, _) in enumerate(all_vals) if rn == resnum)
        wall_scores.append(sc)
if wall_scores:
    mean_wall = sum(wall_scores)/len(wall_scores)
    above_median_wall = sum(1 for s in wall_scores if s > median_all)
    print(f'  Wall residues found: {len(wall_scores)}/{len(wall_residues)}')
    print(f'  Wall mean score: {mean_wall:.4f} ({mean_wall/mean_all:.1f}x vs all)')
    print(f'  Wall above median: {above_median_wall}/{len(wall_scores)}')
