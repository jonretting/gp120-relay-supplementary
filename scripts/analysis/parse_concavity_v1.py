#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse ConCavity results and map to relay pocket residues."""

relay_residues_pdb = {109, 212, 255, 271, 274, 426, 428, 432}

scores = {}
with open(os.path.expanduser('~/concavity_results/hiv_6meo_relay_pocket_residue.pdb')) as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21]
            if chain == 'G':
                resnum = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if atomname == 'CA':
                    scores[resnum] = (resname, bfactor)

print('=== Relay Pocket Residues (Chain G gp120) ===')
all_vals = sorted(scores.items(), key=lambda x: -x[1][1])
total = len(all_vals)
for resnum in sorted(relay_residues_pdb):
    if resnum in scores:
        resname, sc = scores[resnum]
        rank = next(i+1 for i, (rn, (nm, s)) in enumerate(all_vals) if rn == resnum)
        percentile = (1 - rank/total) * 100
        print(f'  {resname}{resnum}: score={sc:.4f}  rank={rank}/{total}  top {100-percentile:.1f}%')
    else:
        print(f'  ?{resnum}: NOT FOUND in chain G')

print()
print('=== Top 30 Highest-Scoring Residues (Chain G) ===')
for resnum, (resname, sc) in all_vals[:30]:
    marker = ' <-- RELAY POCKET' if resnum in relay_residues_pdb else ''
    print(f'  {resname}{resnum}: {sc:.4f}{marker}')

vals = [s for _, s in scores.values()]
print(f'\nTotal residues: {len(scores)}')
print(f'Mean score: {sum(vals)/len(vals):.4f}')
print(f'Max score: {max(vals):.4f}')
print(f'Median score: {sorted(vals)[len(vals)//2]:.4f}')

relay_scores = [scores[r][1] for r in relay_residues_pdb if r in scores]
if relay_scores:
    print(f'\nRelay pocket mean: {sum(relay_scores)/len(relay_scores):.4f}')
    print(f'Relay pocket max: {max(relay_scores):.4f}')
    print(f'Relay pocket min: {min(relay_scores):.4f}')
    above_median = sum(1 for s in relay_scores if s > sorted(vals)[len(vals)//2])
    print(f'Relay residues above median: {above_median}/{len(relay_scores)}')
