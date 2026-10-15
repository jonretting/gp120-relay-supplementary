#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Check 6MEO PDB residue numbering for chain G (gp120)."""

# Find the RPGGGDMR motif in chain G
target_seq = "RPGGGDMR"
residues = {}

with open(os.path.expanduser('~/hiv_6meo.pdb')) as f:
    for line in f:
        if line.startswith('ATOM'):
            chain = line[21]
            if chain == 'G':
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if atomname == 'CA':
                    residues[resnum] = resname

one_letter = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

sorted_res = sorted(residues.items())
seq = ''.join(one_letter.get(residues[r], '?') for r, _ in sorted_res)

print(f"Chain G total residues: {len(sorted_res)}")
print(f"Residue range: {sorted_res[0][0]} - {sorted_res[-1][0]}")
print()

# Find RPGGGDMR in the sequence
idx = seq.find(target_seq)
if idx >= 0:
    start_resnum = sorted_res[idx][0]
    end_resnum = sorted_res[idx + len(target_seq) - 1][0]
    print(f"RPGGGDMR found at position {idx} in chain G sequence")
    print(f"PDB residue numbers: {start_resnum} - {end_resnum}")
    print()
    for i, aa in enumerate(target_seq):
        rn = sorted_res[idx + i][0]
        rname = sorted_res[idx + i][1]
        print(f"  {aa} = {rname}{rn}")
else:
    print(f"RPGGGDMR NOT found in chain G!")
    print("Searching for partial matches...")
    for sublen in range(len(target_seq)-1, 3, -1):
        for start in range(len(target_seq) - sublen + 1):
            sub = target_seq[start:start+sublen]
            idx = seq.find(sub)
            if idx >= 0:
                print(f"  Found '{sub}' at position {idx}, PDB resnum {sorted_res[idx][0]}")
                break

# Also print residues 460-480 range
print()
print("=== PDB residues near expected relay pocket ===")
for rn, rname in sorted_res:
    if 460 <= rn <= 480 or 370 <= rn <= 400:
        ol = one_letter.get(rname, '?')
        print(f"  {rname}{rn} ({ol})")
