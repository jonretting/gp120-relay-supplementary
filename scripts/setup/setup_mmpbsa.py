import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Create index file and input for gmx_MMPBSA on gp120-CD4 interface."""
import subprocess
import os

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
os.makedirs(WORKDIR, exist_ok=True)

TPR = os.path.join(DATA_ROOT, "topology", "md.tpr")
XTC = os.path.join(DATA_ROOT, "topology", "md.xtc")  # Note: large trajectory not in package

# Create index groups: chain G (gp120) as receptor, chain A (CD4) as ligand
# From topology: Protein_chain_G is molblock 0, Protein_chain_A is molblock 1
# We need to create a custom index file with these groups

# Use gmx make_ndx to split by chain
# Chain G atoms come first, then chain A
# From GRO inspection: 82 gp120 residues, 23 CD4 residues = 105 total
# We need atom indices, not residue indices

# Read GRO to get atom counts per chain
with open(os.path.join(DATA_ROOT, "topology", "md.gro")) as f:
    lines = f.readlines()

chain_g_atoms = []  # gp120
chain_a_atoms = []  # CD4
in_chain_a = False
prev_resid = None

for line in lines[2:-1]:
    if len(line) < 20:
        continue
    resname = line[5:10].strip()
    if resname in ('SOL', 'NA', 'CL'):
        break
    try:
        resid = int(line[0:5].strip())
        atomnum = int(line[15:20].strip())
    except ValueError:
        continue

    # Chain boundary: resid goes from 478 -> 23 (chain A starts)
    if prev_resid is not None and resid < prev_resid and resid < 100:
        in_chain_a = True

    if in_chain_a:
        chain_a_atoms.append(atomnum)
    else:
        chain_g_atoms.append(atomnum)

    prev_resid = resid

print(f'Chain G (gp120): {len(chain_g_atoms)} atoms')
print(f'Chain A (CD4):   {len(chain_a_atoms)} atoms')
print(f'Total protein:   {len(chain_g_atoms) + len(chain_a_atoms)} atoms')

# Write index file
ndx_path = os.path.join(WORKDIR, 'mmpbsa.ndx')
with open(ndx_path, 'w') as f:
    f.write('[ receptor ]\n')
    for i, a in enumerate(chain_g_atoms):
        f.write(f'{a:6d}')
        if (i + 1) % 15 == 0:
            f.write('\n')
    f.write('\n')

    f.write('[ ligand ]\n')
    for i, a in enumerate(chain_a_atoms):
        f.write(f'{a:6d}')
        if (i + 1) % 15 == 0:
            f.write('\n')
    f.write('\n')

print(f'Index file written: {ndx_path}')

# Write gmx_MMPBSA input file
# GB model (igb=5) is faster and doesn't need a PDB radius file
# Per-residue decomposition is the key output
inp_path = os.path.join(WORKDIR, 'mmpbsa.in')
with open(inp_path, 'w') as f:
    f.write("""&general
startframe=1, endframe=9999, interval=10,
verbose=2, keep_files=0,
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=3,
print_res="within 6"
/
""")

print(f'Input file written: {inp_path}')
print()
print('Ready to run:')
print(f'  cd {WORKDIR}')
print(f'  gmx_MMPBSA -O -i mmpbsa.in -cs {TPR} -ci mmpbsa.ndx -cg 0 1 -ct {XTC} -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -do FINAL_DECOMP_MMPBSA.dat -deo FINAL_DECOMP_MMPBSA.csv')
