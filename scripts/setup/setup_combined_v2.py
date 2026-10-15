import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Set up combined ternary+ligand system using the correctly positioned
GAFF2-parameterized BMS-626529 coordinates."""
import os, shutil

WORKDIR = os.path.join(DATA_ROOT, "results", "04-md-simulations")
TOPOLOGY_DIR = os.path.join(WORKDIR, 'topology')
LIGAND_DIR = os.path.join(WORKDIR, 'positioned_acpype')
TERNARY_GRO = os.path.join(DATA_ROOT, "topology", "md.gro")

os.chdir(WORKDIR)

# Step 1: Extract ligand ITP from positioned_acpype GMX top
print('=== Step 1: Extract ligand ITP ===')
with open(os.path.join(LIGAND_DIR, 'bms626529_GMX.top')) as f:
    top_content = f.read()

lines = top_content.split('\n')
in_atomtypes = False
in_moleculetype = False
atomtypes_lines = []
itp_lines = []

for line in lines:
    stripped = line.strip()
    if stripped.startswith('[ defaults ]'):
        continue
    if stripped.startswith('[ atomtypes ]'):
        in_atomtypes = True
        atomtypes_lines.append(line)
        continue
    if stripped.startswith('[ moleculetype ]'):
        in_atomtypes = False
        in_moleculetype = True
        itp_lines.append(line)
        continue
    if stripped.startswith('[ system ]') or stripped.startswith('[ molecules ]'):
        in_moleculetype = False
        continue
    if in_atomtypes:
        atomtypes_lines.append(line)
    elif in_moleculetype:
        itp_lines.append(line)

with open(os.path.join(WORKDIR, 'UNL_atomtypes.itp'), 'w') as f:
    f.write('; BMS-626529 atom types - GAFF2\n')
    f.write('\n'.join(atomtypes_lines))
    f.write('\n')

with open(os.path.join(WORKDIR, 'UNL.itp'), 'w') as f:
    f.write('; BMS-626529 topology - GAFF2 with AM1-BCC charges\n\n')
    f.write('\n'.join(itp_lines))
    f.write('\n')

print(f'  UNL_atomtypes.itp: {len(atomtypes_lines)} lines')
print(f'  UNL.itp: {len(itp_lines)} lines')

# Step 2: Read ligand GRO (in ternary frame)
print('\n=== Step 2: Read ligand coordinates ===')
with open(os.path.join(LIGAND_DIR, 'bms626529_GMX.gro')) as f:
    lig_lines = f.readlines()

lig_atoms = lig_lines[2:-1]
n_lig = int(lig_lines[1].strip())
print(f'  Ligand: {n_lig} atoms')

# Verify coordinates are in ternary frame
x_vals = []
for line in lig_atoms:
    x = float(line[20:28])
    x_vals.append(x)
x_min, x_max = min(x_vals), max(x_vals)
print(f'  X range: {x_min:.2f} - {x_max:.2f} nm ({x_min*10:.1f} - {x_max*10:.1f} A)')

# Step 3: Read ternary structure
print('\n=== Step 3: Read ternary structure ===')
with open(TERNARY_GRO) as f:
    ternary_lines = f.readlines()

n_ternary = int(ternary_lines[1].strip())
box_line = ternary_lines[-1]
ternary_atoms = ternary_lines[2:-1]
print(f'  Ternary: {n_ternary} atoms')

# Step 4: Remove overlapping waters
print('\n=== Step 4: Remove overlapping waters ===')
lig_positions = []
for line in lig_atoms:
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    lig_positions.append((x, y, z))

OVERLAP_DIST = 0.2  # nm
filtered_atoms = []
removed_waters = 0
i = 0
while i < len(ternary_atoms):
    line = ternary_atoms[i]
    resname = line[5:10].strip()
    if resname == 'SOL':
        try:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except (ValueError, IndexError):
            filtered_atoms.append(line)
            i += 1
            continue
        overlaps = False
        for lx, ly, lz in lig_positions:
            d2 = (x-lx)**2 + (y-ly)**2 + (z-lz)**2
            if d2 < OVERLAP_DIST**2:
                overlaps = True
                break
        if overlaps:
            removed_waters += 1
            i += 3
            continue
    filtered_atoms.append(line)
    i += 1

print(f'  Removed {removed_waters} overlapping waters ({removed_waters*3} atoms)')

# Insert ligand after protein, before SOL
water_start = None
for idx, line in enumerate(filtered_atoms):
    if line[5:10].strip() == 'SOL':
        water_start = idx
        break
if water_start is None:
    water_start = len(filtered_atoms)

combined = filtered_atoms[:water_start] + lig_atoms + filtered_atoms[water_start:]
n_total = len(combined)
print(f'  Combined: {n_total} atoms')

with open('combined.gro', 'w') as f:
    f.write('Ternary + BMS-626529\n')
    f.write(f'{n_total}\n')
    for line in combined:
        f.write(line if line.endswith('\n') else line + '\n')
    f.write(box_line)

# Step 5: Create topology
print('\n=== Step 5: Write combined topology ===')
with open(os.path.join(TOPOLOGY_DIR, 'topol.top')) as f:
    top_lines = f.readlines()

new_top = []
for line in top_lines:
    if '#include "amber99sb-ildn.ff/forcefield.itp"' in line:
        new_top.append(line)
        new_top.append('\n; BMS-626529 ligand atom types\n')
        new_top.append('#include "UNL_atomtypes.itp"\n')
        continue
    if '#include "topol_Protein_chain_C.itp"' in line:
        new_top.append(line.replace('#include "', '#include "topology/'))
        new_top.append('\n; BMS-626529 ligand topology\n')
        new_top.append('#include "UNL.itp"\n')
        continue
    if '#include "topol_Protein_chain_' in line or '#include "posre_Protein_chain_' in line:
        new_top.append(line.replace('#include "', '#include "topology/'))
        continue
    if line.startswith('SOL'):
        sol_count = 92689 - removed_waters
        new_top.append(f'SOL         {sol_count}\n')
        continue
    new_top.append(line)

# Insert UNL molecule entry before SOL
final_top = []
for line in new_top:
    if line.startswith('SOL'):
        final_top.append('UNL                 1\n')
    final_top.append(line)

with open('topol.top', 'w') as f:
    f.writelines(final_top)
print('  topol.top written')

# Step 6: MDP files (reuse existing or create)
print('\n=== Step 6: Verify MDP files ===')
mdp_dir = os.path.join(WORKDIR, 'mdp')
for name in ['em.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']:
    path = os.path.join(mdp_dir, name)
    exists = os.path.exists(path)
    print(f'  {name}: {"OK" if exists else "MISSING"}')

print('\n=== SETUP COMPLETE ===')
print(f'  System: {n_total} atoms ({n_ternary} ternary + {n_lig} ligand - {removed_waters*3} removed water atoms)')
print(f'  Ready to run: bash run_ligand_md.sh')
