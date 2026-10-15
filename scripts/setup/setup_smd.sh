#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
set -euo pipefail

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

WORKDIR="$DATA_ROOT/steered-md"
MDDIR="$DATA_ROOT/topology"
TOPOL="$MDDIR/topol.top"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ============================================================
# 1. Create chain-specific index groups
# ============================================================
# Chain A (gp120): GROMACS resids 1-398 (atoms 1-6204)
# Chain B (CD4):   GROMACS resids 1-176 (atoms 6205-8977)
# Chain C (CCR5):  GROMACS resids 1-313 (atoms 8978-14149)

cat > make_chains.input <<'EOF'
r 1-398 & 1
name 17 gp120
r 1-176 & 1
name 18 CD4_full
q
EOF

echo "Making chain index groups..."
# We need to split by atom ranges since resids overlap across chains
# Use the structure file to extract chain atom ranges

python3 << 'PYEOF'
import re

gro_file = "$MDDIR/md.gro"
with open(gro_file) as f:
    lines = f.readlines()

n_atoms = int(lines[1].strip())
protein_atoms = []
chain_a = []  # gp120
chain_b = []  # CD4
chain_c = []  # CCR5

current_chain = 'A'
prev_resid = 0
chain_switch = 0

for i, line in enumerate(lines[2:2+n_atoms], start=1):
    if len(line) < 44:
        continue
    resid_str = line[:5].strip()
    if not resid_str.isdigit():
        continue
    resid = int(resid_str)
    atom_name = line[10:15].strip()
    
    # Detect chain switches by resid reset
    if resid < prev_resid and prev_resid > 100:
        chain_switch += 1
        if chain_switch == 1:
            current_chain = 'B'
        elif chain_switch == 2:
            current_chain = 'C'
    prev_resid = resid
    
    if current_chain == 'A':
        chain_a.append(i)
    elif current_chain == 'B':
        chain_b.append(i)
    elif current_chain == 'C':
        chain_c.append(i)

def write_ndx_group(f, name, atoms):
    f.write(f"[ {name} ]\n")
    for j, a in enumerate(atoms, 1):
        f.write(f"{a:7d}")
        if j % 15 == 0:
            f.write("\n")
    if len(atoms) % 15 != 0:
        f.write("\n")

with open("chains.ndx", "w") as f:
    write_ndx_group(f, "gp120", chain_a)
    write_ndx_group(f, "CD4", chain_b)
    write_ndx_group(f, "CCR5", chain_c)

# Relay pocket: chain A resids 372-379
relay = [a for a, line in zip(chain_a, 
         [lines[2+a-1] for a in chain_a]) 
         if int(line[:5].strip()) in range(372, 380)]

write_ndx_group(open("chains.ndx", "a"), "Relay_Pocket", relay)

print(f"Chain A (gp120): {len(chain_a)} atoms ({chain_a[0]}-{chain_a[-1]})")
print(f"Chain B (CD4):   {len(chain_b)} atoms ({chain_b[0]}-{chain_b[-1]})")
print(f"Chain C (CCR5):  {len(chain_c)} atoms ({chain_c[0]}-{chain_c[-1]})")
print(f"Relay pocket:    {len(relay)} atoms")
PYEOF

# Merge with existing index
cat "$MDDIR/index.ndx" chains.ndx > combined.ndx

echo "Chain index groups created."

# ============================================================
# 2. Extract the last frame as starting structure
# ============================================================
echo "Extracting last frame from 20ns trajectory..."
echo "0" | gmx trjconv -s "$MDDIR/md.tpr" -f "$MDDIR/md.xtc" \
    -o start_frame.gro -dump 20000 -n "$MDDIR/index.ndx" 2>/dev/null

echo "Starting structure extracted."

# ============================================================
# 3. Create the steered MD (pulling) MDP file
# ============================================================
cat > pull.mdp << 'MDP'
; Steered MD -- pull CD4 from gp120
integrator  = md
dt          = 0.002      ; 2 fs
nsteps      = 500000     ; 1 ns pulling
nstxout     = 0
nstxout-compressed = 500 ; every 1 ps
nstenergy   = 500
nstlog      = 500

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300

; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5

; Electrostatics
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pme-order   = 4
fourierspacing = 0.12

; Constraints
constraints = h-bonds
constraint-algorithm = LINCS
lincs-iter  = 1
lincs-order = 4

; PBC
pbc         = xyz

; COM motion removal
nstcomm     = 100
comm-mode   = linear
comm-grps   = System

; Pull code
pull                 = yes
pull-ncoords         = 1
pull-ngroups         = 2
pull-group1-name     = gp120
pull-group2-name     = CD4
pull-coord1-type     = umbrella
pull-coord1-geometry = distance
pull-coord1-dim      = Y Y Y
pull-coord1-groups   = 1 2
pull-coord1-start    = yes
pull-coord1-rate     = 0.01     ; 0.01 nm/ps = 10 m/s
pull-coord1-k        = 1000     ; kJ/mol/nm^2
MDP

echo "Pull MDP created."

# ============================================================
# 4. Create the umbrella sampling MDP file
# ============================================================
cat > umbrella.mdp << 'MDP'
; Umbrella sampling window
integrator  = md
dt          = 0.002
nsteps      = 5000000   ; 10 ns per window
nstxout     = 0
nstxout-compressed = 5000  ; every 10 ps
nstenergy   = 5000
nstlog      = 5000

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300

; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5

; Electrostatics
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pme-order   = 4
fourierspacing = 0.12

; Constraints
constraints = h-bonds
constraint-algorithm = LINCS
lincs-iter  = 1
lincs-order = 4

; PBC
pbc         = xyz

; COM motion removal
nstcomm     = 100
comm-mode   = linear
comm-grps   = System

; Pull code -- umbrella restraint (no rate)
pull                 = yes
pull-ncoords         = 1
pull-ngroups         = 2
pull-group1-name     = gp120
pull-group2-name     = CD4
pull-coord1-type     = umbrella
pull-coord1-geometry = distance
pull-coord1-dim      = Y Y Y
pull-coord1-groups   = 1 2
pull-coord1-start    = yes
pull-coord1-rate     = 0.0      ; no pulling, just restraint
pull-coord1-k        = 1000     ; kJ/mol/nm^2
MDP

echo "Umbrella MDP created."

# ============================================================
# 5. Generate TPR for steered MD
# ============================================================
echo "Generating steered MD TPR..."
gmx grompp -f pull.mdp \
    -c start_frame.gro \
    -p "$TOPOL" \
    -n combined.ndx \
    -o pull.tpr \
    -maxwarn 2 2>&1 | tail -10

echo ""
echo "================================================"
echo "Setup complete. Files in $WORKDIR:"
echo "  pull.tpr     - steered MD run input"
echo "  pull.mdp     - steered MD parameters"
echo "  umbrella.mdp - umbrella sampling template"
echo "  combined.ndx - index with chain groups"
echo "  start_frame.gro - starting structure (20ns)"
echo ""
echo "To run steered MD:"
echo "  cd $WORKDIR && gmx mdrun -s pull.tpr -deffnm pull -ntomp 12"
echo ""
echo "After pulling, extract umbrella windows with:"
echo "  python3 extract_umbrella_windows.py"
echo "================================================"
