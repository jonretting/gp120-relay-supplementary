#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Setup Steered MD v3 -- Fixed NPT equilibration.

v2 FAILURE DIAGNOSIS:
  - Parrinello-Rahman barostat during NPT equilibration was unstable
  - Average pressure = 84 bar (reference = 1 bar), total drift = -473 bar
  - Box expanded from 21x21x28 nm to 990x990x3675 nm (catastrophic)
  - Root cause: Parrinello-Rahman is a production barostat, NOT for initial
    relaxation. It oscillates wildly when starting far from equilibrium.

v3 FIXES:
  1. NPT uses C-rescale barostat (stochastic cell rescaling) -- stable for
     equilibration, proper ensemble, replaces deprecated Berendsen
  2. refcoord_scaling = com -- prevents position restraints from fighting
     the barostat (crucial for correct NPT with POSRES)
  3. isotropic coupling for NPT (semiisotropic only for production pull)
  4. Production pull.mdp keeps Parrinello-Rahman semiisotropic (correct
     for sampling, stable once NPT has equilibrated the box)

Usage:
  python3 setup_smd_v3.py
  bash run_smd_v3.sh
"""
import subprocess
import os
import math

DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)

WORK = os.path.join(DATA_ROOT, "results", "05-smd-v3")
SRC_GRO = os.path.join(DATA_ROOT, "results", "05-smd-v2", "em.gro")
SRC_TOP = os.path.join(DATA_ROOT, "results", "05-smd-v2", "topol.top")
SRC_NDX = os.path.join(DATA_ROOT, "results", "05-smd-v2", "index.ndx")

os.makedirs(WORK, exist_ok=True)
os.chdir(WORK)


def run(cmd, input_text=None):
    print(f">>> {cmd}")
    r = subprocess.run(cmd, shell=True, input=input_text,
                       capture_output=True, text=True)
    if r.stdout:
        print(r.stdout[-500:])
    if r.stderr:
        print(r.stderr[-500:])
    if r.returncode != 0:
        print(f"WARNING: exit code {r.returncode}")
    return r


print("=" * 60)
print("SMD v3 -- Fixed NPT (C-rescale barostat)")
print("=" * 60)

# Step 1: Copy topology and index from v2 (EM was fine)
print("\nStep 1: Linking topology and EM structure from v2")
for f in ["topol.top", "index.ndx",
          "topol_Protein_chain_A.itp", "topol_Protein_chain_B.itp",
          "topol_Protein_chain_C.itp", "posre_Protein_chain_A.itp",
          "posre_Protein_chain_B.itp", "posre_Protein_chain_C.itp"]:
    src = os.path.join(DATA_ROOT, "results", "05-smd-v2", f)
    if os.path.exists(src):
        run(f"cp {src} {WORK}/")

# Verify EM structure
with open(SRC_GRO) as fh:
    lines = fh.readlines()
n_atoms = int(lines[1].strip())
box = lines[-1].split()
print(f"EM structure: {n_atoms} atoms, box = {' '.join(box[:3])} nm")

if float(box[0]) > 100:
    print("ERROR: EM structure already has exploded box! Cannot proceed.")
    exit(1)

# Step 2: Write MDP files

# NVT (same as v2 -- this worked fine)
with open("nvt.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 50000      ; 100 ps
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500
continuation = no
gen-vel     = yes
gen-temp    = 300
gen-seed    = -1
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
pcoupl      = no
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
define      = -DPOSRES
""")

# NPT -- FIXED: C-rescale + isotropic + refcoord_scaling
with open("npt.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 250000     ; 500 ps
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500
continuation = yes
gen-vel     = no
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
; CRITICAL FIX: C-rescale for equilibration (not Parrinello-Rahman)
pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 5.0
ref-p       = 1.0
compressibility = 4.5e-5
; CRITICAL FIX: scale reference coordinates with COM of pressure group
refcoord-scaling = com
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
define      = -DPOSRES
""")

# SMD pull MDP -- Parrinello-Rahman semiisotropic (fine for production)
with open("pull.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 2500000    ; 5 ns
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500
continuation = yes
gen-vel     = no
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = semiisotropic
tau-p       = 5.0
ref-p       = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
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
pull-coord1-geometry = direction
pull-coord1-vec      = 0 0 1
pull-coord1-dim      = N N Y
pull-coord1-groups   = 1 2
pull-coord1-start    = yes
pull-coord1-rate     = 0.001    ; 0.001 nm/ps = 1 nm/ns
pull-coord1-k        = 1000

pull-nstxout         = 50
pull-nstfout         = 50
""")

# Umbrella window template
with open("umbrella.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 5000000    ; 10 ns per window
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500
continuation = yes
gen-vel     = no
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = semiisotropic
tau-p       = 5.0
ref-p       = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS

pull                 = yes
pull-ncoords         = 1
pull-ngroups         = 2
pull-group1-name     = gp120
pull-group2-name     = CD4
pull-coord1-type     = umbrella
pull-coord1-geometry = direction
pull-coord1-vec      = 0 0 1
pull-coord1-dim      = N N Y
pull-coord1-groups   = 1 2
pull-coord1-start    = yes
pull-coord1-rate     = 0
pull-coord1-k        = 1000

pull-nstxout         = 50
pull-nstfout         = 50
""")

print("MDP files written: nvt.mdp, npt.mdp, pull.mdp, umbrella.mdp")

# Step 3: Write launch script
with open("run_smd_v3.sh", "w") as f:
    f.write(f"""\
#!/bin/bash
set -e

NTOMP=16
cd "$(dirname "$0")"

echo "=== NVT Equilibration (100 ps) ==="
gmx grompp -f nvt.mdp -c {SRC_GRO} -r {SRC_GRO} -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2
gmx mdrun -deffnm nvt -ntomp $NTOMP -pin on
echo "NVT done: $(tail -1 nvt.gro)"

echo "=== NPT Equilibration (500 ps, C-rescale barostat) ==="
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
gmx mdrun -deffnm npt -ntomp $NTOMP -pin on

# Verify box didn't explode
BOX=$(tail -1 npt.gro | awk '{{print $1}}')
echo "NPT final box X: $BOX nm"
if (( $(echo "$BOX > 50" | bc -l) )); then
    echo "ERROR: Box explosion detected (X > 50 nm)! NPT is unstable."
    echo "Check npt.log and npt.edr for pressure/volume issues."
    exit 1
fi
echo "NPT box OK: $(tail -1 npt.gro)"

echo "=== Checking NPT pressure ==="
echo "Pressure" | gmx energy -f npt.edr 2>&1 | tail -3

echo "=== Steered MD (5 ns, 0.001 nm/ps along Z) ==="
gmx grompp -f pull.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -n index.ndx -o pull.tpr -maxwarn 2
nohup gmx mdrun -s pull.tpr -deffnm pull -ntomp $NTOMP -pin on > pull_run.log 2>&1 &
echo "Pull PID: $!"
echo "Monitoring: tail -f pull_run.log"
""")

os.chmod("run_smd_v3.sh", 0o755)

print()
print("=" * 60)
print("v3 SETUP COMPLETE")
print(f"  EM source: {SRC_GRO} ({n_atoms} atoms, box {' '.join(box[:3])} nm)")
print("  Key fix: NPT uses C-rescale barostat + refcoord_scaling = com")
print("  Launch: bash run_smd_v3.sh")
print("=" * 60)
