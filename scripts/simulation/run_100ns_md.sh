#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
set -e

WORK=$DATA_ROOT
MD=$WORK/md
MDP=$WORK/mdp

cd $MD

echo "============================================================"
echo "100ns MD SIMULATION -- gp120 relay pocket"
echo "Machine: $(hostname), $(nproc) cores"
echo "============================================================"

# Update md.mdp for 100ns
cat > $MDP/md.mdp << 'MDPEOF'
; Production MD -- 100 ns
integrator  = md
nsteps      = 50000000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = yes
gen_vel     = no
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
MDPEOF

echo ""
echo "--- Step 1: pdb2gmx ---"
echo '6' | gmx pdb2gmx -f reduced_system.pdb -o processed.gro -water tip3p -ignh 2>&1 | tail -5

echo ""
echo "--- Step 2: editconf ---"
gmx editconf -f processed.gro -o boxed.gro -c -d 1.5 -bt cubic 2>&1 | tail -5

echo ""
echo "--- Step 3: solvate ---"
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | tail -3

echo ""
echo "--- Step 4: add ions ---"
gmx grompp -f $MDP/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 10 2>&1 | tail -3
echo 'SOL' | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 2>&1 | tail -3

echo ""
echo "--- Step 5: energy minimization ---"
gmx grompp -f $MDP/em.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 10 2>&1 | tail -3
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 16 2>&1 | tail -5

echo ""
echo "--- Step 6: NVT equilibration (100 ps) ---"
gmx grompp -f $MDP/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 10 2>&1 | tail -3
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 16 2>&1 | tail -5

echo ""
echo "--- Step 7: NPT equilibration (100 ps) ---"
gmx grompp -f $MDP/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 10 2>&1 | tail -3
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 16 2>&1 | tail -3

echo ""
echo "--- Step 8: Production MD (100 ns) ---"
gmx grompp -f $MDP/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 10 2>&1 | tail -3

echo "Starting 100ns production run..."
echo "Estimated time: ~3-5 days on 16 CPU cores"
nohup gmx mdrun -deffnm md -ntmpi 1 -ntomp 16 > md_production.log 2>&1 &
MD_PID=$!
echo "Production MD started with PID: $MD_PID"
echo $MD_PID > md_pid.txt

echo ""
echo "============================================================"
echo "SETUP COMPLETE"
echo "  Monitor: tail -f $MD/md_production.log"
echo "  Progress: gmx check -f $MD/md.xtc"
echo "  PID: $MD_PID"
echo "============================================================"
