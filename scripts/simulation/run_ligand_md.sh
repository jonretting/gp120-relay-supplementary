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

cd $DATA_ROOT/results/04-md-simulations

echo "=== Step 1: Energy Minimization ==="
gmx grompp -f mdp/em.mdp -c combined.gro -p topol.top -o em.tpr -maxwarn 5 2>&1 | tail -20
gmx mdrun -deffnm em -ntmpi 1 -ntomp 24 2>&1 | tail -10
echo "EM complete: $(grep 'Potential Energy' em.log | tail -1)"

echo "=== Step 2: NVT Equilibration (100 ps) ==="
gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 5 2>&1 | tail -20
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 24 2>&1 | tail -10
echo "NVT complete"

echo "=== Step 3: NPT Equilibration (100 ps) ==="
gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 5 2>&1 | tail -10
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 24 2>&1 | tail -10
echo "NPT complete"

echo "=== Step 4: Production MD (10 ns) ==="
gmx grompp -f mdp/md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 5 2>&1 | tail -10
gmx mdrun -deffnm md -ntmpi 1 -ntomp 24 -pin on 2>&1 | tail -10
echo "=== PRODUCTION MD COMPLETE ==="
echo "Output: md.xtc, md.gro, md.edr, md.log"
