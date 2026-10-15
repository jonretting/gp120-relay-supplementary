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
LOG=$WORK/pipeline.log

exec > "$LOG" 2>&1

echo "========================================"
echo "FULL TERNARY COMPLEX MD PIPELINE"
echo "Machine: $(hostname), $(nproc) cores"
echo "Started: $(date)"
echo "========================================"

cd "$MD"

echo ""
echo "--- Step 1: Energy Minimization ---"
echo "$(date)"
gmx grompp -f "$MDP/em.mdp" -c ionized.gro -p topol.top -o em.tpr -maxwarn 10
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 20
echo "EM complete: $(date)"

echo ""
echo "--- Step 2: NVT Equilibration (100 ps) ---"
echo "$(date)"
gmx grompp -f "$MDP/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 10
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 20
echo "NVT complete: $(date)"

echo ""
echo "--- Step 3: NPT Equilibration (100 ps) ---"
echo "$(date)"
gmx grompp -f "$MDP/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 10
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 20
echo "NPT complete: $(date)"

echo ""
echo "--- Step 4: Production MD (20 ns) ---"
echo "$(date)"
gmx grompp -f "$MDP/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 10
gmx mdrun -deffnm md -ntmpi 1 -ntomp 20
echo "Production complete: $(date)"

echo ""
echo "========================================"
echo "PIPELINE COMPLETE: $(date)"
echo "========================================"
