#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
# Validate 20ns full ternary complex trajectory
# Checks: trajectory integrity, energy stability, RMSD, RMSF


# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
set -e
cd $DATA_ROOT/topology
OUTDIR=$DATA_ROOT/results/00-summary
mkdir -p $OUTDIR

echo "=== TRAJECTORY VALIDATION: 20ns Full Ternary Complex ==="
echo "Date: $(date)"
echo ""

# 1. Trajectory integrity
echo "--- Step 1: Trajectory integrity ---"
gmx check -f md.xtc 2>&1 | tail -10 > $OUTDIR/traj_check.txt
cat $OUTDIR/traj_check.txt
echo ""

# 2. Energy stability (temperature, pressure, potential)
echo "--- Step 2: Energy stability ---"
echo "11 12 13 0" | gmx energy -f md.edr -o $OUTDIR/energy.xvg 2>&1 | tail -20 > $OUTDIR/energy_summary.txt
cat $OUTDIR/energy_summary.txt
echo ""

# 3. Temperature check
echo "--- Step 3: Temperature ---"
echo "16 0" | gmx energy -f md.edr -o $OUTDIR/temperature.xvg 2>&1 | grep -E "Temperature|Average|RMSD"
echo ""

# 4. RMSD of backbone vs initial
echo "--- Step 4: Backbone RMSD ---"
echo "4 4" | gmx rms -s md.tpr -f md.xtc -o $OUTDIR/rmsd_backbone.xvg -tu ns 2>&1 | tail -5
echo ""

# 5. RMSD per chain (protein backbone)
# Need to figure out index groups
echo "--- Step 5: Index groups available ---"
echo "q" | gmx make_ndx -f md.tpr 2>&1 | grep -E "^\s+[0-9]+" | head -20
echo ""

# 6. RMSF per residue
echo "--- Step 6: Per-residue RMSF ---"
echo "4" | gmx rmsf -s md.tpr -f md.xtc -o $OUTDIR/rmsf_backbone.xvg -res 2>&1 | tail -5
echo ""

# 7. Radius of gyration
echo "--- Step 7: Radius of gyration ---"
echo "1" | gmx gyrate -s md.tpr -f md.xtc -o $OUTDIR/gyrate.xvg 2>&1 | tail -5
echo ""

# 8. Number of hydrogen bonds
echo "--- Step 8: Hydrogen bonds (protein-protein) ---"
echo "1 1" | gmx hbond -s md.tpr -f md.xtc -num $OUTDIR/hbonds.xvg 2>&1 | tail -5
echo ""

# 9. Extract final frame as PDB for inspection
echo "--- Step 9: Extract final frame ---"
echo "1" | gmx trjconv -s md.tpr -f md.xtc -o $OUTDIR/final_frame.pdb -dump 20000 -pbc mol 2>&1 | tail -3
echo ""

# 10. Box dimensions over time
echo "--- Step 10: Box dimensions ---"
echo "24 0" | gmx energy -f md.edr -o $OUTDIR/box.xvg 2>&1 | grep -E "Box|Average" | head -5
echo ""

echo "=== VALIDATION COMPLETE ==="
echo "Output files in: $OUTDIR"
ls -la $OUTDIR/
