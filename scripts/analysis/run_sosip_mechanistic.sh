#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
# Run SOSIP mechanistic cross-validation (DCCM + PCA + ensemble docking)
# on a completed 500ns SOSIP trajectory.
#
# Usage: bash run_sosip_mechanistic.sh <PDB_ID> <SOSIP_DIR>
# Example: bash run_sosip_mechanistic.sh 5U7O $DATA_ROOT/results/04-md-simulations/5u7o
#
# Prerequisites: GROMACS 2025.2, Python3 with numpy/networkx, AutoDock Vina 1.2.7
# Input files: md.xtc, md_500ns.tpr (or md.tpr), md.gro, pocket.ndx (from analyze_sosip_500ns.py)


# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
set -e

PDB_ID=${1:?Usage: $0 <PDB_ID> <SOSIP_DIR>}
SOSIP_DIR=${2:?Usage: $0 <PDB_ID> <SOSIP_DIR>}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORK_DIR="${SOSIP_DIR}/mechanistic"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

TPR="${SOSIP_DIR}/md_500ns.tpr"
if [ ! -f "$TPR" ]; then
    TPR="${SOSIP_DIR}/md.tpr"
fi
XTC="${SOSIP_DIR}/md.xtc"
GRO="${SOSIP_DIR}/md.gro"

echo "================================================================"
echo "  SOSIP Mechanistic Cross-Validation: ${PDB_ID}"
echo "  Working directory: ${WORK_DIR}"
echo "  TPR: ${TPR}"
echo "  XTC: ${XTC}"
echo "================================================================"

# Ensure pocket.ndx exists (from prior analyze_sosip_500ns.py run)
if [ ! -f "${SOSIP_DIR}/pocket.ndx" ] && [ ! -f "pocket.ndx" ]; then
    echo "ERROR: pocket.ndx not found. Run analyze_sosip_500ns.py first."
    exit 1
fi
cp -n "${SOSIP_DIR}/pocket.ndx" pocket.ndx 2>/dev/null || true

# ============================================================
# PHASE 1: DCCM (Dynamic Cross-Correlation Matrix)
# ============================================================
echo ""
echo "=== PHASE 1: DCCM ==="

# Run gmx covar to get covariance matrix
echo "Running gmx covar (dt=50ps, All_CA fit and analysis)..."
echo 'All_CA All_CA' | gmx covar \
    -s "$TPR" -f "$XTC" -n pocket.ndx \
    -ascii covar_ascii.dat \
    -o eigenvalues_covar.xvg \
    -v eigenvectors_covar.trr \
    -dt 50 2>&1 | tail -5

# Compute DCCM from covariance
echo "Computing DCCM..."
python3 "${SCRIPT_DIR}/compute_sosip_dccm_500ns.py" "$PDB_ID" covar_ascii.dat "$GRO"

echo "DCCM complete: ${PDB_ID}_500ns_dccm.json"

# ============================================================
# PHASE 2: PCA / Essential Dynamics
# ============================================================
echo ""
echo "=== PHASE 2: PCA ==="

# PCA on relay pocket atoms only
echo "Running gmx covar on relay pocket..."
echo 'Relay_Pocket Relay_Pocket' | gmx covar \
    -s "$TPR" -f "$XTC" -n pocket.ndx \
    -o eigenvalues_pocket.xvg \
    -v eigenvectors_pocket.trr \
    -dt 50 2>&1 | tail -5

# Extract eigenvalues
echo "Running gmx anaeig..."
echo 'Relay_Pocket' | gmx anaeig \
    -s "$TPR" -f "$XTC" -n pocket.ndx \
    -v eigenvectors_pocket.trr \
    -entropy \
    -eig eigenvalues_pocket.xvg 2>&1 | tail -10

# Extract PC projections
for pc in 1 2 3; do
    echo 'Relay_Pocket' | gmx anaeig \
        -s "$TPR" -f "$XTC" -n pocket.ndx \
        -v eigenvectors_pocket.trr \
        -eig eigenvalues_pocket.xvg \
        -proj "proj_pc${pc}.xvg" \
        -first $pc -last $pc 2>&1 > /dev/null
done

# 2D free energy landscape (PC1 vs PC2)
echo 'Relay_Pocket' | gmx anaeig \
    -s "$TPR" -f "$XTC" -n pocket.ndx \
    -v eigenvectors_pocket.trr \
    -eig eigenvalues_pocket.xvg \
    -2d 2d_pc1_pc2.xvg \
    -first 1 -last 2 2>&1 > /dev/null

# Parse eigenvalues into JSON
python3 -c "
import json
vals = []
with open('eigenvalues_pocket.xvg') as f:
    for line in f:
        if line.startswith(('#','@')): continue
        parts = line.split()
        if len(parts) >= 2:
            vals.append(float(parts[1]))
total = sum(vals) if vals else 1.0
pcs = []
cum = 0
for i, v in enumerate(vals[:10]):
    pct = 100.0 * v / total
    cum += pct
    pcs.append({'pc': i+1, 'eigenvalue_nm2': round(v, 4), 'pct_variance': round(pct, 1), 'cumulative_pct': round(cum, 1)})
result = {
    'pdb_id': '${PDB_ID}',
    'method': 'PCA via gmx covar + gmx anaeig on relay pocket, 500ns, dt=50ps',
    'total_variance_nm2': round(total, 4),
    'top_3_pct': round(sum(p['pct_variance'] for p in pcs[:3]), 1),
    'top_5_pct': round(sum(p['pct_variance'] for p in pcs[:5]), 1),
    'principal_components': pcs,
    'comparison_6meo_100ns': {
        'top_3_pct': 51.0,
        'top_5_pct': 65.5,
        'total_variance_nm2': 2.4776,
        'note': '6MEO was 100ns/105 atoms (full reduced system); SOSIP is 500ns/pocket atoms only'
    }
}
with open('${PDB_ID}_500ns_pca.json', 'w') as f:
    json.dump(result, f, indent=2)
print(f'PCA: top 3 PCs = {result[\"top_3_pct\"]}%, total variance = {result[\"total_variance_nm2\"]:.4f} nm^2')
"

echo "PCA complete: ${PDB_ID}_500ns_pca.json"

# ============================================================
# PHASE 3: Ensemble Docking (if vina available)
# ============================================================
echo ""
echo "=== PHASE 3: Ensemble Docking ==="

if command -v vina &> /dev/null; then
    # Cluster the trajectory
    echo "GROMOS clustering (2.5 A cutoff)..."
    echo 'Relay_Pocket_CA Relay_Pocket' | gmx cluster \
        -s "$TPR" -f "$XTC" -n pocket.ndx \
        -method gromos -cutoff 0.25 \
        -cl clusters.pdb -sz cluster_sizes.xvg \
        -dt 500 2>&1 | tail -5

    # Count clusters
    N_CLUSTERS=$(grep -c "^MODEL" clusters.pdb 2>/dev/null || echo "0")
    echo "Found ${N_CLUSTERS} clusters"

    if [ "$N_CLUSTERS" -gt 0 ]; then
        mkdir -p ensemble_frames
        # Split cluster PDB into individual frames
        python3 -c "
import os
frames = []
current = []
model_num = 0
with open('clusters.pdb') as f:
    for line in f:
        if line.startswith('MODEL'):
            model_num += 1
            current = [line]
        elif line.startswith('ENDMDL'):
            current.append(line)
            outf = f'ensemble_frames/cluster_{model_num:02d}.pdb'
            with open(outf, 'w') as o:
                o.writelines(current)
            frames.append(outf)
            current = []
        else:
            current.append(line)
print(f'Split {len(frames)} cluster PDBs')
"
        # Prepare PDBQT and dock (first 20 clusters max)
        N_DOCK=$((N_CLUSTERS < 20 ? N_CLUSTERS : 20))
        echo "Docking BMS-626529 into ${N_DOCK} cluster representatives..."

        # BMS-626529 PDBQT would need to be available
        BMS_PDBQT="${SCRIPT_DIR}/../../structures/bms626529.pdbqt"
        if [ -f "$BMS_PDBQT" ]; then
            mkdir -p ensemble_docking
            for i in $(seq 1 $N_DOCK); do
                FRAME=$(printf "ensemble_frames/cluster_%02d.pdb" $i)
                RECEPTOR=$(printf "ensemble_docking/receptor_%02d.pdbqt" $i)
                OUTPUT=$(printf "ensemble_docking/dock_%02d.pdbqt" $i)

                # Convert to PDBQT
                obabel "$FRAME" -O "$RECEPTOR" -xr 2>/dev/null || continue

                # Get pocket center from the frame
                CENTER=$(python3 -c "
with open('$FRAME') as f:
    xs, ys, zs = [], [], []
    for line in f:
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
if xs:
    print(f'{sum(xs)/len(xs):.1f} {sum(ys)/len(ys):.1f} {sum(zs)/len(zs):.1f}')
else:
    print('0 0 0')
" 2>/dev/null)

                CX=$(echo $CENTER | awk '{print $1}')
                CY=$(echo $CENTER | awk '{print $2}')
                CZ=$(echo $CENTER | awk '{print $3}')

                vina --receptor "$RECEPTOR" --ligand "$BMS_PDBQT" \
                    --center_x $CX --center_y $CY --center_z $CZ \
                    --size_x 22 --size_y 22 --size_z 22 \
                    --exhaustiveness 16 --num_modes 5 \
                    --out "$OUTPUT" 2>&1 | grep "^   1" | head -1
            done

            # Collect results
            python3 -c "
import os, json, re
results = []
for i in range(1, ${N_DOCK}+1):
    outf = f'ensemble_docking/dock_{i:02d}.pdbqt'
    if not os.path.exists(outf): continue
    with open(outf) as f:
        for line in f:
            m = re.search(r'VINA RESULT:\s+([\-\d\.]+)', line)
            if m:
                results.append({'cluster': i, 'affinity_kcal_mol': float(m.group(1))})
                break
if results:
    affs = [r['affinity_kcal_mol'] for r in results]
    mean_aff = sum(affs)/len(affs)
    std_aff = (sum((a-mean_aff)**2 for a in affs)/len(affs))**0.5
    report = {
        'pdb_id': '${PDB_ID}',
        'method': 'AutoDock Vina 1.2.7 ensemble docking, BMS-626529, GROMOS 2.5A clusters from 500ns',
        'n_frames': len(results),
        'mean_affinity_kcal_mol': round(mean_aff, 3),
        'std_affinity_kcal_mol': round(std_aff, 3),
        'best_affinity_kcal_mol': round(min(affs), 3),
        'worst_affinity_kcal_mol': round(max(affs), 3),
        'per_cluster': results,
        'comparison_6meo_100ns': {
            'mean': -6.808, 'std': 0.750, 'n_frames': 20
        }
    }
    with open('${PDB_ID}_500ns_ensemble_docking.json', 'w') as f:
        json.dump(report, f, indent=2)
    print(f'Ensemble docking: {len(results)} frames, mean={mean_aff:.3f} +/- {std_aff:.3f} kcal/mol')
else:
    print('No docking results produced')
"
        else
            echo "BMS-626529 PDBQT not found at ${BMS_PDBQT} -- skipping docking"
            echo "(Clustering done; docking can be run later with ligand file)"
        fi
    fi
else
    echo "AutoDock Vina not found -- skipping ensemble docking"
    echo "(DCCM and PCA results are complete)"
fi

echo ""
echo "================================================================"
echo "  ${PDB_ID} Mechanistic Cross-Validation COMPLETE"
echo "  Results in: ${WORK_DIR}"
echo "================================================================"
