#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
# Ensemble docking on SOSIP 500ns trajectory
# Usage: bash run_ensemble_dock.sh <PDB_ID> <SOSIP_DIR> <LIGAND_PDBQT>
set -e

PDB_ID=${1:?Usage: $0 PDB_ID SOSIP_DIR LIGAND}
SOSIP_DIR=${2:?}
LIGAND=${3:?}
WORK="${SOSIP_DIR}/mechanistic"
cd "$WORK"

TPR="${SOSIP_DIR}/md_500ns.tpr"
XTC="${SOSIP_DIR}/md.xtc"

echo "=== Ensemble Docking: ${PDB_ID} ==="

# Cluster
echo "GROMOS clustering (2.5 A cutoff, dt=500ps)..."
echo 'Relay_Pocket_CA Relay_Pocket' | gmx cluster \
    -s "$TPR" -f "$XTC" -n pocket.ndx \
    -method gromos -cutoff 0.25 \
    -cl clusters.pdb -sz cluster_sizes.xvg \
    -dt 500 2>&1 | tail -3

N_CLUSTERS=$(grep -c "^MODEL" clusters.pdb 2>/dev/null || echo "0")
echo "Clusters found: ${N_CLUSTERS}"

# Split into individual PDBs
mkdir -p ensemble_frames ensemble_docking
python3 -c "
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

# Dock first 20 clusters max
N_DOCK=$((N_CLUSTERS < 20 ? N_CLUSTERS : 20))
echo "Docking BMS-626529 into ${N_DOCK} frames..."

for i in $(seq 1 $N_DOCK); do
    FRAME=$(printf "ensemble_frames/cluster_%02d.pdb" $i)
    RECEPTOR=$(printf "ensemble_docking/receptor_%02d.pdbqt" $i)
    OUTPUT=$(printf "ensemble_docking/dock_%02d.pdbqt" $i)
    LOG=$(printf "ensemble_docking/dock_%02d.log" $i)

    if [ ! -f "$FRAME" ]; then continue; fi

    # Convert receptor
    obabel "$FRAME" -O "$RECEPTOR" -xr 2>/dev/null

    # Get pocket center
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
")

    CX=$(echo $CENTER | awk '{print $1}')
    CY=$(echo $CENTER | awk '{print $2}')
    CZ=$(echo $CENTER | awk '{print $3}')

    vina --receptor "$RECEPTOR" --ligand "$LIGAND" \
        --center_x $CX --center_y $CY --center_z $CZ \
        --size_x 22 --size_y 22 --size_z 22 \
        --exhaustiveness 16 --num_modes 5 \
        --out "$OUTPUT" > "$LOG" 2>&1

    SCORE=$(grep "^   1" "$LOG" | head -1 | awk '{print $2}')
    printf "  Cluster %02d: %s kcal/mol\n" $i "$SCORE"
done

# Collect results into JSON
python3 -c "
import os, json, re
results = []
for i in range(1, ${N_DOCK}+1):
    log = f'ensemble_docking/dock_{i:02d}.log'
    if not os.path.exists(log): continue
    with open(log) as f:
        for line in f:
            if line.strip().startswith('1'):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        aff = float(parts[1])
                        results.append({'cluster': i, 'affinity_kcal_mol': aff})
                        break
                    except: pass
if results:
    affs = [r['affinity_kcal_mol'] for r in results]
    mean_aff = sum(affs)/len(affs)
    std_aff = (sum((a-mean_aff)**2 for a in affs)/len(affs))**0.5
    report = {
        'pdb_id': '${PDB_ID}',
        'method': 'AutoDock Vina 1.2.7, BMS-626529, GROMOS 2.5A clusters from 500ns MD',
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
    print(f'RESULT: {len(results)} frames, mean={mean_aff:.3f} +/- {std_aff:.3f} kcal/mol')
    print(f'  Best: {min(affs):.3f}, Worst: {max(affs):.3f}')
else:
    print('No docking results')
"

echo "=== ${PDB_ID} Ensemble Docking COMPLETE ==="
