#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
set -e

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

BD="$DATA_ROOT"
WD="$DATA_ROOT/topology"
cd "$WD"

echo "=== Ensemble Docking Analysis ==="
echo "Started: $(date)"
echo ""

# Clean up old results
rm -rf "$DATA_ROOT/results/02-docking/ensemble_frames" "$DATA_ROOT/results/02-docking/ensemble_poses"
rm -f "$WD/clusters.pdb" "$WD/cluster_index.ndx" "$WD/cluster.log" "$WD/cluster_sizes.xvg" "$WD/cluster_rmsd.xvg"

echo "--- Step 1: Clustering backbone (GROMOS method, 2.5 A cutoff for fewer clusters) ---"
printf '4\n1\n' | gmx cluster -s "$WD/md.tpr" -f "$WD/md.xtc" \
  -method gromos -cutoff 0.25 \
  -cl "$WD/clusters.pdb" -clndx "$WD/cluster_index.ndx" \
  -g "$WD/cluster.log" -sz "$WD/cluster_sizes.xvg" \
  -dist "$WD/cluster_rmsd.xvg" 2>&1 | tail -10

NCLUST=$(grep -c 'ENDMDL' "$WD/clusters.pdb" 2>/dev/null || echo "?")
echo "Found $NCLUST cluster representatives"
echo ""

mkdir -p "$DATA_ROOT/results/02-docking/ensemble_frames" "$DATA_ROOT/results/02-docking/ensemble_poses"

echo "--- Step 2: Extracting cluster PDBs (max 20) ---"
python3 -c "
import os
frames_dir = '$DATA_ROOT/results/02-docking/ensemble_frames'
current_model = []
model_num = 0
max_models = 20
with open('$WD/clusters.pdb') as f:
    for line in f:
        if line.startswith('MODEL'):
            model_num += 1
            if model_num > max_models:
                break
            current_model = [line]
        elif line.startswith('ENDMDL'):
            current_model.append(line)
            fname = os.path.join(frames_dir, f'cluster_{model_num:02d}.pdb')
            with open(fname, 'w') as out:
                out.writelines(current_model)
            current_model = []
        elif current_model:
            current_model.append(line)
print(f'Extracted {min(model_num, max_models)} cluster PDBs (of {model_num} total)')
"

echo ""
echo "--- Step 3: Computing pocket center from cluster PDBs ---"
POCKET_CENTER=$(python3 -c "
import os, glob
xs, ys, zs = [], [], []
for pdb in sorted(glob.glob('$DATA_ROOT/results/02-docking/ensemble_frames/cluster_01.pdb'))[:1]:
    for line in open(pdb):
        if not line.startswith('ATOM'):
            continue
        resnum = int(line[22:26].strip())
        # Pocket residues in this system: MET94 starts at index 69 in topology
        # but PDB renumbering from cluster extraction may differ
        # Just use ALL protein atoms for center estimation, then refine
        xs.append(float(line[30:38]))
        ys.append(float(line[38:46]))
        zs.append(float(line[46:54]))
if xs:
    print(f'{sum(xs)/len(xs):.1f},{sum(ys)/len(ys):.1f},{sum(zs)/len(zs):.1f}')
else:
    print('0,0,0')
")

CENTER_X=$(echo $POCKET_CENTER | cut -d, -f1)
CENTER_Y=$(echo $POCKET_CENTER | cut -d, -f2)
CENTER_Z=$(echo $POCKET_CENTER | cut -d, -f3)
echo "Pocket center: ($CENTER_X, $CENTER_Y, $CENTER_Z)"
SIZE=25

echo ""
echo "--- Step 4: Converting cluster PDBs to PDBQT ---"
for pdb in "$DATA_ROOT"/results/02-docking/ensemble_frames/cluster_*.pdb; do
    base=$(basename "$pdb" .pdb)
    obabel "$pdb" -O "$DATA_ROOT/results/02-docking/ensemble_frames/${base}.pdbqt" -xr 2>/dev/null
    echo "  Converted: $base"
done

echo ""
echo "--- Step 5: Docking into each cluster ---"

for pdbqt in "$DATA_ROOT"/results/02-docking/ensemble_frames/cluster_*.pdbqt; do
    base=$(basename "$pdbqt" .pdbqt)
    echo "  Docking $base..."

    vina --receptor "$pdbqt" --ligand "$DATA_ROOT/structures/bms626529.pdbqt" \
         --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \
         --size_x $SIZE --size_y $SIZE --size_z $SIZE \
         --exhaustiveness 16 \
         --out "$DATA_ROOT/results/02-docking/ensemble_poses/${base}_bms626529.pdbqt" 2>&1 || echo "    BMS-626529 docking failed for $base"

    vina --receptor "$pdbqt" --ligand "$DATA_ROOT/structures/vicriviroc.pdbqt" \
         --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \
         --size_x $SIZE --size_y $SIZE --size_z $SIZE \
         --exhaustiveness 16 \
         --out "$DATA_ROOT/results/02-docking/ensemble_poses/${base}_vicriviroc.pdbqt" 2>&1 || echo "    Vicriviroc docking failed for $base"

    bms_score=$(grep 'REMARK VINA RESULT' "$DATA_ROOT/results/02-docking/ensemble_poses/${base}_bms626529.pdbqt" 2>/dev/null | head -1 | awk '{print $4}')
    vic_score=$(grep 'REMARK VINA RESULT' "$DATA_ROOT/results/02-docking/ensemble_poses/${base}_vicriviroc.pdbqt" 2>/dev/null | head -1 | awk '{print $4}')
    echo "    BMS-626529: ${bms_score:-N/A} kcal/mol | Vicriviroc: ${vic_score:-N/A} kcal/mol"
done

echo ""
echo "--- Step 6: Summary ---"
python3 -c "
import os, glob, statistics

results = {'bms626529': [], 'vicriviroc': []}
resdir = '$DATA_ROOT/results/02-docking/ensemble_poses'

for pf in sorted(glob.glob(os.path.join(resdir, '*_bms626529.pdbqt'))):
    with open(pf) as f:
        for line in f:
            if 'REMARK VINA RESULT' in line:
                results['bms626529'].append(float(line.split()[3]))
                break

for pf in sorted(glob.glob(os.path.join(resdir, '*_vicriviroc.pdbqt'))):
    with open(pf) as f:
        for line in f:
            if 'REMARK VINA RESULT' in line:
                results['vicriviroc'].append(float(line.split()[3]))
                break

print('ENSEMBLE DOCKING RESULTS')
print('=' * 60)
for compound, scores in results.items():
    if not scores:
        print(f'{compound}: No results')
        continue
    mean = statistics.mean(scores)
    sd = statistics.stdev(scores) if len(scores) > 1 else 0
    print(f'{compound}:')
    print(f'  Frames docked: {len(scores)}')
    print(f'  Mean affinity:  {mean:.3f} kcal/mol')
    print(f'  SD:             {sd:.3f} kcal/mol')
    print(f'  Min (best):     {min(scores):.3f} kcal/mol')
    print(f'  Max (worst):    {max(scores):.3f} kcal/mol')
    print(f'  Range:          {max(scores)-min(scores):.3f} kcal/mol')
    print(f'  All scores:     {\", \".join(f\"{s:.3f}\" for s in scores)}')
    print()

print('=' * 60)
if results['bms626529'] and results['vicriviroc']:
    bms_mean = statistics.mean(results['bms626529'])
    vic_mean = statistics.mean(results['vicriviroc'])
    print(f'BMS-626529 ensemble mean: {bms_mean:.3f} kcal/mol')
    print(f'Vicriviroc ensemble mean: {vic_mean:.3f} kcal/mol')
    print(f'Original single-frame BMS-626529: -6.692 kcal/mol')
    print(f'Original single-frame Vicriviroc: -8.797 kcal/mol')
"

echo ""
echo "=== ENSEMBLE DOCKING COMPLETE ==="
echo "Finished: $(date)"
