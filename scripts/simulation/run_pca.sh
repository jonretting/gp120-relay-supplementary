#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
set -e

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

WD="$DATA_ROOT/topology"
cd "$WD"

echo "=== PCA / Essential Dynamics Analysis ==="
echo "Started: $(date)"
echo ""

echo "--- Step 1: Covariance matrix (Relay Pocket) ---"
printf '17\n17\n' | gmx covar -s "$WD/md.tpr" -f "$WD/md.xtc" -n "$WD/analysis.ndx" -o "$WD/eigenvalues.xvg" -v "$WD/eigenvectors.trr" -l "$WD/covar.log" 2>&1 | tail -20

echo ""
echo "--- Step 2: PC projections ---"
printf '17\n17\n' | gmx anaeig -s "$WD/md.tpr" -f "$WD/md.xtc" -n "$WD/analysis.ndx" -v "$WD/eigenvectors.trr" -eig "$WD/eigenvalues.xvg" -proj "$WD/proj_pc1.xvg" -first 1 -last 1 2>&1 | tail -5

printf '17\n17\n' | gmx anaeig -s "$WD/md.tpr" -f "$WD/md.xtc" -n "$WD/analysis.ndx" -v "$WD/eigenvectors.trr" -eig "$WD/eigenvalues.xvg" -proj "$WD/proj_pc2.xvg" -first 2 -last 2 2>&1 | tail -5

printf '17\n17\n' | gmx anaeig -s "$WD/md.tpr" -f "$WD/md.xtc" -n "$WD/analysis.ndx" -v "$WD/eigenvectors.trr" -eig "$WD/eigenvalues.xvg" -proj "$WD/proj_pc3.xvg" -first 3 -last 3 2>&1 | tail -5

echo ""
echo "--- Step 3: 2D projection PC1 vs PC2 ---"
printf '17\n17\n' | gmx anaeig -s "$WD/md.tpr" -f "$WD/md.xtc" -n "$WD/analysis.ndx" -v "$WD/eigenvectors.trr" -eig "$WD/eigenvalues.xvg" -2d "$WD/2d_pc1_pc2.xvg" -first 1 -last 2 2>&1 | tail -5

echo ""
echo "--- Step 4: Cosine content ---"
for pc in 1 2 3; do
    printf '17\n17\n' | gmx anaeig -s "$WD/md.tpr" -f "$WD/md.xtc" -n "$WD/analysis.ndx" -v "$WD/eigenvectors.trr" -eig "$WD/eigenvalues.xvg" -first $pc -last $pc 2>&1 | grep -i "cosine" || true
done

echo ""
echo "--- Step 5: Variance explained ---"
python3 -c "
vals = []
for line in open('$WD/eigenvalues.xvg'):
    if line.startswith('#') or line.startswith('@'):
        continue
    parts = line.split()
    if len(parts) >= 2:
        vals.append(float(parts[1]))

total = sum(vals)
cumulative = 0
print(f'Total variance: {total:.4f} nm^2')
print(f'{\"PC\":>4} {\"Eigenvalue\":>12} {\"% Variance\":>12} {\"Cumulative %\":>14}')
print('-' * 46)
for i, v in enumerate(vals[:10]):
    cumulative += v
    pct = v / total * 100
    cum_pct = cumulative / total * 100
    print(f'{i+1:>4} {v:>12.4f} {pct:>11.1f}% {cum_pct:>13.1f}%')
print(f'')
print(f'Top 3 PCs explain {sum(vals[:3])/total*100:.1f}% of variance')
print(f'Top 5 PCs explain {sum(vals[:5])/total*100:.1f}% of variance')
"

echo ""
echo "=== PCA COMPLETE ==="
echo "Finished: $(date)"
