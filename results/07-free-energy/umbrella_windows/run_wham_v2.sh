#!/bin/sh
echo "=== WHAM v2 — Full 29-Window PMF ==="
echo "Date: $(date -u)"
echo "Windows: 29 (w00-w28, including clean w24 re-run)"
echo "Range: -1.5 to 4.4 nm, 120 bins, tol=1e-4, skip first 500 ps"
echo ""
echo "--- Step 1: Standard WHAM (200 bootstraps) ---"
gmx wham -it /umbrella/tpr_files_v2.dat -ix /umbrella/pullx_files_v2.dat -o /umbrella/pmf_v2.xvg -hist /umbrella/histo_v2.xvg -bsres /umbrella/pmf_v2_bootstrap.xvg -bsprof /umbrella/pmf_v2_bsprofiles.xvg -unit kCal -temp 310 -nBootstrap 200 -bs-method hist -tol 1e-4 -min -1.5 -max 4.4 -bins 120 -b 500 2>&1 | tee /umbrella/wham_v2.log
echo ""
echo "--- Step 2: Convergence check (tau=10) ---"
gmx wham -it /umbrella/tpr_files_v2.dat -ix /umbrella/pullx_files_v2.dat -o /umbrella/pmf_v2_ac10.xvg -unit kCal -temp 310 -nBootstrap 200 -bs-method hist -bs-tau 10 -tol 1e-4 -min -1.5 -max 4.4 -bins 120 -b 500 2>&1 | tee /umbrella/wham_v2_ac10.log
echo ""
echo "--- Step 3: Time-convergence first half (500-5000 ps) ---"
gmx wham -it /umbrella/tpr_files_v2.dat -ix /umbrella/pullx_files_v2.dat -o /umbrella/pmf_v2_first.xvg -unit kCal -temp 310 -nBootstrap 100 -bs-method hist -tol 1e-4 -min -1.5 -max 4.4 -bins 120 -b 500 -e 5000 2>&1 | tee /umbrella/wham_v2_first.log
echo ""
echo "--- Step 4: Time-convergence second half (5000-10000 ps) ---"
gmx wham -it /umbrella/tpr_files_v2.dat -ix /umbrella/pullx_files_v2.dat -o /umbrella/pmf_v2_second.xvg -unit kCal -temp 310 -nBootstrap 100 -bs-method hist -tol 1e-4 -min -1.5 -max 4.4 -bins 120 -b 5000 -e 10000 2>&1 | tee /umbrella/wham_v2_second.log
echo ""
echo "=== WHAM v2 Complete ==="
ls -lh /umbrella/pmf_v2*.xvg /umbrella/histo_v2.xvg /umbrella/wham_v2*.log 2>/dev/null
echo ""
echo "=== PMF Summary (pmf_v2.xvg) ==="
grep -v "^[#@]" /umbrella/pmf_v2.xvg
echo ""
echo "=== Bootstrap Error Summary ==="
grep -v "^[#@]" /umbrella/pmf_v2_bootstrap.xvg
