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
cd $DATA_ROOT/results/04-md-simulations/positioned_acpype
source ~/miniforge3/bin/activate

echo '=== Step 1: antechamber - GAFF2 types + BCC charges ==='
antechamber -dr no \
  -i bms626529_positioned.pdb -fi pdb \
  -o bms626529_gaff2.mol2 -fo mol2 \
  -c bcc -nc 0 -m 1 -s 2 -df 2 -at gaff2 -pf n
echo 'antechamber done'

echo '=== Step 2: parmchk2 - generate frcmod ==='
parmchk2 -i bms626529_gaff2.mol2 -f mol2 -o bms626529.frcmod -s 2
echo 'parmchk2 done'

echo '=== Step 3: tleap - build topology ==='
cat > tleap.in << 'TLEAP'
source leaprc.gaff2
UNL = loadmol2 bms626529_gaff2.mol2
loadamberparams bms626529.frcmod
check UNL
saveamberparm UNL bms626529.prmtop bms626529.inpcrd
savepdb UNL bms626529_full.pdb
quit
TLEAP
tleap -f tleap.in
echo 'tleap done'

echo '=== Step 4: parmed - convert to GROMACS ==='
python3 << 'PYEOF'
import parmed
parm = parmed.load_file('bms626529.prmtop', 'bms626529.inpcrd')
parm.save('bms626529_GMX.top', overwrite=True)
parm.save('bms626529_GMX.gro', overwrite=True)
print(f'GROMACS files: {len(parm.atoms)} atoms, {len(parm.residues)} residues')
PYEOF
echo '=== ALL DONE ==='
ls -la bms626529*
