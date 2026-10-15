#!/bin/sh
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
set -e

# SOSIP Reduced-System MD Pipeline
# Usage: bash run_sosip_md.sh <PDB_ID> [NTOMP]
# Expects <PDB_ID>_reduced.pdb in current directory (from setup_sosip_md.py)

PDB_ID="$1"
NTOMP="${2:-16}"
if [ -z "$PDB_ID" ]; then
    echo "Usage: sh run_sosip_md.sh <PDB_ID> [NTOMP]"
    exit 1
fi
REDUCED="${PDB_ID}_reduced.pdb"

if [ ! -f "$REDUCED" ]; then
    echo "ERROR: $REDUCED not found. Run setup_sosip_md.py first."
    exit 1
fi

echo "=== SOSIP MD Pipeline: $PDB_ID ==="
echo "  Input: $REDUCED"
echo "  Threads: $NTOMP"
echo "  Start: $(date)"

# Write MDP files
cat > em.mdp << 'MDP'
integrator  = steep
emtol       = 500.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
MDP

cat > nvt.mdp << 'MDP'
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = no
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 20
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = no
pbc         = xyz
DispCorr    = EnerPres
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
MDP

cat > npt.mdp << 'MDP'
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 20
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
pbc         = xyz
DispCorr    = EnerPres
gen_vel     = no
MDP

cat > md.mdp << 'MDP'
integrator  = md
nsteps      = 10000000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 20
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
DispCorr    = EnerPres
gen_vel     = no
MDP

echo ""
echo "--- Step 1: pdb2gmx ---"
gmx pdb2gmx -f "$REDUCED" -o processed.gro -water tip3p -ff amber99sb-ildn -ignh 2>&1 | tail -5

echo ""
echo "--- Step 2: editconf (box) ---"
gmx editconf -f processed.gro -o box.gro -c -d 1.5 -bt cubic 2>&1 | tail -3

echo ""
echo "--- Step 3: solvate ---"
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top 2>&1 | tail -3

echo ""
echo "--- Step 4: ions ---"
gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 | tail -3
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 2>&1 | tail -5

echo ""
echo "--- Step 5: Energy Minimization ---"
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2 2>&1 | tail -3
gmx mdrun -deffnm em -ntomp "$NTOMP" 2>&1 | tail -3
echo "  EM potential energy:"
grep "Potential Energy" em.log | tail -1 || echo "  (check em.log)"

echo ""
echo "--- Step 6: NVT Equilibration (100 ps) ---"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2 2>&1 | tail -3
gmx mdrun -deffnm nvt -ntomp "$NTOMP" 2>&1 | tail -3

echo ""
echo "--- Step 7: NPT Equilibration (100 ps) ---"
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2 2>&1 | tail -3
gmx mdrun -deffnm npt -ntomp "$NTOMP" 2>&1 | tail -3

echo ""
echo "--- Step 8: Production MD (20 ns) ---"
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2 2>&1 | tail -3
echo "  Launching production run..."
nohup gmx mdrun -deffnm md -ntomp "$NTOMP" > md_production.log 2>&1 &
MD_PID=$!
echo "  Production MD PID: $MD_PID"
echo "  Log: md_production.log"
echo "  Trajectory: md.xtc (20 ns, 10,000,000 steps)"
echo ""
echo "=== Pipeline complete. Production MD running in background. ==="
echo "  Monitor: tail -f md_production.log"
echo "  End: $(date)"
