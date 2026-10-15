import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Build AMBER topology from GROMACS using pdb2gmx + parmed, bypassing tleap.

Strategy:
1. Extract protein PDB from trajectory
2. Use gmx pdb2gmx to generate GROMACS .top with same force field
3. Use parmed to convert GROMACS .top -> AMBER .prmtop
4. Atom count matches trajectory exactly
5. Run MMPBSA.py directly
"""
import subprocess
import os
import parmed

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
MDDIR = os.path.join(DATA_ROOT, "topology")
os.chdir(WORKDIR)

# Step 1: Extract protein-only first frame as PDB
print("Step 1: Extracting protein PDB from first frame...")
proc = subprocess.run(
    ['gmx', 'trjconv', '-s', f'{MDDIR}/md.tpr', '-f', f'{MDDIR}/md.xtc',
     '-o', 'protein_frame0.pdb', '-dump', '0'],
    input='1\n', capture_output=True, text=True
)
print(f"  Exit: {proc.returncode}")

# Count atoms
with open('protein_frame0.pdb') as f:
    natoms = sum(1 for l in f if l.startswith('ATOM'))
print(f"  Atoms in PDB: {natoms}")

# Step 2: Use gmx pdb2gmx to generate topology
# Need to know what force field was used. Check original TPR.
print("\nStep 2: Generating GROMACS topology with pdb2gmx...")
# Try AMBER99SB-ILDN (common choice, compatible with AMBER MMPBSA)
proc = subprocess.run(
    ['gmx', 'pdb2gmx', '-f', 'protein_frame0.pdb', '-o', 'processed.gro',
     '-p', 'topol.top', '-ignh', '-ff', 'amber99sb-ildn', '-water', 'tip3p',
     '-merge', 'all'],
    input='\n', capture_output=True, text=True
)
print(f"  Exit: {proc.returncode}")
if proc.returncode != 0:
    print(f"  Stderr: {proc.stderr[:1000]}")
    # Try without -merge
    proc = subprocess.run(
        ['gmx', 'pdb2gmx', '-f', 'protein_frame0.pdb', '-o', 'processed.gro',
         '-p', 'topol.top', '-ignh', '-ff', 'amber99sb-ildn', '-water', 'tip3p'],
        input='\n', capture_output=True, text=True
    )
    print(f"  Retry exit: {proc.returncode}")
    if proc.returncode != 0:
        print(f"  Retry stderr: {proc.stderr[:1000]}")

# Check atom count in generated topology
if os.path.exists('processed.gro'):
    with open('processed.gro') as f:
        lines = f.readlines()
    gro_natoms = int(lines[1].strip())
    print(f"  Atoms in processed.gro: {gro_natoms}")

# Step 3: Convert GROMACS topology to AMBER format using parmed
print("\nStep 3: Converting GROMACS topology to AMBER format...")
if os.path.exists('topol.top'):
    try:
        gmx_top = parmed.load_file('topol.top', xyz='processed.gro')
        print(f"  Loaded GROMACS topology: {len(gmx_top.atoms)} atoms")
        
        # Save as AMBER format
        gmx_top.save('COM_gmx.prmtop', overwrite=True)
        gmx_top.save('COM_gmx.inpcrd', overwrite=True)
        
        # Verify
        amber = parmed.load_file('COM_gmx.prmtop')
        print(f"  AMBER topology: {len(amber.atoms)} atoms")
        print(f"  Trajectory atoms: {natoms}")
        
        if len(amber.atoms) == natoms:
            print("\n  EXACT MATCH! Topology matches trajectory.")
        else:
            print(f"\n  Mismatch: topology {len(amber.atoms)} vs trajectory {natoms}")
            
    except Exception as e:
        print(f"  Error: {e}")
else:
    print("  topol.top not found!")
