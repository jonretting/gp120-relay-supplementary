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
Launch MMPBSA.py using parmed-converted topologies that match trajectory exactly.

Uses the COM_gmx.prmtop (1715 atoms) which matches the trajectory.
Splits into receptor (gp120, chain G) and ligand (CD4, chain A) topologies.
"""
import parmed
import subprocess
import os

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
MDDIR = os.path.join(DATA_ROOT, "topology")
os.chdir(WORKDIR)

# Load the matched topology
com = parmed.load_file('COM_gmx.prmtop', 'COM_gmx.inpcrd')
print(f"Complex: {len(com.atoms)} atoms, {len(com.residues)} residues")

# Identify chain boundary
# Chain G (gp120) = first 82 residues (resid 94-478 in original)
# Chain A (CD4) = last 23 residues (resid 23-62 in original)
# From GRO inspection: residue numbering resets at chain boundary

res_list = [(r.number, r.name) for r in com.residues]
print(f"\nResidue list (first 5): {res_list[:5]}")
print(f"Residue list (boundary): {res_list[80:84]}")
print(f"Residue list (last 5): {res_list[-5:]}")

# Split at residue index 82 (0-based)
REC_NRES = 82
LIG_NRES = len(com.residues) - REC_NRES

# Count atoms per chain
rec_atom_count = sum(len(list(com.residues[i].atoms)) for i in range(REC_NRES))
lig_atom_count = sum(len(list(com.residues[i].atoms)) for i in range(REC_NRES, len(com.residues)))
print(f"\nReceptor (gp120): {REC_NRES} residues, {rec_atom_count} atoms")
print(f"Ligand (CD4):     {LIG_NRES} residues, {lig_atom_count} atoms")
print(f"Total:            {rec_atom_count + lig_atom_count} atoms")

# Create receptor topology
rec = com.copy(parmed.Structure)
# Strip ligand atoms (everything after receptor)
lig_indices = []
for i in range(REC_NRES, len(com.residues)):
    for a in com.residues[i].atoms:
        lig_indices.append(a.idx + 1)
rec_mask = '@' + ','.join(str(i) for i in lig_indices)

rec_top = parmed.load_file('COM_gmx.prmtop', 'COM_gmx.inpcrd')
rec_top.strip(rec_mask)
rec_top.save('REC_gmx.prmtop', overwrite=True)
print(f"\nREC_gmx.prmtop: {len(rec_top.atoms)} atoms")

# Create ligand topology
rec_indices = []
for i in range(REC_NRES):
    for a in com.residues[i].atoms:
        rec_indices.append(a.idx + 1)
lig_mask = '@' + ','.join(str(i) for i in rec_indices)

lig_top = parmed.load_file('COM_gmx.prmtop', 'COM_gmx.inpcrd')
lig_top.strip(lig_mask)
lig_top.save('LIG_gmx.prmtop', overwrite=True)
print(f"LIG_gmx.prmtop: {len(lig_top.atoms)} atoms")

# Extract protein-only trajectory
print("\nExtracting protein-only trajectory...")
proc = subprocess.run(
    ['gmx', 'trjconv', '-s', f'{MDDIR}/md.tpr', '-f', f'{MDDIR}/md.xtc',
     '-o', 'protein_traj.xtc', '-pbc', 'nojump'],
    input='1\n', capture_output=True, text=True
)
print(f"  trjconv exit: {proc.returncode}")

# Set PB radii in topology (needed for MMPBSA)
print("\nSetting PB radii...")
for top_file in ['COM_gmx.prmtop', 'REC_gmx.prmtop', 'LIG_gmx.prmtop']:
    t = parmed.load_file(top_file)
    action = parmed.tools.actions.changeRadii(t, 'mbondi2')
    action.execute()
    t.save(top_file, overwrite=True)
    print(f"  {top_file}: radii set to mbondi2")

# Create MMPBSA input file
print("\nCreating MMPBSA.py input file...")
mmpbsa_in = """Per-residue MM/GBSA decomposition of gp120-CD4 relay interface
&general
startframe=1, endframe=9999, interval=10,
verbose=2, keep_files=0,
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=3,
print_res="within 6"
/
"""
with open('mmpbsa_direct.in', 'w') as f:
    f.write(mmpbsa_in)

# Build the cpptraj input to create receptor and ligand trajectories
# We need separate trajectories for COM, REC, LIG
print("\nCreating receptor/ligand trajectory extraction scripts...")

# Receptor atom mask (1-based, first rec_atom_count atoms)
rec_atom_range = f'1-{rec_atom_count}'
lig_atom_range = f'{rec_atom_count+1}-{rec_atom_count+lig_atom_count}'

cpptraj_script = f"""parm COM_gmx.prmtop
trajin protein_traj.xtc 1 last 10
strip !(@{rec_atom_range})
trajout rec_traj.nc netcdf
go

clear all
parm COM_gmx.prmtop
trajin protein_traj.xtc 1 last 10
strip !(@{lig_atom_range})
trajout lig_traj.nc netcdf
go

clear all
parm COM_gmx.prmtop
trajin protein_traj.xtc 1 last 10
trajout com_traj.nc netcdf
go
"""
with open('extract_trajs.cpptraj', 'w') as f:
    f.write(cpptraj_script)

# Run cpptraj to extract trajectories
print("Extracting receptor/ligand/complex trajectories with cpptraj...")
proc = subprocess.run(
    ['cpptraj', '-i', 'extract_trajs.cpptraj'],
    capture_output=True, text=True
)
print(f"  cpptraj exit: {proc.returncode}")
if proc.returncode != 0:
    print(f"  stderr: {proc.stderr[:500]}")
    print(f"  stdout: {proc.stdout[-500:]}")

# Check output files
for f in ['com_traj.nc', 'rec_traj.nc', 'lig_traj.nc']:
    if os.path.exists(f):
        size = os.path.getsize(f) / (1024*1024)
        print(f"  {f}: {size:.1f} MB")
    else:
        print(f"  {f}: MISSING")

# Now run MMPBSA.py
print("\n=== Ready to launch MMPBSA.py ===")
cmd = [
    'MMPBSA.py', '-O',
    '-i', 'mmpbsa_direct.in',
    '-cp', 'COM_gmx.prmtop',
    '-rp', 'REC_gmx.prmtop',
    '-lp', 'LIG_gmx.prmtop',
    '-y', 'com_traj.nc',
    '-yr', 'rec_traj.nc',
    '-yl', 'lig_traj.nc',
    '-o', 'FINAL_RESULTS_MMPBSA.dat',
    '-eo', 'FINAL_RESULTS_MMPBSA.csv',
    '-do', 'FINAL_DECOMP_MMPBSA.dat',
    '-deo', 'FINAL_DECOMP_MMPBSA.csv',
]
print(f"Command: {' '.join(cmd)}")
print("\nLaunching MMPBSA.py...")

with open('mmpbsa_direct.log', 'w') as logf:
    proc = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT)
    print(f"PID: {proc.pid}")
    print("Running in foreground... (check mmpbsa_direct.log for progress)")
    proc.wait()
    print(f"Exit code: {proc.returncode}")

# Check results
if proc.returncode == 0:
    print("\n=== MMPBSA.py COMPLETED SUCCESSFULLY ===")
    if os.path.exists('FINAL_RESULTS_MMPBSA.dat'):
        with open('FINAL_RESULTS_MMPBSA.dat') as f:
            print(f.read()[:2000])
else:
    print("\n=== MMPBSA.py FAILED ===")
    with open('mmpbsa_direct.log') as f:
        print(f.read()[-2000:])
