import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Relaunch MMPBSA.py with corrected input file."""
import subprocess
import os

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
os.chdir(WORKDIR)

# Fix the input file: use explicit residue numbers instead of "within 6"
# Residues 1-105 in MMPBSA numbering
# Relay pocket = residues 69-76 (gp120 465-472)
# CD4 interface = residues 83-105
# Surrounding gp120 = key contact residues
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
print_res="1-105"
/
"""
with open('mmpbsa_direct.in', 'w') as f:
    f.write(mmpbsa_in)
print("Fixed input file written")

# Launch MMPBSA.py
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
print(f"Launching: {' '.join(cmd)}")

with open('mmpbsa_direct.log', 'w') as logf:
    proc = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT)
    print(f"PID: {proc.pid}")
    proc.wait()
    print(f"Exit code: {proc.returncode}")

with open('mmpbsa_direct.log') as f:
    log = f.read()

if proc.returncode == 0:
    print("\n=== SUCCESS ===")
    if os.path.exists('FINAL_RESULTS_MMPBSA.dat'):
        with open('FINAL_RESULTS_MMPBSA.dat') as f:
            print(f.read()[:3000])
else:
    print("\n=== FAILED ===")
    print(log[-2000:])
