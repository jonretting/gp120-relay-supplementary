#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Phase 2: Molecular Dynamics -- Reduced system around relay pocket
PDB: 6MEO | Focus: Chain G residues 465-472 + 15 Å shell

Strategy: Extract a sphere around the relay pocket, solvate, minimize,
equilibrate, and run 20ns production MD.
"""
import subprocess, os, sys
from pathlib import Path

WORK = Path.home() / "md-relay-pocket"
INPUT = WORK / "input"
MD = WORK / "md"
MDP = WORK / "mdp"

RELAY_RESIDUES = list(range(465, 473))
SHELL_RADIUS = 15.0  # Å around pocket center

def run(cmd, cwd=None, **kw):
    print(f"  >> {cmd}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd or MD, **kw)
    if r.returncode != 0:
        print(f"  STDERR: {r.stderr[:800]}")
    if r.stdout.strip():
        last_lines = r.stdout.strip().split("\n")[-5:]
        for l in last_lines:
            print(f"     {l}")
    return r

def compute_pocket_center(pdb_path):
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            resnum = int(line[22:26].strip())
            if chain == "G" and resnum in RELAY_RESIDUES:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    print(f"  Pocket center: ({cx:.2f}, {cy:.2f}, {cz:.2f}) from {len(xs)} atoms")
    return cx, cy, cz

def extract_reduced_system(pdb_path, cx, cy, cz):
    """Extract all residues with at least one atom within SHELL_RADIUS of pocket center.
    Keep whole residues to avoid dangling bonds."""
    # First pass: find residues with atoms in range
    close_residues = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            dist = ((x - cx)**2 + (y - cy)**2 + (z - cz)**2) ** 0.5
            if dist < SHELL_RADIUS:
                chain = line[21]
                resnum = int(line[22:26].strip())
                close_residues.add((chain, resnum))

    # Second pass: extract full residues
    reduced_pdb = MD / "reduced_system.pdb"
    n_atoms = 0
    with open(pdb_path) as f, open(reduced_pdb, "w") as out:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                resnum = int(line[22:26].strip())
                if (chain, resnum) in close_residues:
                    out.write(line)
                    n_atoms += 1
            elif line.startswith("TER") or line.startswith("END"):
                out.write(line)
        out.write("END\n")

    chains = set(cr[0] for cr in close_residues)
    print(f"  Reduced system: {n_atoms} atoms, {len(close_residues)} residues")
    print(f"  Chains present: {sorted(chains)}")
    return reduced_pdb

def write_mdp_files():
    """Write GROMACS .mdp parameter files."""
    os.makedirs(MDP, exist_ok=True)

    # Energy minimization
    with open(MDP / "em.mdp", "w") as f:
        f.write("""; Energy minimization
integrator  = steep
emtol       = 500.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
""")

    # NVT equilibration (100 ps)
    with open(MDP / "nvt.mdp", "w") as f:
        f.write("""; NVT equilibration
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 500
nstlog      = 500
continuation = no
gen_vel     = yes
gen_temp    = 300
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
define      = -DPOSRES
""")

    # NPT equilibration (100 ps)
    with open(MDP / "npt.mdp", "w") as f:
        f.write("""; NPT equilibration
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 500
nstlog      = 500
continuation = yes
gen_vel     = no
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
define      = -DPOSRES
""")

    # Production MD (20 ns)
    with open(MDP / "md.mdp", "w") as f:
        f.write("""; Production MD -- 20 ns
integrator  = md
nsteps      = 10000000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = yes
gen_vel     = no
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
""")

    # Ions MDP
    with open(MDP / "ions.mdp", "w") as f:
        f.write("""; Ions
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
""")

    print("  MDP files written")

def setup_gromacs(reduced_pdb):
    """Run GROMACS system preparation."""

    # pdb2gmx -- generate topology
    # Use CHARMM27 (available in standard GROMACS) with TIP3P water
    # Select force field interactively is hard, so try AMBER99SB-ILDN which is well-tested
    r = run(f"echo '6\n1' | gmx pdb2gmx -f {reduced_pdb} -o processed.gro -water tip3p -ignh -merge all",
            cwd=MD)
    if r.returncode != 0:
        # Try with CHARMM27
        r = run(f"echo '2\n1' | gmx pdb2gmx -f {reduced_pdb} -o processed.gro -water tip3p -ignh -merge all",
                cwd=MD)
    if r.returncode != 0:
        print("  ERROR: pdb2gmx failed -- may need manual force field selection")
        print("  Available force fields:")
        run("gmx pdb2gmx -f /dev/null -ff none 2>&1 | head -20", cwd=MD)
        return False

    # Create box
    r = run("gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt dodecahedron", cwd=MD)
    if r.returncode != 0:
        return False

    # Solvate
    r = run("gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top", cwd=MD)
    if r.returncode != 0:
        return False

    # Add ions
    r = run(f"gmx grompp -f {MDP}/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 5", cwd=MD)
    if r.returncode != 0:
        return False

    r = run("echo 'SOL' | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15",
            cwd=MD)
    if r.returncode != 0:
        return False

    print("  System preparation complete")
    return True

def run_minimization():
    """Energy minimization."""
    print("\n--- Energy Minimization ---")
    r = run(f"gmx grompp -f {MDP}/em.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 5", cwd=MD)
    if r.returncode != 0:
        return False
    r = run("gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 8", cwd=MD)
    return r.returncode == 0

def run_nvt():
    """NVT equilibration."""
    print("\n--- NVT Equilibration (100 ps) ---")
    r = run(f"gmx grompp -f {MDP}/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 5", cwd=MD)
    if r.returncode != 0:
        return False
    r = run("gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8", cwd=MD)
    return r.returncode == 0

def run_npt():
    """NPT equilibration."""
    print("\n--- NPT Equilibration (100 ps) ---")
    r = run(f"gmx grompp -f {MDP}/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 5", cwd=MD)
    if r.returncode != 0:
        return False
    r = run("gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8", cwd=MD)
    return r.returncode == 0

def run_production():
    """Production MD -- 20 ns (runs in background)."""
    print("\n--- Production MD (20 ns) ---")
    r = run(f"gmx grompp -f {MDP}/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 5", cwd=MD)
    if r.returncode != 0:
        return False
    print("  Starting 20 ns production run (background)...")
    print("  Estimated time: 12-36 hours on 8 CPU cores")
    r = run("nohup gmx mdrun -deffnm md -ntmpi 1 -ntomp 8 > md_production.log 2>&1 &", cwd=MD)
    return True

def main():
    pdb_path = INPUT / "hiv_6meo.pdb"
    if not pdb_path.exists():
        print(f"ERROR: PDB not found at {pdb_path}")
        sys.exit(1)

    os.makedirs(MD, exist_ok=True)

    print("=" * 60)
    print("PHASE 2: MOLECULAR DYNAMICS SETUP")
    print("  Reduced system: 15 Å shell around relay pocket")
    print("  Target: 20 ns production MD")
    print("=" * 60)

    # Compute pocket center
    print("\n--- Computing pocket center ---")
    cx, cy, cz = compute_pocket_center(pdb_path)

    # Extract reduced system
    print("\n--- Extracting reduced system ---")
    reduced_pdb = extract_reduced_system(pdb_path, cx, cy, cz)

    # Write MDP files
    print("\n--- Writing MDP parameter files ---")
    write_mdp_files()

    # GROMACS setup
    print("\n--- GROMACS system preparation ---")
    if not setup_gromacs(reduced_pdb):
        print("\nERROR: System preparation failed")
        sys.exit(1)

    # Minimization
    if not run_minimization():
        print("\nERROR: Minimization failed")
        sys.exit(1)

    # NVT
    if not run_nvt():
        print("\nERROR: NVT failed")
        sys.exit(1)

    # NPT
    if not run_npt():
        print("\nERROR: NPT failed")
        sys.exit(1)

    # Production
    run_production()

    print("\n" + "=" * 60)
    print("MD SETUP COMPLETE")
    print("  Production run started in background")
    print("  Monitor: tail -f <md_dir>/md_production.log")
    print("  Check progress: gmx check -f <md_dir>/md.xtc")
    print("=" * 60)

if __name__ == "__main__":
    main()
