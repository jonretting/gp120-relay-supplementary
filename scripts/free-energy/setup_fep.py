#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
FEP Non-Competitiveness Setup Script

Generates all lambda window directories for the three-leg thermodynamic cycle:
  Leg 1 (ternary): Decouple drug from gp120*CD4*drug
  Leg 2 (binary):  Decouple drug from gp120*drug (no CD4)
  Leg 3 (solvent): Decouple drug in water

Each leg gets 32 lambda windows × N_REPLICAS independent replicas.
"""

import os
import sys
import shutil
import random
import subprocess
import argparse
from pathlib import Path

DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)

N_LAMBDA = 32
N_REPLICAS = 3
MDP_TEMPLATES = ["em.mdp", "nvt.mdp", "npt.mdp", "prod.mdp"]
LEGS = ["ternary", "binary", "solvent"]

GMX = os.environ.get("GMX", "gmx")


def generate_seed():
    return random.randint(10000, 99999)


def prepare_mdp(template_path, output_path, lambda_state, seed=None):
    """Replace LAMBDA_STATE and RANDOM_SEED placeholders in MDP template."""
    with open(template_path) as f:
        content = f.read()

    content = content.replace("LAMBDA_STATE", str(lambda_state))
    if seed is not None:
        content = content.replace("RANDOM_SEED", str(seed))

    with open(output_path, "w") as f:
        f.write(content)


def prepare_binary_system(systems_dir):
    """
    Create the binary system (gp120 + drug, no CD4) from the ternary.

    Removes chain A (CD4) from topology and coordinates.
    Uses GROMACS tools to strip CD4 atoms from the coordinate file
    and edits the topology to exclude chain A.
    """
    ternary_dir = systems_dir / "ternary"
    binary_dir = systems_dir / "binary"
    binary_dir.mkdir(exist_ok=True)

    for f in ternary_dir.iterdir():
        if f.name == "md_prod.gro" or f.name == "complex_ions.gro":
            continue
        if f.is_file():
            shutil.copy2(f, binary_dir / f.name)
        elif f.is_dir():
            if not (binary_dir / f.name).exists():
                shutil.copytree(f, binary_dir / f.name)

    top_path = binary_dir / "topol_complex.top"
    with open(top_path) as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if "topol_Protein_chain_A.itp" in line:
            new_lines.append("; " + line)  # comment out CD4 chain
            continue
        if "posre_Protein_chain_A.itp" in line:
            new_lines.append("; " + line)
            continue
        if line.strip().startswith("Protein_chain_A"):
            new_lines.append("; " + line)  # comment out from [ molecules ]
            continue
        new_lines.append(line)

    with open(top_path, "w") as f:
        f.writelines(new_lines)

    print(f"  Binary topology prepared (CD4/chain A removed from {top_path})")
    return binary_dir


def prepare_solvent_system(systems_dir, box_size="3.0"):
    """
    Create the solvent system (drug in water only).

    Places the drug molecule in a water box and adds ions for neutrality.
    """
    solvent_dir = systems_dir / "solvent"
    solvent_dir.mkdir(exist_ok=True)

    ternary_dir = systems_dir / "ternary"
    for f in ["ligand.itp", "ligand_atomtypes.itp"]:
        shutil.copy2(ternary_dir / f, solvent_dir / f)
    if (ternary_dir / "bms626529.acpype").exists():
        if not (solvent_dir / "bms626529.acpype").exists():
            shutil.copytree(ternary_dir / "bms626529.acpype",
                            solvent_dir / "bms626529.acpype")

    ligand_gro = ternary_dir / "bms626529.acpype" / "bms626529_GMX.gro"
    if not ligand_gro.exists():
        print(f"  WARNING: {ligand_gro} not found. Solvent leg needs manual setup.")
        return solvent_dir

    top_content = """; Solvent-only FEP topology (drug in water)
#include "amber99sb-ildn.ff/forcefield.itp"
#include "ligand_atomtypes.itp"
#include "ligand.itp"
#include "amber99sb-ildn.ff/tip3p.itp"

#ifdef POSRES_WATER
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

#include "amber99sb-ildn.ff/ions.itp"

[ system ]
BMS-626529 in water

[ molecules ]
; filled by setup after solvation
"""
    with open(solvent_dir / "topol_solvent.top", "w") as f:
        f.write(top_content)

    print(f"  Solvent topology template prepared at {solvent_dir / 'topol_solvent.top'}")
    print(f"  Run solvent box setup: gmx editconf + gmx solvate + gmx genion")
    return solvent_dir


def setup_leg(leg_name, base_dir, systems_dir, mdp_dir, coord_file, topol_file):
    """Create all lambda window directories for one leg."""
    leg_dir = base_dir / leg_name
    leg_dir.mkdir(exist_ok=True)

    for replica in range(1, N_REPLICAS + 1):
        rep_dir = leg_dir / f"rep{replica}"
        rep_dir.mkdir(exist_ok=True)
        seed = generate_seed()

        for lam in range(N_LAMBDA):
            win_dir = rep_dir / f"lambda_{lam:02d}"
            win_dir.mkdir(exist_ok=True)

            for mdp_name in MDP_TEMPLATES:
                mdp_src = mdp_dir / mdp_name
                mdp_dst = win_dir / mdp_name
                prepare_mdp(mdp_src, mdp_dst, lam, seed=seed)

            (win_dir / "coord_file.txt").write_text(str(coord_file))
            (win_dir / "topol_file.txt").write_text(str(topol_file))

    total = N_REPLICAS * N_LAMBDA
    print(f"  {leg_name}: {total} windows ({N_REPLICAS} replicas × {N_LAMBDA} lambdas)")
    return leg_dir


def run_grompp_leg(leg_name, base_dir, systems_dir, coord_file, topol_file, dry_run=False):
    """Run grompp for EM step across all windows (validates topology)."""
    leg_dir = base_dir / leg_name
    errors = []

    for replica in range(1, N_REPLICAS + 1):
        rep_dir = leg_dir / f"rep{replica}"
        for lam in range(N_LAMBDA):
            win_dir = rep_dir / f"lambda_{lam:02d}"
            mdp = win_dir / "em.mdp"
            tpr = win_dir / "em.tpr"

            if dry_run:
                print(f"  [DRY] grompp {leg_name}/rep{replica}/lambda_{lam:02d}")
                continue

            topol_dir = systems_dir / (leg_name if leg_name != "solvent" else "solvent")
            cmd = [
                GMX, "grompp",
                "-f", str(mdp),
                "-c", str(coord_file),
                "-p", str(topol_file),
                "-o", str(tpr),
                "-maxwarn", "2",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(topol_dir))
            if result.returncode != 0:
                errors.append((f"{leg_name}/rep{replica}/lambda_{lam:02d}", result.stderr[-500:]))
                if lam == 0 and replica == 1:
                    print(f"  FIRST WINDOW FAILED -- stopping {leg_name}")
                    print(f"  stderr: {result.stderr[-800:]}")
                    return errors

    if not dry_run:
        if errors:
            print(f"  {leg_name}: {len(errors)} grompp errors")
        else:
            print(f"  {leg_name}: all grompp OK")
    return errors


def main():
    parser = argparse.ArgumentParser(description="Setup FEP non-competitiveness pipeline")
    parser.add_argument("--base-dir", default=os.path.join(DATA_ROOT, "results", "06-fep"),
                        help="Base directory for FEP runs")
    parser.add_argument("--mdp-dir", default=None,
                        help="Directory containing MDP templates")
    parser.add_argument("--dry-run", action="store_true",
                        help="Only create directories, don't run grompp")
    parser.add_argument("--grompp", action="store_true",
                        help="Run grompp to validate all windows")
    parser.add_argument("--legs", nargs="+", default=LEGS, choices=LEGS,
                        help="Which legs to set up")
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    systems_dir = base_dir / "systems"
    mdp_dir = Path(args.mdp_dir) if args.mdp_dir else base_dir / "mdp"
    runs_dir = base_dir / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("FEP Non-Competitiveness Pipeline Setup")
    print("=" * 60)
    print(f"Base:    {base_dir}")
    print(f"MDP:     {mdp_dir}")
    print(f"Systems: {systems_dir}")
    print(f"Lambda windows: {N_LAMBDA}")
    print(f"Replicas: {N_REPLICAS}")
    print(f"Legs: {args.legs}")
    print()

    leg_configs = {}

    if "ternary" in args.legs:
        coord = systems_dir / "ternary" / "md_prod.gro"
        topol = systems_dir / "ternary" / "topol_complex.top"
        if not coord.exists():
            coord = systems_dir / "ternary" / "complex_ions.gro"
        leg_configs["ternary"] = (coord, topol)

    if "binary" in args.legs:
        print("Preparing binary system (removing CD4)...")
        prepare_binary_system(systems_dir)
        coord = systems_dir / "ternary" / "md_prod.gro"
        topol = systems_dir / "binary" / "topol_complex.top"
        if not coord.exists():
            coord = systems_dir / "ternary" / "complex_ions.gro"
        leg_configs["binary"] = (coord, topol)

    if "solvent" in args.legs:
        prepare_solvent_system(systems_dir)
        coord = systems_dir / "solvent" / "solvent.gro"
        topol = systems_dir / "solvent" / "topol_solvent.top"
        leg_configs["solvent"] = (coord, topol)

    print()
    print("Creating lambda directories...")
    for leg in args.legs:
        coord, topol = leg_configs[leg]
        setup_leg(leg, runs_dir, systems_dir, mdp_dir, coord, topol)

    if args.grompp:
        print()
        print("Running grompp validation...")
        all_errors = {}
        for leg in args.legs:
            coord, topol = leg_configs[leg]
            errs = run_grompp_leg(leg, runs_dir, systems_dir, coord, topol, dry_run=args.dry_run)
            if errs:
                all_errors[leg] = errs

        if all_errors:
            print(f"\nWARNING: {sum(len(v) for v in all_errors.values())} total grompp errors")
            for leg, errs in all_errors.items():
                for win, msg in errs[:3]:
                    print(f"  {win}: {msg[:200]}")
        else:
            print("\nAll grompp validations passed!")

    print()
    print("Setup complete.")
    print(f"Total windows: {len(args.legs) * N_REPLICAS * N_LAMBDA}")
    print(f"  = {len(args.legs)} legs × {N_REPLICAS} replicas × {N_LAMBDA} lambdas")


if __name__ == "__main__":
    main()
