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
Task 6.2: Non-Competitiveness MM-GBSA
Compare gp120-CD4 binding free energy in drug-bound vs apo ternary complex.
If DeltaDeltaG ~ 0, drug does not compete with CD4 binding.

Apo reference: DeltaG = -38.02 kcal/mol (Task 1.1b)

Receptor: gp120 (chain 0)
Ligand: CD4 (chain 1)
Complex context: drug-bound ternary (gp120+CD4+CCR5+BMS-626529)
"""
import os
import subprocess
import json
import re
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(os.path.join(DATA_ROOT, "results", "04-md-simulations"))
OUT_DIR = WORK_DIR / "mmpbsa_noncompetitive"
OUT_DIR.mkdir(exist_ok=True)
os.chdir(str(WORK_DIR))


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def run(cmd, cwd=None, timeout=7200):
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        timeout=timeout, cwd=str(cwd or WORK_DIR)
    )
    return result.stdout + result.stderr, result.returncode


log("=" * 70)
log("Task 6.2: Non-Competitiveness MM-GBSA (drug-bound ternary)")
log("Measuring gp120-CD4 binding in presence of BMS-626529")
log("Apo reference: DeltaG = -38.02 kcal/mol")
log("=" * 70)

# Step 1: Create chain index file
log("Step 1: Creating chain index groups from GRO file")
gro_file = None
for candidate in ["npt.gro", "nvt.gro", "em.gro"]:
    p = WORK_DIR / candidate
    if p.exists():
        gro_file = p
        break

if gro_file is None:
    log("  No GRO found, generating from TPR...")
    run("echo 0 | gmx editconf -f md.tpr -o system.gro")
    gro_file = WORK_DIR / "system.gro"

log(f"  Using: {gro_file}")

with open(gro_file) as f:
    lines = f.readlines()

n_atoms = int(lines[1].strip())
log(f"  Total atoms in GRO: {n_atoms}")

prev_resnum = -1
chain_id = 0
atoms_by_chain = {}
solvent_names = {"SOL", "NA", "CL", "K", "WAT", "HOH"}

for i, line in enumerate(lines[2:2 + n_atoms]):
    resnum = int(line[0:5].strip())
    resname = line[5:10].strip()
    atom_num = i + 1

    if resname in solvent_names:
        continue

    if resnum < prev_resnum and prev_resnum > 0:
        chain_id += 1

    if chain_id not in atoms_by_chain:
        atoms_by_chain[chain_id] = {"atoms": [], "resnames": set(),
                                     "first_res": resnum}
    atoms_by_chain[chain_id]["atoms"].append(atom_num)
    atoms_by_chain[chain_id]["resnames"].add(resname)
    atoms_by_chain[chain_id]["last_res"] = resnum
    prev_resnum = resnum

chain_names = {}
for cid, info in sorted(atoms_by_chain.items()):
    n_at = len(info["atoms"])
    n_res = info["last_res"]
    has_unl = "UNL" in info["resnames"]
    log(f"  Chain {cid}: {n_at} atoms, res 1-{n_res}, UNL={has_unl}")

    if has_unl and n_at < 200:
        chain_names[cid] = "UNL"
    elif n_res > 350 and cid == 0:
        chain_names[cid] = "gp120"
    elif 100 < n_res < 250 and not has_unl:
        chain_names[cid] = "CD4"
    elif n_res > 250:
        chain_names[cid] = "CCR5"
    else:
        chain_names[cid] = f"chain_{cid}"

log(f"  Chain mapping: {chain_names}")

ndx_file = OUT_DIR / "noncomp.ndx"
with open(ndx_file, "w") as f:
    for cid, info in sorted(atoms_by_chain.items()):
        name = chain_names.get(cid, f"chain_{cid}")
        f.write(f"[ {name} ]\n")
        for j, a in enumerate(info["atoms"]):
            f.write(f"{a:7d}")
            if (j + 1) % 15 == 0:
                f.write("\n")
        if len(info["atoms"]) % 15 != 0:
            f.write("\n")
        f.write("\n")

group_names = [chain_names[cid] for cid in sorted(atoms_by_chain.keys())]
gp120_grp = group_names.index("gp120")
cd4_grp = group_names.index("CD4")
log(f"  gp120 group index: {gp120_grp}, CD4 group index: {cd4_grp}")

# Step 2: PBC-unwrapped trajectory
log("Step 2: Checking PBC-unwrapped trajectory")
unwrapped = WORK_DIR / "md_pbc_whole.xtc"
if not unwrapped.exists():
    log("  PBC unwrapping needed...")
    run("echo '1\n0' | gmx trjconv -f md.xtc -s md.tpr "
        "-o md_pbc_whole.xtc -pbc mol -center")

if unwrapped.exists():
    size_mb = unwrapped.stat().st_size / 1e6
    log(f"  {unwrapped.name}: {size_mb:.1f} MB")
else:
    log("  ERROR: No unwrapped trajectory!")
    exit(1)

# Step 3: Extract trajectory subset (every 5th frame)
log("Step 3: Extracting trajectory subset")
subset = OUT_DIR / "md_subset.xtc"
if not subset.exists():
    run(f"echo '0' | gmx trjconv -f {unwrapped} -s md.tpr "
        f"-o {subset} -skip 5")
log("  Subset trajectory created")

# Step 4: Write MMPBSA input (overall DeltaG only, no decomp)
log("Step 4: Writing MMPBSA input")
mmpbsa_in = OUT_DIR / "noncomp.in"
with open(mmpbsa_in, "w") as f:
    f.write("&general\n"
            "  startframe=1,\n"
            "  endframe=200,\n"
            "  interval=1,\n"
            "  verbose=2,\n"
            "  keep_files=0,\n"
            "/\n"
            "&gb\n"
            "  igb=5,\n"
            "  saltcon=0.150,\n"
            "/\n")

# Step 5: Run gmx_MMPBSA
log("Step 5: Running gmx_MMPBSA (receptor=gp120, ligand=CD4)")
log(f"  Index groups: receptor={gp120_grp} (gp120), ligand={cd4_grp} (CD4)")

cmd = (
    f"cd {OUT_DIR} ; "
    f"gmx_MMPBSA -O "
    f"-i noncomp.in "
    f"-cs {WORK_DIR}/md.tpr "
    f"-ci noncomp.ndx "
    f"-cg {gp120_grp} {cd4_grp} "
    f"-ct md_subset.xtc "
    f"-cp {WORK_DIR}/topol.top "
    f"-nogui "
    f"> mmpbsa_run.log 2>&1"
)
log(f"  Command: {cmd}")
output, rc = run(cmd, cwd=OUT_DIR, timeout=14400)
log(f"  gmx_MMPBSA exit code: {rc}")

# Step 6: Parse results
log("Step 6: Parsing results")
results_file = None
for name in ["FINAL_RESULTS_MMPBSA.dat", "FINAL_RESULTS.dat"]:
    p = OUT_DIR / name
    if p.exists():
        results_file = p
        break

results = {
    "task": "6.2",
    "description": "Non-Competitiveness MM-GBSA (drug-bound ternary)",
    "receptor": "gp120",
    "ligand": "CD4",
    "context": "drug-bound (BMS-626529)",
    "apo_reference_kcal": -38.02,
}

if results_file and results_file.exists():
    with open(results_file) as f:
        content = f.read()
    log(f"  Results file: {results_file.name}")
    log(content[-1500:])

    delta_match = re.search(
        r"ΔTOTAL\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)", content
    )
    if not delta_match:
        delta_match = re.search(
            r"TOTAL\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)", content
        )
    if delta_match:
        dg = float(delta_match.group(1))
        sd_prop = float(delta_match.group(2))
        sd = float(delta_match.group(3))
        results["deltaG_drug_bound"] = dg
        results["sd"] = sd
        results["sd_prop"] = sd_prop
        results["deltaG_apo"] = -38.02
        results["deltadeltaG"] = round(dg - (-38.02), 2)
        log(f"  DeltaG (drug-bound): {dg:.2f} +/- {sd:.2f} kcal/mol")
        log(f"  DeltaG (apo ref): -38.02 kcal/mol")
        log(f"  DeltaDeltaG: {dg - (-38.02):.2f} kcal/mol")
        if abs(dg - (-38.02)) < 10:
            log("  RESULT: Drug does NOT significantly compete with "
                "CD4 binding (non-competitive mechanism)")
        else:
            log(f"  RESULT: Significant difference detected "
                f"({dg - (-38.02):.1f} kcal/mol)")
else:
    log("  WARNING: No results file found!")
    run_log = OUT_DIR / "mmpbsa_run.log"
    if run_log.exists():
        with open(run_log) as f:
            log(f.read()[-2000:])

json_out = OUT_DIR / "noncompetitive_results.json"
with open(json_out, "w") as f:
    json.dump(results, f, indent=2)
log(f"Results saved to: {json_out}")
log("DONE")
