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
Transfer BMS-626529 (JYV) from 6MU6 to the equilibrated ternary complex
by superimposing the relay pocket backbone of 6MU6 onto our structure.

Strategy:
1. Extract gp120 C-alpha from 6MU6 near the ligand binding site
2. Extract corresponding C-alpha from our equilibrated structure
3. Compute optimal superposition (Kabsch)
4. Apply the transformation to the ligand coordinates
5. Write the combined receptor+ligand PDB
"""
import os
import json
import numpy as np
from pathlib import Path
from datetime import datetime

DOCK_DIR = Path(os.path.join(DATA_ROOT, "results", "02-docking"))

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def parse_pdb_atoms(filename, record_types=("ATOM",), chain=None, resrange=None):
    """Parse PDB file and return list of atom dicts."""
    atoms = []
    with open(filename) as f:
        for line in f:
            if not line.startswith(record_types):
                continue
            rec = {
                "type": line[0:6].strip(),
                "serial": int(line[6:11].strip()),
                "name": line[12:16].strip(),
                "altloc": line[16],
                "resname": line[17:20].strip(),
                "chain": line[21],
                "resnum": int(line[22:26].strip()),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "line": line,
            }
            if chain and rec["chain"] != chain:
                continue
            if resrange and not (resrange[0] <= rec["resnum"] <= resrange[1]):
                continue
            atoms.append(rec)
    return atoms

def kabsch(P, Q):
    """Kabsch algorithm: find rotation R that minimizes ||P - Q*R||."""
    C = P.T @ Q
    U, S, Vt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return R

log("=" * 70)
log("LIGAND TRANSFER: 6MU6 JYV -> Ternary Complex")
log("=" * 70)

# Load structures
pdb_6mu6 = DOCK_DIR / "6mu6.pdb"
receptor_pdb = DOCK_DIR / "receptor_equilibrated.pdb"

# Step 1: Find the ligand and its surrounding protein in 6MU6
log("Step 1: Analyzing 6MU6 binding site")

# Get JYV ligand atoms
ligand_atoms = parse_pdb_atoms(pdb_6mu6, record_types=("HETATM",))
jyv_atoms = [a for a in ligand_atoms if a["resname"] == "JYV"]

# Group by chain
jyv_chains = {}
for a in jyv_atoms:
    if a["chain"] not in jyv_chains:
        jyv_chains[a["chain"]] = []
    jyv_chains[a["chain"]].append(a)

log(f"  JYV found in chains: {sorted(jyv_chains.keys())}")
for ch, atoms in jyv_chains.items():
    log(f"  Chain {ch}: {len(atoms)} atoms")

# Use the first JYV instance
lig_chain = sorted(jyv_chains.keys())[0]
lig_atoms = jyv_chains[lig_chain]
lig_coords = np.array([[a["x"], a["y"], a["z"]] for a in lig_atoms])
lig_center = lig_coords.mean(axis=0)
log(f"  Using chain {lig_chain}: {len(lig_atoms)} atoms, center = ({lig_center[0]:.1f}, {lig_center[1]:.1f}, {lig_center[2]:.1f})")

# Step 2: Find the gp120 chain in 6MU6 that contains the ligand binding site
log("Step 2: Finding gp120 chain in 6MU6")

all_6mu6_atoms = parse_pdb_atoms(pdb_6mu6, record_types=("ATOM",))

# Find which chain has residues near the ligand
ca_near_lig = []
for a in all_6mu6_atoms:
    if a["name"] != "CA":
        continue
    dist = np.sqrt((a["x"] - lig_center[0])**2 + (a["y"] - lig_center[1])**2 + (a["z"] - lig_center[2])**2)
    if dist < 20.0:
        ca_near_lig.append({**a, "dist": dist})

from collections import Counter
chain_counts = Counter(a["chain"] for a in ca_near_lig)
ref_chain = chain_counts.most_common(1)[0][0]
log(f"  Ligand-proximal chain: {ref_chain}")

# Get ALL CA from this chain for sequence analysis
ref_ca = [a for a in all_6mu6_atoms if a["name"] == "CA" and a["chain"] == ref_chain]
log(f"  Chain {ref_chain}: {len(ref_ca)} CA atoms total")

# Step 3: Find corresponding residues in the equilibrated ternary complex
log("Step 3: Mapping residues between 6MU6 and ternary complex")

# The ternary complex has gp120 as the first chain (no chain ID in GROMACS PDB)
# We need to match by sequence

# Get reference residue sequence from 6MU6
aa3to1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
           "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
           "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

ref_sequence = [(a["resnum"], aa3to1.get(a["resname"], "?")) for a in ref_ca]
log(f"  6MU6 reference residues: {ref_sequence[:10]}...")

# Get target gp120 sequence from receptor PDB
# In GROMACS PDB output, residue numbers restart for each chain.
# gp120 is the first chain. We identify it by finding the first drop in resnum.
target_atoms = parse_pdb_atoms(receptor_pdb, record_types=("ATOM",))

# Find chain boundary: scan CA atoms until resnum drops
all_target_ca = [a for a in target_atoms if a["name"] == "CA"]
gp120_end_idx = len(all_target_ca)
for i in range(1, len(all_target_ca)):
    if all_target_ca[i]["resnum"] <= all_target_ca[i-1]["resnum"]:
        gp120_end_idx = i
        break

target_ca = all_target_ca[:gp120_end_idx]
log(f"  Target gp120 CA: {len(target_ca)} residues (stopped at chain break)")

# Build target sequence
target_sequence = [(a["resnum"], aa3to1.get(a["resname"], "?")) for a in target_ca]
target_seq_str = "".join(r[1] for r in target_sequence)
log(f"  Target gp120 sequence length: {len(target_sequence)}")

# Build reference sequence string for matching
ref_seq_str = "".join(r[1] for r in ref_sequence)
log(f"  Matching sequence fragment: {ref_seq_str[:30]}...")

# Find the offset by looking for common subsequences
# The relay pocket in our structure is at res 372-379 (RPGGGDMR)
# In 6MU6 numbering, the pocket is at different residue numbers
# But the sequence should match

# Find RPGGGDMR in both
pocket_motif = "RPGGGDMR"
idx_target = target_seq_str.find(pocket_motif)
log(f"  RPGGGDMR in target: position {idx_target} (res {target_sequence[idx_target][0] if idx_target >= 0 else '?'})")

# Find in reference chain
ref_chain_all_ca = [a for a in all_6mu6_atoms if a["name"] == "CA" and a["chain"] == ref_chain]
ref_full_seq = "".join(aa3to1.get(a["resname"], "?") for a in ref_chain_all_ca)
idx_ref = ref_full_seq.find(pocket_motif)
log(f"  RPGGGDMR in 6MU6 chain {ref_chain}: position {idx_ref} (res {ref_chain_all_ca[idx_ref]['resnum'] if idx_ref >= 0 else '?'})")

if idx_target < 0 or idx_ref < 0:
    log("  ERROR: Cannot find RPGGGDMR in one of the sequences")
    log(f"  Target seq near end: ...{target_seq_str[-30:]}")
    log(f"  Ref seq near end: ...{ref_full_seq[-30:]}")
    exit(1)

# Compute sequence offset
# target_seq[idx_target] = ref_seq[idx_ref]
# target_resnum = target_sequence[idx_target][0]
# ref_resnum = ref_chain_all_ca[idx_ref]["resnum"]
offset = target_sequence[idx_target][0] - ref_chain_all_ca[idx_ref]["resnum"]

# Actually, we should do a proper residue-by-residue alignment
# Use the CA near ligand from 6MU6 and find matching CA in target
log(f"  Sequence offset (target - ref): {offset}")

# Step 4: Align using ALL BACKBONE ATOMS (N, CA, C, O) of the pocket and flanking residues
# The RPGGGDMR motif is 100% conserved -- use it as the alignment anchor
log("Step 4: Building POCKET-CENTERED alignment using all backbone atoms")

# Pocket residues: ref = 469-476, target = 372-379
# Use pocket ± 10 flanking residues for more points
ref_pocket_start = ref_chain_all_ca[idx_ref]["resnum"]
ref_pocket_end = ref_chain_all_ca[min(idx_ref + 7, len(ref_chain_all_ca)-1)]["resnum"]
flank = 10

log(f"  Ref pocket: res {ref_pocket_start}-{ref_pocket_end}")
log(f"  Using flanking range: ±{flank} residues")

# Get ALL backbone atoms (not just CA) from pocket ± flank
ref_bb_range = range(ref_pocket_start - flank, ref_pocket_end + flank + 1)
tgt_bb_range = range(372 - flank, 379 + flank + 1)

ref_bb_atoms = parse_pdb_atoms(pdb_6mu6, record_types=("ATOM",), chain=ref_chain)
ref_bb = [a for a in ref_bb_atoms if a["name"] in ("N", "CA", "C", "O") and a["resnum"] in ref_bb_range]

tgt_bb_atoms = parse_pdb_atoms(receptor_pdb, record_types=("ATOM",))
# Filter to gp120 only (first 6204 atoms by serial)
tgt_gp120_max_serial = max(a["serial"] for a in target_ca)
tgt_bb = [a for a in tgt_bb_atoms if a["name"] in ("N", "CA", "C", "O") 
          and a["resnum"] in tgt_bb_range and a["serial"] <= tgt_gp120_max_serial]

log(f"  Ref backbone atoms in range: {len(ref_bb)}")
log(f"  Target backbone atoms in range: {len(tgt_bb)}")

# Match by (resnum-offset, atomname) pairs
paired_ref = []
paired_target = []

for r in ref_bb:
    ref_key = (r["resnum"] + offset, r["name"])
    for t in tgt_bb:
        if (t["resnum"], t["name"]) == ref_key:
            # Verify same residue type
            if r["resname"] == t["resname"]:
                paired_ref.append([r["x"], r["y"], r["z"]])
                paired_target.append([t["x"], t["y"], t["z"]])
            break

log(f"  Matched backbone atom pairs: {len(paired_ref)}")

if len(paired_ref) < 8:
    log("  WARNING: Few backbone pairs, checking residue type matching...")
    # Try without residue type check (backbone atoms are universal)
    paired_ref = []
    paired_target = []
    for r in ref_bb:
        ref_key = (r["resnum"] + offset, r["name"])
        for t in tgt_bb:
            if (t["resnum"], t["name"]) == ref_key:
                paired_ref.append([r["x"], r["y"], r["z"]])
                paired_target.append([t["x"], t["y"], t["z"]])
                break
    log(f"  Matched (type-relaxed): {len(paired_ref)}")

if len(paired_ref) < 4:
    log("  ERROR: Too few matching atoms for alignment")
    exit(1)

P = np.array(paired_ref)
Q = np.array(paired_target)

# Center both
P_center = P.mean(axis=0)
Q_center = Q.mean(axis=0)
P_centered = P - P_center
Q_centered = Q - Q_center

# Kabsch
R = kabsch(Q_centered, P_centered)

# Test RMSD
P_aligned = (P_centered @ R) + Q_center
rmsd = np.sqrt(np.mean(np.sum((P_aligned - Q)**2, axis=1)))
log(f"  Alignment RMSD: {rmsd:.2f} Å ({len(paired_ref)} CA pairs)")

# Step 5: Transform ligand coordinates
log("Step 5: Transforming ligand to target frame")

lig_transformed = ((lig_coords - P_center) @ R) + Q_center
log(f"  Original ligand center: ({lig_center[0]:.1f}, {lig_center[1]:.1f}, {lig_center[2]:.1f})")
new_center = lig_transformed.mean(axis=0)
log(f"  Transformed center: ({new_center[0]:.1f}, {new_center[1]:.1f}, {new_center[2]:.1f})")

# Check distance to pocket center in target
pocket_ca = [a for a in target_ca if 372 <= a["resnum"] <= 379]
pocket_center = np.array([[a["x"], a["y"], a["z"]] for a in pocket_ca]).mean(axis=0)
dist_to_pocket = np.linalg.norm(new_center - pocket_center)
log(f"  Distance to target pocket center: {dist_to_pocket:.1f} Å")

# Step 6: Write combined PDB
log("Step 6: Writing combined receptor+ligand PDB")

output_pdb = DOCK_DIR / "ternary_with_ligand.pdb"
with open(output_pdb, "w") as f:
    # Write receptor atoms
    with open(receptor_pdb) as rf:
        for line in rf:
            if line.startswith(("ATOM", "TER")):
                f.write(line)
    
    # Write transformed ligand
    atom_serial = 14150
    for i, atom in enumerate(lig_atoms):
        x, y, z = lig_transformed[i]
        f.write(
            f"HETATM{atom_serial:5d} {atom['name']:>4s} JYV {lig_chain}{999:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
        )
        atom_serial += 1
    f.write("END\n")

log(f"  Written: {output_pdb}")

# Also write just the ligand for reference
ligand_out = DOCK_DIR / "bms626529_positioned.pdb"
with open(ligand_out, "w") as f:
    for i, atom in enumerate(lig_atoms):
        x, y, z = lig_transformed[i]
        f.write(
            f"HETATM{i+1:5d} {atom['name']:>4s} JYV {lig_chain}{999:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
        )
    f.write("END\n")

# Save results
results = {
    "method": "Kabsch superposition from 6MU6 to equilibrated ternary complex",
    "reference_structure": "6MU6",
    "reference_chain": ref_chain,
    "ligand_resname": "JYV (BMS-814508)",
    "ligand_atoms": len(lig_atoms),
    "alignment_ca_pairs": len(paired_ref),
    "alignment_rmsd_angstrom": round(float(rmsd), 3),
    "sequence_offset": int(offset),
    "ligand_center_original": [round(float(c), 1) for c in lig_center],
    "ligand_center_transferred": [round(float(c), 1) for c in new_center],
    "distance_to_pocket_center": round(float(dist_to_pocket), 1),
    "output_files": {
        "receptor_with_ligand": str(output_pdb),
        "ligand_only": str(ligand_out),
    },
}

with open(DOCK_DIR / "ligand_transfer_results.json", "w") as f:
    json.dump(results, f, indent=2)

log(f"\nResults: {json.dumps(results, indent=2)}")
log("=" * 70)
log("DONE")
