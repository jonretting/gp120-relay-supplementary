#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Virtual screening: 120+ diverse fragments against gp120 relay pocket.
Docks each compound at BOTH relay pocket and control site for selectivity.

Usage:
  python3 virtual_screen.py                    # Run all compounds
  python3 virtual_screen.py --start 0 --end 60 # Run subset (for splitting across machines)
"""
import subprocess, os, sys, json, argparse, time
from pathlib import Path

WORK = Path.home() / "md-relay-pocket"
INPUT = WORK / "input"
SCREEN = WORK / "virtual_screen"
os.makedirs(SCREEN, exist_ok=True)

RELAY_RESIDUES = list(range(465, 473))

# ============================================================================
# FRAGMENT LIBRARY -- 4 Tiers
# ============================================================================

LIBRARY = [
    # --- TIER A: HIV therapeutics (positive controls) ---
    {"name": "BMS-626529", "smiles": "COc1cnc(C(=O)N2CCN(c3cccc4[nH]nc(C)c34)CC2)c(OC[C@@H]2CC(F)(F)CO2)c1", "tier": "A", "class": "attachment_inhibitor"},
    {"name": "Maraviroc", "smiles": "CC(C)c1nnc(C)n1[C@@H]1C[C@@H](n2cc(C(=O)NC3CC4CCCCC4C3)nn2)C1", "tier": "A", "class": "ccr5_antagonist"},
    {"name": "Cenicriviroc", "smiles": "CC(C)OC(=O)Oc1ccc(Cn2c(=S)[nH]c3ccc(CC(=O)Nc4ccc(C)cc4OC)cc32)cc1", "tier": "A", "class": "ccr5_antagonist"},
    {"name": "Vicriviroc", "smiles": "CC(=O)N1CC(c2ccc3c(c2)OCO3)CC1c1nc(C(F)(F)F)cc(=O)[nH]1", "tier": "A", "class": "ccr5_antagonist"},
    {"name": "Aplaviroc", "smiles": "OC(=O)CCC(=O)N(Cc1ccc(Cl)cc1)C1CCN(Cc2ccccc2)CC1", "tier": "A", "class": "ccr5_antagonist"},
    {"name": "Temsavir", "smiles": "CC(=O)OC(C)OC(=O)N1CCN(c2cccc3[nH]nc(C)c23)CC1c1ncc(OC)cc1OC[C@@H]1CC(F)(F)CO1", "tier": "A", "class": "attachment_inhibitor_prodrug"},
    {"name": "NBD-556", "smiles": "O=[N+]([O-])c1ccc(NC(=O)c2cc(Cl)c(Cl)cc2Cl)c(O)c1", "tier": "A", "class": "cd4_mimic"},
    {"name": "BMS-378806", "smiles": "COc1cnc(C(=O)N2CCN(c3cccc4[nH]nc(C)c34)CC2)c(OCC2CCCO2)c1", "tier": "A", "class": "attachment_inhibitor"},

    # --- TIER B: PPI inhibitors (mechanistically relevant) ---
    {"name": "ABT-263_Navitoclax", "smiles": "CC1(C)CCC(c2ccc(Cl)cc2)=C(c2ccc(OCCN3CCN(c4ccc(OC(=O)NS(=O)(=O)c5ccc(NCC6CCOCC6)c([N+](=O)[O-])c5)cc4)CC3)cc2)C1", "tier": "B", "class": "ppi_inhibitor"},
    {"name": "ABT-199_Venetoclax", "smiles": "CC1(C)CCC(c2ccc(Cl)cc2)=C(c2ccc(OCCN3CCN(c4ccc(OC(=O)NS(=O)(=O)c5ccc(NCC6CCOCC6)c([N+](=O)[O-])c5)cc4)CC3)cc2)C1", "tier": "B", "class": "ppi_inhibitor"},
    {"name": "Nutlin-3a", "smiles": "COc1ccc(C2=NC(c3ccc(Cl)cc3)C(c3ccc(Cl)cc3)N2C(=O)N2CCCC2)c(OC(C)C)c1", "tier": "B", "class": "ppi_inhibitor"},
    {"name": "AT-406", "smiles": "CNC(=O)c1cc(Oc2ccc3[nH]c(C(=O)NC4CCC(N5CCN(C)CC5)CC4)cc3c2)cc(C(C)(C)C)c1O", "tier": "B", "class": "ppi_inhibitor"},
    {"name": "GDC-0199", "smiles": "CC(C)c1cc(-c2cc3cc(NS(=O)(=O)c4ccc(NCC5CCOCC5)c([N+](=O)[O-])c4)ccc3[nH]2)no1", "tier": "B", "class": "ppi_inhibitor"},
    {"name": "Idasanutlin", "smiles": "COc1cc(F)c(C2=NC(CC(=O)O)(c3ccc(Cl)cc3)C(c3ccc(Cl)cc3F)N2C(=O)C2CCC2)c(F)c1", "tier": "B", "class": "ppi_inhibitor"},

    # --- TIER C: Diverse drug-like fragments (ZINC/ChEMBL) ---
    {"name": "Frag_Indole", "smiles": "c1ccc2[nH]ccc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Benzimidazole", "smiles": "c1ccc2[nH]cnc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Pyrazole", "smiles": "c1cc[nH]n1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Imidazole", "smiles": "c1c[nH]cn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Triazole", "smiles": "c1cn[nH]n1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Piperazine", "smiles": "C1CNCCN1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Morpholine", "smiles": "C1COCCN1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Piperidine", "smiles": "C1CCNCC1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Pyrrolidine", "smiles": "C1CCNC1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Thiophene", "smiles": "c1ccsc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Furan", "smiles": "c1ccoc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Pyridine", "smiles": "c1ccncc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Quinoline", "smiles": "c1ccc2ncccc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Isoquinoline", "smiles": "c1ccc2cnccc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Naphthalene", "smiles": "c1ccc2ccccc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Phenol", "smiles": "Oc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Aniline", "smiles": "Nc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Benzoic_acid", "smiles": "OC(=O)c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Benzamide", "smiles": "NC(=O)c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Acetanilide", "smiles": "CC(=O)Nc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Nicotinamide", "smiles": "NC(=O)c1cccnc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Isonicotinamide", "smiles": "NC(=O)c1ccncc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Sulfonamide", "smiles": "NS(=O)(=O)c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Catechol", "smiles": "Oc1ccccc1O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Guanidine", "smiles": "NC(=N)N", "tier": "C", "class": "fragment"},
    {"name": "Frag_Urea_phenyl", "smiles": "NC(=O)Nc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Hydantoin", "smiles": "O=C1NC(=O)CN1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Barbiturate", "smiles": "O=C1CC(=O)NC(=O)N1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Oxazole", "smiles": "c1cocn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Thiazole", "smiles": "c1cscn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Tetrazole", "smiles": "c1nn[nH]n1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Aminopyridine", "smiles": "Nc1ccccn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Hydroxypyridine", "smiles": "Oc1ccccn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Methylimidazole", "smiles": "Cc1c[nH]cn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Dimethylpyrazole", "smiles": "Cc1cc(C)[nH]n1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Benzothiazole", "smiles": "c1ccc2scnc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Benzoxazole", "smiles": "c1ccc2ocnc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Indazole", "smiles": "c1ccc2[nH]ncc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Purine", "smiles": "c1nc2[nH]cnc2c(=O)[nH]1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Coumarin", "smiles": "O=c1ccc2ccccc2o1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Chromone", "smiles": "O=c1ccoc2ccccc12", "tier": "C", "class": "fragment"},
    {"name": "Frag_Acridine", "smiles": "c1ccc2nc3ccccc3cc2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Carbazole", "smiles": "c1ccc2c(c1)[nH]c1ccccc12", "tier": "C", "class": "fragment"},
    {"name": "Frag_Phenothiazine", "smiles": "c1ccc2c(c1)Nc1ccccc1S2", "tier": "C", "class": "fragment"},
    {"name": "Frag_Flavone", "smiles": "O=c1cc(-c2ccccc2)oc2ccccc12", "tier": "C", "class": "fragment"},
    {"name": "Frag_Xanthine", "smiles": "O=c1[nH]c(=O)c2[nH]cnc2[nH]1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Caffeine", "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Theophylline", "smiles": "Cn1c(=O)c2[nH]cnc2n(C)c1=O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Salicylic_acid", "smiles": "OC(=O)c1ccccc1O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Anthranilic_acid", "smiles": "Nc1ccccc1C(=O)O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Picolinic_acid", "smiles": "OC(=O)c1ccccn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Nicotinic_acid", "smiles": "OC(=O)c1cccnc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Tryptamine", "smiles": "NCCc1c[nH]c2ccccc12", "tier": "C", "class": "fragment"},
    {"name": "Frag_Histamine", "smiles": "NCCc1c[nH]cn1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Tyramine", "smiles": "NCCc1ccc(O)cc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Dopamine", "smiles": "NCCc1ccc(O)c(O)c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Serotonin", "smiles": "NCCc1c[nH]c2ccc(O)cc12", "tier": "C", "class": "fragment"},
    {"name": "Frag_Melatonin", "smiles": "COc1ccc2[nH]cc(CCNC(C)=O)c2c1", "tier": "C", "class": "fragment"},
    {"name": "Frag_GABA", "smiles": "NCCCC(=O)O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Glycine", "smiles": "NCC(=O)O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Proline", "smiles": "OC(=O)C1CCCN1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Phenylalanine", "smiles": "NC(Cc1ccccc1)C(=O)O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Tryptophan", "smiles": "NC(Cc1c[nH]c2ccccc12)C(=O)O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Arginine_analog", "smiles": "NC(CCCNC(=N)N)C(=O)O", "tier": "C", "class": "fragment"},
    {"name": "Frag_Biphenyl", "smiles": "c1ccc(-c2ccccc2)cc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Fluorobenzene", "smiles": "Fc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Chlorobenzene", "smiles": "Clc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Toluene", "smiles": "Cc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Anisole", "smiles": "COc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Nitrobenzene", "smiles": "O=[N+]([O-])c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Cyanobenzene", "smiles": "N#Cc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Benzaldehyde", "smiles": "O=Cc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Acetophenone", "smiles": "CC(=O)c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Methyl_benzoate", "smiles": "COC(=O)c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Phenylacetic_acid", "smiles": "OC(=O)Cc1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Cinnamic_acid", "smiles": "OC(=O)/C=C/c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Stilbene", "smiles": "C(=C/c1ccccc1)\\c1ccccc1", "tier": "C", "class": "fragment"},
    {"name": "Frag_Chalcone", "smiles": "O=C(/C=C/c1ccccc1)c1ccccc1", "tier": "C", "class": "fragment"},

    # --- TIER D: Negative controls (should NOT score well) ---
    {"name": "Neg_ATP", "smiles": "c1nc(N)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C3O)c2n1", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Glucose", "smiles": "OCC1OC(O)C(O)C(O)C1O", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Cholesterol", "smiles": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Sucrose", "smiles": "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Citric_acid", "smiles": "OC(=O)CC(O)(CC(=O)O)C(=O)O", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Glycerol", "smiles": "OCC(O)CO", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Ethanol", "smiles": "CCO", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Acetone", "smiles": "CC(C)=O", "tier": "D", "class": "negative_control"},
    {"name": "Neg_DMSO", "smiles": "CS(C)=O", "tier": "D", "class": "negative_control"},
    {"name": "Neg_Urea", "smiles": "NC(N)=O", "tier": "D", "class": "negative_control"},
]

def run_cmd(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r

def parse_vina_output(log_path):
    results = []
    in_results = False
    with open(log_path) as f:
        for line in f:
            if "-----+-----" in line:
                in_results = True
                continue
            if in_results and line.strip():
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        results.append({
                            "mode": int(parts[0]),
                            "affinity_kcal_mol": float(parts[1]),
                            "rmsd_lb": float(parts[2]),
                            "rmsd_ub": float(parts[3])
                        })
                    except ValueError:
                        if results:
                            break
    return results

def prepare_receptor_once(pdb_path, center, label, radius=20.0):
    cx, cy, cz = center
    rec_dir = SCREEN / "receptors"
    os.makedirs(rec_dir, exist_ok=True)
    rec_pdb = rec_dir / "{}.pdb".format(label)
    rec_pdbqt = rec_dir / "{}.pdbqt".format(label)

    if rec_pdbqt.exists():
        return rec_pdbqt

    lines = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            if ((x - cx)**2 + (y - cy)**2 + (z - cz)**2)**0.5 < radius:
                lines.append(line)
    lines.append("END\n")
    with open(rec_pdb, "w") as f:
        f.writelines(lines)
    run_cmd("obabel {} -O {} -xr -h".format(rec_pdb, rec_pdbqt))
    return rec_pdbqt

def pocket_center(pdb_path, chain, res_range):
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] == chain and int(line[22:26].strip()) in res_range:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))

def prepare_ligand(compound, lig_dir):
    name = compound["name"].replace(" ", "_").replace("/", "_")
    smi_file = lig_dir / "{}.smi".format(name)
    mol2_file = lig_dir / "{}.mol2".format(name)
    pdbqt_file = lig_dir / "{}.pdbqt".format(name)

    if pdbqt_file.exists():
        return pdbqt_file

    with open(smi_file, "w") as f:
        f.write("{} {}\n".format(compound["smiles"], name))

    r = run_cmd("obabel {} -O {} --gen3d -h 2>/dev/null".format(smi_file, mol2_file))
    if not mol2_file.exists():
        return None

    r = run_cmd("obabel {} -O {} 2>/dev/null".format(mol2_file, pdbqt_file))
    if not pdbqt_file.exists():
        return None

    return pdbqt_file

def dock_one(receptor_pdbqt, ligand_pdbqt, center, out_dir, label):
    cx, cy, cz = center
    out_pdbqt = out_dir / "{}_docked.pdbqt".format(label)
    log_file = out_dir / "{}.log".format(label)

    cmd = ("vina --receptor {} --ligand {} --out {} "
           "--center_x {:.3f} --center_y {:.3f} --center_z {:.3f} "
           "--size_x 22 --size_y 22 --size_z 22 "
           "--exhaustiveness 16 --num_modes 5 2>&1 | tee {}").format(
        receptor_pdbqt, ligand_pdbqt, out_pdbqt, cx, cy, cz, log_file)

    run_cmd(cmd)
    if log_file.exists():
        return parse_vina_output(str(log_file))
    return []

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=len(LIBRARY))
    parser.add_argument("--pdb", default=str(INPUT / "hiv_6meo.pdb"))
    args = parser.parse_args()

    pdb = args.pdb
    compounds = LIBRARY[args.start:args.end]
    total = len(compounds)

    print("=" * 60)
    print("VIRTUAL SCREENING: {} compounds (index {}-{})".format(total, args.start, args.end))
    print("=" * 60)

    # Compute pocket centers
    relay_center = pocket_center(pdb, "G", RELAY_RESIDUES)
    ctrl_center = pocket_center(pdb, "G", list(range(415, 426)))

    # Prepare receptors (once)
    relay_rec = prepare_receptor_once(pdb, relay_center, "relay")
    ctrl_rec = prepare_receptor_once(pdb, ctrl_center, "control")
    print("  Relay receptor: {}".format(relay_rec))
    print("  Control receptor: {}".format(ctrl_rec))

    # Prepare ligand directory
    lig_dir = SCREEN / "ligands"
    os.makedirs(lig_dir, exist_ok=True)

    dock_dir = SCREEN / "docking"
    os.makedirs(dock_dir, exist_ok=True)

    results = []
    failed = []

    for i, compound in enumerate(compounds):
        idx = args.start + i
        name = compound["name"]
        print("\n--- [{}/{}] {} (tier {}) ---".format(i + 1, total, name, compound["tier"]))

        lig_pdbqt = prepare_ligand(compound, lig_dir)
        if lig_pdbqt is None:
            print("  FAILED: Could not prepare ligand")
            failed.append(name)
            continue

        safe_name = name.replace(" ", "_").replace("/", "_")

        # Dock at relay pocket
        relay_r = dock_one(relay_rec, lig_pdbqt, relay_center, dock_dir,
                          "{}_relay".format(safe_name))

        # Dock at control site
        ctrl_r = dock_one(ctrl_rec, lig_pdbqt, ctrl_center, dock_dir,
                         "{}_control".format(safe_name))

        relay_best = relay_r[0]["affinity_kcal_mol"] if relay_r else None
        ctrl_best = ctrl_r[0]["affinity_kcal_mol"] if ctrl_r else None
        delta = None
        if relay_best is not None and ctrl_best is not None:
            delta = round(relay_best - ctrl_best, 3)

        entry = {
            "index": idx,
            "name": name,
            "tier": compound["tier"],
            "class": compound["class"],
            "smiles": compound["smiles"],
            "relay_best_kcal_mol": relay_best,
            "control_best_kcal_mol": ctrl_best,
            "delta_kcal_mol": delta,
            "relay_modes": relay_r[:3],
            "control_modes": ctrl_r[:3],
            "selective": delta is not None and delta <= -1.0,
            "druggable": relay_best is not None and relay_best <= -6.0
        }
        results.append(entry)

        status = ""
        if relay_best is not None:
            status += "relay={:.1f}".format(relay_best)
        if ctrl_best is not None:
            status += " ctrl={:.1f}".format(ctrl_best)
        if delta is not None:
            status += " delta={:.1f}".format(delta)
            if delta <= -1.0:
                status += " SELECTIVE"
        if relay_best is not None and relay_best <= -6.0:
            status += " DRUGGABLE"
        print("  {}".format(status))

    # Save results
    output = {
        "experiment": "Virtual screening against gp120 C5 relay pocket",
        "date": "2026-07-16",
        "software": "AutoDock Vina v1.2.7",
        "structure": "PDB 6MEO",
        "total_compounds": total,
        "successful": len(results),
        "failed": failed,
        "results": sorted(results, key=lambda x: x["relay_best_kcal_mol"] or 0),
        "summary": {
            "druggable_count": sum(1 for r in results if r["druggable"]),
            "selective_count": sum(1 for r in results if r["selective"]),
            "tier_breakdown": {}
        }
    }

    for tier in ["A", "B", "C", "D"]:
        tier_r = [r for r in results if r["tier"] == tier]
        if tier_r:
            affs = [r["relay_best_kcal_mol"] for r in tier_r if r["relay_best_kcal_mol"] is not None]
            output["summary"]["tier_breakdown"][tier] = {
                "count": len(tier_r),
                "mean_affinity": round(sum(affs)/len(affs), 2) if affs else None,
                "best_affinity": round(min(affs), 2) if affs else None,
                "druggable": sum(1 for r in tier_r if r["druggable"]),
                "selective": sum(1 for r in tier_r if r["selective"])
            }

    results_path = SCREEN / "virtual_screen_results.json"
    with open(results_path, "w") as f:
        json.dump(output, f, indent=2)

    print("\n" + "=" * 60)
    print("VIRTUAL SCREENING COMPLETE")
    print("=" * 60)
    print("  Total: {}  Successful: {}  Failed: {}".format(total, len(results), len(failed)))
    print("  Druggable (relay < -6.0): {}".format(output["summary"]["druggable_count"]))
    print("  Selective (delta < -1.0): {}".format(output["summary"]["selective_count"]))
    print("\n  Top 10 by relay affinity:")
    for r in output["results"][:10]:
        d = "d={:.1f}".format(r["delta_kcal_mol"]) if r["delta_kcal_mol"] is not None else ""
        print("    {:.1f} kcal/mol  {}  {} [{}]".format(
            r["relay_best_kcal_mol"], d, r["name"], r["tier"]))
    print("\n  Results: {}".format(results_path))

if __name__ == "__main__":
    main()
