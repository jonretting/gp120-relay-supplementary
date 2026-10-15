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
Computational Alanine Scanning of HIV-1 gp120 Relay Pocket (465-472)

For each of the 8 relay pocket residues:
  1. Generates a mutant PDB (residue -> Ala) using BioPython
  2. Computes contact-weighted GNM coupling (relay <-> CD4) using ProDy
     - Spring constants are proportional to inter-residue heavy-atom contacts
     - Ala mutants have fewer contacts -> weaker springs -> altered dynamics
  3. Re-docks BMS-626529 and Vicriviroc using AutoDock Vina
  4. Compares to wild-type to quantify fitness cost

The "high-fitness target" argument: if mutating a pocket residue to escape
a drug also destroys the CD4-to-gp41 relay signal, the virus is trapped.
"""

import json
import os
import sys
import subprocess
import shutil
import numpy as np
from pathlib import Path

PDB_FILE = os.path.expanduser("~/hiv_6meo.pdb")
OUTPUT_DIR = os.path.expanduser("~/alanine_scanning")
RELAY_RESIDUES = [465, 466, 467, 468, 469, 470, 471, 472]
RELAY_NAMES = ["ARG465", "PRO466", "GLY467", "GLY468", "GLY469", "ASP470", "MET471", "ARG472"]
GP120_CHAIN = "G"
CD4_CHAIN = "A"
CD4_BINDING_RESIDUES = [25, 39, 40, 43, 45, 46]
POCKET_CENTER_APPROX = [134.6, 184.6, 136.0]
VINA_BOX_SIZE = 22
CONTACT_CUTOFF = 4.5  # heavy-atom contact distance in Angstroms

os.makedirs(OUTPUT_DIR, exist_ok=True)


def mutate_to_ala(pdb_path, chain_id, resid, output_path):
    """Mutate a residue to alanine by stripping side-chain atoms beyond CB."""
    from Bio.PDB import PDBParser, PDBIO, Select

    class AlanineMutantSelect(Select):
        def accept_atom(self, atom):
            res = atom.get_parent()
            ch = res.get_parent()
            if ch.id == chain_id and res.id[1] == resid:
                backbone_atoms = {"N", "CA", "C", "O", "CB"}
                return atom.name in backbone_atoms
            return True

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("wt", pdb_path)

    target_res = None
    for chain in structure[0]:
        if chain.id == chain_id:
            for res in chain:
                if res.id[1] == resid:
                    target_res = res
                    break

    if target_res is None:
        print(f"  WARNING: Residue {chain_id}:{resid} not found!")
        return False, ""

    original_resname = target_res.resname
    n_atoms_before = len(list(target_res.get_atoms()))
    kept = sum(1 for a in target_res.get_atoms() if a.name in {"N", "CA", "C", "O", "CB"})
    removed = n_atoms_before - kept
    print(f"  Mutating {chain_id}:{original_resname}{resid} -> ALA ({removed} side-chain atoms removed)")

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, AlanineMutantSelect())
    fix_resname_to_ala(output_path, chain_id, resid)
    return True, original_resname


def fix_resname_to_ala(pdb_path, chain_id, resid):
    """Post-process PDB to change residue name to ALA."""
    lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                ch = line[21]
                try:
                    rid = int(line[22:26].strip())
                except ValueError:
                    rid = -1
                if ch == chain_id and rid == resid:
                    line = line[:17] + "ALA" + line[20:]
            lines.append(line)
    with open(pdb_path, "w") as f:
        f.writelines(lines)


def count_interresidue_contacts(structure, cutoff=CONTACT_CUTOFF):
    """
    Count heavy-atom contacts between all residue pairs using KDTree.
    Returns dict: (chain_i, resid_i, chain_j, resid_j) -> contact_count
    """
    from scipy.spatial import cKDTree

    coords = []
    atom_info = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                for atom in res:
                    if atom.element == "H":
                        continue
                    coords.append(atom.get_vector().get_array())
                    atom_info.append((chain.id, res.id[1]))
        break  # first model only

    coords = np.array(coords)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(cutoff)

    contacts = {}
    for i, j in pairs:
        ch_i, rid_i = atom_info[i]
        ch_j, rid_j = atom_info[j]
        if ch_i == ch_j and rid_i == rid_j:
            continue
        key = (ch_i, rid_i, ch_j, rid_j)
        rev_key = (ch_j, rid_j, ch_i, rid_i)
        if rev_key in contacts:
            contacts[rev_key] += 1
        else:
            contacts[key] = contacts.get(key, 0) + 1

    return contacts


def compute_contact_weighted_gnm(pdb_path, label=""):
    """
    Compute GNM with contact-weighted spring constants.
    k_ij = number of heavy-atom contacts between residues i and j.
    This captures side-chain removal effects that standard Cα GNM misses.
    """
    import prody
    from Bio.PDB import PDBParser

    prody.confProDy(verbosity="none")

    bp_parser = PDBParser(QUIET=True)
    bp_struct = bp_parser.get_structure("s", pdb_path)

    print(f"    Counting inter-residue heavy-atom contacts...")
    contacts = count_interresidue_contacts(bp_struct, CONTACT_CUTOFF)
    print(f"    Found {len(contacts)} residue-residue contact pairs")

    structure = prody.parsePDB(pdb_path)
    if structure is None:
        return None

    calphas = structure.select("calpha")
    if calphas is None:
        return None

    n = len(calphas)
    chain_ids = calphas.getChids()
    resnums = calphas.getResnums()

    res_to_idx = {}
    for i in range(n):
        res_to_idx[(chain_ids[i], resnums[i])] = i

    kirchhoff = np.zeros((n, n))
    for (ch_i, rid_i, ch_j, rid_j), count in contacts.items():
        idx_i = res_to_idx.get((ch_i, rid_i))
        idx_j = res_to_idx.get((ch_j, rid_j))
        if idx_i is not None and idx_j is not None:
            kirchhoff[idx_i, idx_j] = -count
            kirchhoff[idx_j, idx_i] = -count

    for i in range(n):
        kirchhoff[i, i] = -np.sum(kirchhoff[i, :])

    eigenvalues, eigenvectors = np.linalg.eigh(kirchhoff)

    n_modes = min(20, n - 1)
    nonzero = eigenvalues[1:n_modes+1]
    modes = eigenvectors[:, 1:n_modes+1]

    inv_eigenvalues = 1.0 / nonzero
    cross_corr = np.zeros((n, n))
    for k in range(n_modes):
        outer = np.outer(modes[:, k], modes[:, k])
        cross_corr += inv_eigenvalues[k] * outer

    diag = np.sqrt(np.diag(cross_corr))
    norm = np.outer(diag, diag)
    norm[norm == 0] = 1.0
    cross_corr_norm = cross_corr / norm

    sq_fluct = np.diag(cross_corr)

    relay_indices = []
    cd4_indices = []
    for i in range(n):
        if chain_ids[i] == GP120_CHAIN and resnums[i] in RELAY_RESIDUES:
            relay_indices.append(i)
        if chain_ids[i] == CD4_CHAIN and resnums[i] in CD4_BINDING_RESIDUES:
            cd4_indices.append(i)

    if not relay_indices or not cd4_indices:
        print(f"  WARNING: relay={len(relay_indices)}, cd4={len(cd4_indices)}")
        return None

    couplings = []
    for ri in relay_indices:
        for ci in cd4_indices:
            couplings.append(abs(cross_corr_norm[ri, ci]))

    intra_pocket_corr = []
    for i, ri in enumerate(relay_indices):
        for j, rj in enumerate(relay_indices):
            if j > i:
                intra_pocket_corr.append(cross_corr_norm[ri, rj])

    mean_coupling = float(np.mean(couplings))
    mean_all_fluct = float(np.mean(sq_fluct))
    mean_relay_fluct = float(np.mean([sq_fluct[i] for i in relay_indices]))

    total_relay_contacts = 0
    for (ch_i, rid_i, ch_j, rid_j), count in contacts.items():
        if (ch_i == GP120_CHAIN and rid_i in RELAY_RESIDUES) or \
           (ch_j == GP120_CHAIN and rid_j in RELAY_RESIDUES):
            total_relay_contacts += count

    return {
        "mean_coupling": mean_coupling,
        "max_coupling": float(np.max(couplings)),
        "min_coupling": float(np.min(couplings)),
        "n_relay": len(relay_indices),
        "n_cd4": len(cd4_indices),
        "mean_relay_fluctuation": mean_relay_fluct,
        "mean_all_fluctuation": mean_all_fluct,
        "rigidity_ratio": mean_relay_fluct / mean_all_fluct if mean_all_fluct > 0 else None,
        "total_relay_contacts": total_relay_contacts,
        "mean_intra_pocket_corr": float(np.mean(intra_pocket_corr)) if intra_pocket_corr else None,
        "total_contact_pairs": len(contacts),
    }


def extract_clean_receptor(pdb_path, output_path, chains=("A", "G")):
    """Extract specified chains, remove water/heteroatoms for docking."""
    from Bio.PDB import PDBParser, PDBIO, Select

    class ChainSelect(Select):
        def accept_chain(self, chain):
            return chain.id in chains
        def accept_residue(self, residue):
            return residue.id[0] == " "

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, ChainSelect())


def prepare_pdbqt(pdb_path, pdbqt_path):
    """Convert PDB to PDBQT for Vina docking."""
    clean_pdb = pdbqt_path.replace(".pdbqt", "_clean.pdb")
    extract_clean_receptor(pdb_path, clean_pdb)
    cmd = f"obabel {clean_pdb} -O {pdbqt_path} -xr 2>&1"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if not os.path.exists(pdbqt_path) or os.path.getsize(pdbqt_path) < 100:
        print(f"    obabel stderr: {result.stderr[:200]}")
        return False
    size = os.path.getsize(pdbqt_path)
    print(f"    PDBQT created: {size:,} bytes")
    return True


def run_vina_docking(receptor_pdbqt, ligand_pdbqt, out_pdbqt, label=""):
    """Run AutoDock Vina docking and return best affinity."""
    cmd = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(POCKET_CENTER_APPROX[0]),
        "--center_y", str(POCKET_CENTER_APPROX[1]),
        "--center_z", str(POCKET_CENTER_APPROX[2]),
        "--size_x", str(VINA_BOX_SIZE),
        "--size_y", str(VINA_BOX_SIZE),
        "--size_z", str(VINA_BOX_SIZE),
        "--exhaustiveness", "16",
        "--out", out_pdbqt,
        "--num_modes", "5",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    best_affinity = None
    for line in result.stdout.split("\n"):
        line = line.strip()
        parts = line.split()
        if len(parts) >= 4:
            try:
                mode = int(parts[0])
                affinity = float(parts[1])
                if mode == 1:
                    best_affinity = affinity
                    break
            except (ValueError, IndexError):
                continue

    if best_affinity is None:
        print(f"    Vina parse issue for {label}.")
        print(f"    stdout (last 5 lines): {chr(10).join(result.stdout.strip().split(chr(10))[-5:])}")
        if result.stderr:
            print(f"    stderr: {result.stderr[:300]}")

    return best_affinity


def main():
    print("=" * 80)
    print("COMPUTATIONAL ALANINE SCANNING -- HIV-1 gp120 RELAY POCKET")
    print("=" * 80)
    print(f"PDB: {PDB_FILE}")
    print(f"Relay residues: {RELAY_NAMES}")
    print(f"Contact cutoff: {CONTACT_CUTOFF} A")
    print(f"Output: {OUTPUT_DIR}")
    print()

    bms_pdbqt = os.path.join(OUTPUT_DIR, "bms626529.pdbqt")
    vic_pdbqt = os.path.join(OUTPUT_DIR, "vicriviroc.pdbqt")

    for name, dst in [("bms626529", bms_pdbqt), ("vicriviroc", vic_pdbqt)]:
        if not os.path.exists(dst):
            src = os.path.expanduser(fos.path.join(DATA_ROOT, "{name}.pdbqt"))
            if os.path.exists(src):
                shutil.copy2(src, dst)
                print(f"  Copied ligand: {src}")

    has_ligands = os.path.exists(bms_pdbqt) and os.path.exists(vic_pdbqt)
    if not has_ligands:
        print("  WARNING: Ligand PDBQT files not found. Skipping docking.")

    # --- STEP 1: Wild-type contact-weighted GNM ---
    print("\n--- STEP 1: Wild-Type Contact-Weighted GNM ---")
    wt_gnm = compute_contact_weighted_gnm(PDB_FILE, label="wild-type")
    if wt_gnm:
        print(f"  WT relay-CD4 coupling: {wt_gnm['mean_coupling']:.4f}")
        print(f"  WT rigidity ratio: {wt_gnm['rigidity_ratio']:.4f}")
        print(f"  WT total relay contacts: {wt_gnm['total_relay_contacts']}")
        print(f"  WT intra-pocket correlation: {wt_gnm['mean_intra_pocket_corr']:.4f}")
    else:
        print("  FATAL: Wild-type GNM failed!")
        sys.exit(1)

    # --- STEP 2: Wild-type docking ---
    wt_bms_score = None
    wt_vic_score = None

    if has_ligands:
        print("\n--- STEP 2: Wild-Type Docking ---")
        wt_pdbqt = os.path.join(OUTPUT_DIR, "wt_receptor.pdbqt")

        existing = os.path.join(DATA_ROOT, "structures", "gp120_gp41_relay_receptor.pdbqt")
        if os.path.exists(existing):
            shutil.copy2(existing, wt_pdbqt)
            print(f"  Using existing receptor PDBQT ({os.path.getsize(existing):,} bytes)")
        else:
            print("  Preparing receptor PDBQT from PDB...")
            if not prepare_pdbqt(PDB_FILE, wt_pdbqt):
                print("  FAILED to prepare WT receptor PDBQT")
                has_ligands = False

        if has_ligands:
            wt_bms_out = os.path.join(OUTPUT_DIR, "wt_bms626529_out.pdbqt")
            print("  Docking BMS-626529 against WT...")
            wt_bms_score = run_vina_docking(wt_pdbqt, bms_pdbqt, wt_bms_out, "WT-BMS")
            print(f"  WT BMS-626529: {wt_bms_score} kcal/mol")

            wt_vic_out = os.path.join(OUTPUT_DIR, "wt_vicriviroc_out.pdbqt")
            print("  Docking Vicriviroc against WT...")
            wt_vic_score = run_vina_docking(wt_pdbqt, vic_pdbqt, wt_vic_out, "WT-Vic")
            print(f"  WT Vicriviroc: {wt_vic_score} kcal/mol")
    else:
        print("\n--- STEP 2: Skipping docking ---")

    # --- STEP 3: Mutant scanning ---
    print("\n--- STEP 3: Alanine Scanning (8 Mutants) ---")

    results = {
        "wild_type": {
            "gnm": wt_gnm,
            "docking": {"bms626529": wt_bms_score, "vicriviroc": wt_vic_score},
        },
        "mutants": {},
    }

    for resid, resname in zip(RELAY_RESIDUES, RELAY_NAMES):
        print(f"\n  === {resname} -> ALA ===")

        mutant_pdb = os.path.join(OUTPUT_DIR, f"mutant_{resname}_to_ALA.pdb")
        success, orig_name = mutate_to_ala(PDB_FILE, GP120_CHAIN, resid, mutant_pdb)
        if not success:
            print(f"  FAILED to generate mutant")
            continue

        print(f"  Computing contact-weighted GNM...")
        mut_gnm = compute_contact_weighted_gnm(mutant_pdb, label=f"{resname}->Ala")

        coupling_delta = None
        coupling_pct = None
        contacts_delta = None

        if mut_gnm:
            coupling_delta = mut_gnm["mean_coupling"] - wt_gnm["mean_coupling"]
            coupling_pct = (coupling_delta / wt_gnm["mean_coupling"]) * 100
            contacts_delta = mut_gnm["total_relay_contacts"] - wt_gnm["total_relay_contacts"]
            rig_delta = mut_gnm["rigidity_ratio"] - wt_gnm["rigidity_ratio"]
            ipc_delta = (mut_gnm["mean_intra_pocket_corr"] - wt_gnm["mean_intra_pocket_corr"]) if mut_gnm["mean_intra_pocket_corr"] else 0

            print(f"  Coupling: {mut_gnm['mean_coupling']:.4f} (Δ: {coupling_delta:+.4f}, {coupling_pct:+.1f}%)")
            print(f"  Relay contacts: {mut_gnm['total_relay_contacts']} (Δ: {contacts_delta:+d})")
            print(f"  Rigidity ratio: {mut_gnm['rigidity_ratio']:.4f} (Δ: {rig_delta:+.4f})")
            print(f"  Intra-pocket corr: {mut_gnm['mean_intra_pocket_corr']:.4f} (Δ: {ipc_delta:+.4f})")

        mut_bms_score = None
        mut_vic_score = None
        bms_delta = None
        vic_delta = None

        if has_ligands:
            mut_pdbqt = os.path.join(OUTPUT_DIR, f"mutant_{resname}_receptor.pdbqt")
            print(f"  Preparing mutant receptor PDBQT...")
            if prepare_pdbqt(mutant_pdb, mut_pdbqt):
                bms_out = os.path.join(OUTPUT_DIR, f"mutant_{resname}_bms626529_out.pdbqt")
                print(f"  Docking BMS-626529...")
                mut_bms_score = run_vina_docking(mut_pdbqt, bms_pdbqt, bms_out, f"{resname}-BMS")
                bms_delta = (mut_bms_score - wt_bms_score) if (mut_bms_score and wt_bms_score) else None
                print(f"  BMS: {mut_bms_score} kcal/mol" +
                      (f" (Δ: {bms_delta:+.3f})" if bms_delta is not None else ""))

                vic_out = os.path.join(OUTPUT_DIR, f"mutant_{resname}_vicriviroc_out.pdbqt")
                print(f"  Docking Vicriviroc...")
                mut_vic_score = run_vina_docking(mut_pdbqt, vic_pdbqt, vic_out, f"{resname}-Vic")
                vic_delta = (mut_vic_score - wt_vic_score) if (mut_vic_score and wt_vic_score) else None
                print(f"  Vic: {mut_vic_score} kcal/mol" +
                      (f" (Δ: {vic_delta:+.3f})" if vic_delta is not None else ""))

        results["mutants"][resname] = {
            "residue_id": resid,
            "original_resname": orig_name,
            "mutation": f"{orig_name}{resid}Ala",
            "gnm": mut_gnm,
            "gnm_coupling_delta": coupling_delta,
            "gnm_coupling_pct_change": coupling_pct,
            "contacts_lost": -contacts_delta if contacts_delta else None,
            "docking": {
                "bms626529": mut_bms_score,
                "bms626529_delta": bms_delta,
                "vicriviroc": mut_vic_score,
                "vicriviroc_delta": vic_delta,
            }
        }

    # --- STEP 4: Summary ---
    print("\n" + "=" * 80)
    print("ALANINE SCANNING RESULTS")
    print("=" * 80)

    header = f"{'Mutation':<18} {'Contacts':>10} {'ΔContacts':>10} {'Coupling':>10} {'ΔCoupling':>10} {'Δ%':>8}"
    if has_ligands and wt_bms_score:
        header += f" {'BMS':>8} {'ΔBMS':>8} {'Vic':>8} {'ΔVic':>8}"
    print(f"\n{header}")
    print("-" * len(header))

    wt_line = f"{'Wild-type':<18} {wt_gnm['total_relay_contacts']:>10d} {'--':>10} {wt_gnm['mean_coupling']:>10.4f} {'--':>10} {'--':>8}"
    if has_ligands and wt_bms_score:
        wt_line += f" {wt_bms_score:>8.3f} {'--':>8} {wt_vic_score:>8.3f} {'--':>8}"
    print(wt_line)

    for resname in RELAY_NAMES:
        if resname not in results["mutants"]:
            continue
        m = results["mutants"][resname]
        gnm = m["gnm"]
        dock = m["docking"]

        if gnm is None:
            continue

        contacts = gnm["total_relay_contacts"]
        dc = m.get("contacts_lost", 0)
        coupling = gnm["mean_coupling"]
        dcoup = m["gnm_coupling_delta"] or 0
        dpct = m["gnm_coupling_pct_change"] or 0

        marker = ""
        if abs(dpct) > 10:
            marker = " *** HIGH COST ***"
        elif abs(dpct) > 5:
            marker = " ** moderate **"

        line = f"{resname+'Ala':<18} {contacts:>10d} {-dc if dc else 0:>+10d} {coupling:>10.4f} {dcoup:>+10.4f} {dpct:>+7.1f}%"
        if has_ligands and wt_bms_score:
            bms = dock.get("bms626529") or 0
            dbms = dock.get("bms626529_delta") or 0
            vic = dock.get("vicriviroc") or 0
            dvic = dock.get("vicriviroc_delta") or 0
            line += f" {bms:>8.3f} {dbms:>+8.3f} {vic:>8.3f} {dvic:>+8.3f}"
        line += marker
        print(line)

    # --- Fitness cost analysis ---
    print("\n--- FITNESS COST ANALYSIS ---")
    high_cost = []
    for resname in RELAY_NAMES:
        if resname not in results["mutants"]:
            continue
        m = results["mutants"][resname]
        pct = m.get("gnm_coupling_pct_change")
        lost = m.get("contacts_lost")
        if pct is not None:
            high_cost.append((resname, pct, lost or 0))

    high_cost.sort(key=lambda x: abs(x[1]), reverse=True)

    print("\nRanked by coupling disruption:")
    for name, pct, lost in high_cost:
        severity = "HIGH" if abs(pct) > 10 else ("MODERATE" if abs(pct) > 5 else "LOW")
        print(f"  {name}->Ala: {pct:+.1f}% coupling, {lost:+d} contacts lost -- {severity} fitness cost")

    fitness_trapped = [x for x in high_cost if abs(x[1]) > 5]
    if fitness_trapped:
        print(f"\n  {len(fitness_trapped)} of 8 relay residues show >5% coupling disruption upon mutation.")
        print("  The virus cannot mutate these positions without impairing its entry mechanism.")
    else:
        print("\n  Coupling changes are modest, suggesting the relay signal propagates")
        print("  through backbone topology rather than side-chain contacts alone.")

    # --- Save results ---
    output_file = os.path.join(OUTPUT_DIR, "alanine_scanning_results.json")
    results["metadata"] = {
        "pdb": "6MEO",
        "chain": GP120_CHAIN,
        "relay_residues": RELAY_NAMES,
        "method": "Contact-weighted GNM (spring constants = heavy-atom contact count)",
        "contact_cutoff_angstroms": CONTACT_CUTOFF,
        "gnm_modes": 20,
        "docking_method": "AutoDock Vina 1.2.7, exhaustiveness 16, box 22x22x22 A",
        "cd4_reference_residues": CD4_BINDING_RESIDUES,
    }

    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to: {output_file}")

    print("\n" + "=" * 80)
    print("ALANINE SCANNING COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
