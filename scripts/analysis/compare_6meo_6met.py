#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Compare relay pocket structure between PDB 6MEO and 6MET."""
import math, sys
from pathlib import Path

WORK = Path.home() / "md-relay-pocket" / "input"
RELAY = list(range(465, 473))
CTRL = list(range(415, 426))

def read_atoms(pdb_path, chain, res_range):
    atoms = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            c = line[21]
            rn = int(line[22:26].strip())
            if c == chain and rn in res_range:
                aname = line[12:16].strip()
                resname = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                bf = float(line[60:66].strip())
                atoms[(resname, rn, aname)] = (x, y, z, bf)
    return atoms

def compare_region(label, a_meo, a_met, res_range):
    print("\n=== {} ===".format(label))
    print("  6MEO atoms: {}  6MET atoms: {}".format(len(a_meo), len(a_met)))

    # CA comparison
    print("\n  CA atom positions:")
    print("  {:<12} {:>8} {:>8} {:>8}   {:>8} {:>8} {:>8}   {:>6}  {:>6} {:>6}".format(
        "Residue", "MEO_x", "MEO_y", "MEO_z", "MET_x", "MET_y", "MET_z", "dist", "BF_meo", "BF_met"))

    ca_dists = []
    for rn in res_range:
        for key in sorted(a_meo.keys()):
            if key[1] == rn and key[2] == "CA":
                if key in a_met:
                    m = a_meo[key]
                    t = a_met[key]
                    d = math.sqrt(sum((a - b)**2 for a, b in zip(m[:3], t[:3])))
                    ca_dists.append(d)
                    print("  {}{:<4} CA   {:8.2f} {:8.2f} {:8.2f}   {:8.2f} {:8.2f} {:8.2f}   {:6.3f}  {:6.1f} {:6.1f}".format(
                        key[0][:3], rn, m[0], m[1], m[2], t[0], t[1], t[2], d, m[3], t[3]))
                else:
                    print("  {}{:<4} CA   MISSING IN 6MET".format(key[0][:3], rn))

    if ca_dists:
        print("\n  CA RMSD: {:.3f} A  (mean: {:.3f}, max: {:.3f})".format(
            math.sqrt(sum(d**2 for d in ca_dists) / len(ca_dists)),
            sum(ca_dists) / len(ca_dists),
            max(ca_dists)))

    # All-atom comparison
    common = set(a_meo.keys()) & set(a_met.keys())
    if common:
        all_dists = []
        for k in common:
            m = a_meo[k]
            t = a_met[k]
            d = math.sqrt(sum((a - b)**2 for a, b in zip(m[:3], t[:3])))
            all_dists.append(d)
        print("\n  All-atom: {} common atoms".format(len(common)))
        print("  Mean displacement: {:.3f} A".format(sum(all_dists) / len(all_dists)))
        print("  RMSD: {:.3f} A".format(math.sqrt(sum(d**2 for d in all_dists) / len(all_dists))))
        print("  Max displacement:  {:.3f} A (atom: {})".format(
            max(all_dists),
            sorted(common, key=lambda k: math.sqrt(sum((a-b)**2 for a,b in zip(a_meo[k][:3], a_met[k][:3]))))[-1]))

    # Check 20A shell differences
    return ca_dists

def check_shell(pdb_path, center, radius=20.0):
    """Count atoms in shell and check composition."""
    cx, cy, cz = center
    chains = {}
    total = 0
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            if ((x-cx)**2 + (y-cy)**2 + (z-cz)**2)**0.5 < radius:
                c = line[21]
                chains[c] = chains.get(c, 0) + 1
                total += 1
    return total, chains

def main():
    meo_path = WORK / "hiv_6meo.pdb"
    met_path = WORK / "6MET.pdb"

    print("=" * 70)
    print("STRUCTURAL COMPARISON: PDB 6MEO vs 6MET")
    print("=" * 70)

    # Relay pocket
    meo_relay = read_atoms(meo_path, "G", RELAY)
    met_relay = read_atoms(met_path, "G", RELAY)
    relay_dists = compare_region("RELAY POCKET (G:465-472)", meo_relay, met_relay, RELAY)

    # Control site
    meo_ctrl = read_atoms(meo_path, "G", CTRL)
    met_ctrl = read_atoms(met_path, "G", CTRL)
    ctrl_dists = compare_region("CONTROL SITE (G:415-425)", meo_ctrl, met_ctrl, CTRL)

    # Shell composition
    print("\n=== 20 A SHELL COMPARISON ===")
    # Use relay pocket center from 6MEO
    xs, ys, zs = [], [], []
    with open(meo_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == "G" and int(line[22:26].strip()) in RELAY:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    relay_center = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))

    meo_total, meo_chains = check_shell(meo_path, relay_center)
    met_total, met_chains = check_shell(met_path, relay_center)

    print("  Relay center: ({:.2f}, {:.2f}, {:.2f})".format(*relay_center))
    print("  6MEO shell: {} atoms  chains: {}".format(meo_total, dict(sorted(meo_chains.items()))))
    print("  6MET shell: {} atoms  chains: {}".format(met_total, dict(sorted(met_chains.items()))))

    if meo_total != met_total:
        print("  DIFFERENCE: {} atoms ({:.1f}%)".format(
            abs(meo_total - met_total),
            abs(meo_total - met_total) / max(meo_total, met_total) * 100))

    # Check if 6MET has different chain naming
    print("\n=== CHAIN INVENTORY ===")
    for label, path in [("6MEO", meo_path), ("6MET", met_path)]:
        chains = {}
        with open(path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    c = line[21]
                    chains[c] = chains.get(c, 0) + 1
        print("  {}: {}".format(label, dict(sorted(chains.items()))))

    print("\n=== DIAGNOSIS ===")
    if relay_dists and max(relay_dists) < 0.01:
        print("  WARNING: Relay pocket coordinates are IDENTICAL between 6MEO and 6MET.")
        print("  This means either:")
        print("    1. The two PDBs share the same coordinates for chain G relay region")
        print("    2. The script read the wrong file")
        print("  Check: are 6MEO.pdb and 6MET.pdb actually different files?")
    elif relay_dists and max(relay_dists) < 0.5:
        print("  Relay pocket is nearly identical (max CA shift < 0.5 A)")
        print("  Docking difference likely due to surrounding shell atom positions")
    elif relay_dists:
        print("  Relay pocket shows conformational differences (max CA shift: {:.2f} A)".format(max(relay_dists)))
        print("  This explains docking score variation between structures")

if __name__ == "__main__":
    main()
