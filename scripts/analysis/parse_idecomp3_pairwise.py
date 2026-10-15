import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse idecomp=3 (pairwise) decomposition results from gmx_MMPBSA.

Extracts pairwise interaction energies between receptor (gp120) and
ligand (CD4) residues, identifies the strongest interacting pairs,
and highlights relay pocket involvement.
"""
import sys, json, re
from collections import defaultdict

DECOMP_FILE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA_ROOT, "mmpbsa_idecomp3/FINAL_DECOMP_id3.dat")

RELAY_RESIDUES = set(range(372, 380))  # gp120 RPGGGDMR (GROMACS numbering)
GP120_MAX = 398
CD4_START = 1   # CD4 restarts numbering at 1

def parse_pairwise_decomp(filepath):
    """Parse the DELTAS section of idecomp=3 output."""
    pairs = []
    in_deltas = False
    header_found = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()

            if 'Pairwise' in line and 'DELTAS' in line:
                in_deltas = True
                continue

            if in_deltas and not header_found:
                if 'Resid 1' in line or 'Total' in line:
                    header_found = True
                continue

            if not in_deltas or not header_found:
                continue

            if line.strip() == '' or line.startswith('---') or line.startswith('==='):
                if pairs:
                    break
                continue

            parts = line.split(',') if ',' in line else line.split()
            if len(parts) < 6:
                continue

            try:
                res1_str = parts[0].strip()
                res2_str = parts[1].strip()

                match1 = re.match(r'([A-Z]{3})\s*(\d+)', res1_str)
                match2 = re.match(r'([A-Z]{3})\s*(\d+)', res2_str)
                if not match1 or not match2:
                    if len(parts) >= 10:
                        res1_name = parts[0].strip()
                        res1_num = int(parts[1].strip())
                        res2_name = parts[2].strip()
                        res2_num = int(parts[3].strip())
                        internal = float(parts[4].strip())
                        vdw = float(parts[5].strip())
                        elec = float(parts[6].strip())
                        polar = float(parts[7].strip())
                        nonpolar = float(parts[8].strip())
                        total = float(parts[9].strip())
                    else:
                        continue
                else:
                    res1_name = match1.group(1)
                    res1_num = int(match1.group(2))
                    res2_name = match2.group(1)
                    res2_num = int(match2.group(2))

                    vals = [float(x) for x in parts[2:]]
                    if len(vals) >= 6:
                        internal, vdw, elec, polar, nonpolar, total = vals[:6]
                    elif len(vals) >= 4:
                        vdw, elec, polar, total = vals[0], vals[1], vals[2], vals[3]
                        internal, nonpolar = 0.0, 0.0
                    else:
                        continue

                pairs.append({
                    'res1': res1_name,
                    'res1_num': res1_num,
                    'res2': res2_name,
                    'res2_num': res2_num,
                    'vdw': vdw,
                    'elec': elec,
                    'polar': polar,
                    'total': total,
                })
            except (ValueError, IndexError):
                continue

    return pairs


def classify_pair(p):
    """Determine which residue is receptor (gp120) and which is ligand (CD4)."""
    n1, n2 = p['res1_num'], p['res2_num']
    if n1 <= GP120_MAX and n2 <= GP120_MAX:
        return 'intra-gp120'
    elif n1 > GP120_MAX and n2 > GP120_MAX:
        return 'intra-CD4'
    else:
        return 'interface'


def main():
    print("=" * 70)
    print("PAIRWISE DECOMPOSITION ANALYSIS (idecomp=3)")
    print("=" * 70)

    pairs = parse_pairwise_decomp(DECOMP_FILE)
    if not pairs:
        print(f"ERROR: No pairwise data parsed from {DECOMP_FILE}")
        print("Trying alternative parsing...")
        with open(DECOMP_FILE, 'r') as f:
            content = f.read()
        print(f"File size: {len(content)} bytes")
        print(f"First 2000 chars:\n{content[:2000]}")
        return

    print(f"\nTotal pairwise entries parsed: {len(pairs)}")

    interface = [p for p in pairs if classify_pair(p) == 'interface']
    print(f"Interface pairs (gp120-CD4): {len(interface)}")

    interface_sorted = sorted(interface, key=lambda p: p['total'])

    print("\n--- Top 30 STRONGEST interface pairs (most negative ΔG) ---")
    print(f"{'Rank':>4} {'gp120 Res':>12} {'CD4 Res':>12} {'VdW':>8} {'Elec':>8} {'Polar':>8} {'Total':>8} {'Relay?':>6}")
    print("-" * 75)
    for i, p in enumerate(interface_sorted[:30]):
        if p['res1_num'] <= GP120_MAX:
            gp120_str = f"{p['res1']}{p['res1_num']}"
            cd4_str = f"{p['res2']}{p['res2_num']}"
            gp120_num = p['res1_num']
        else:
            gp120_str = f"{p['res2']}{p['res2_num']}"
            cd4_str = f"{p['res1']}{p['res1_num']}"
            gp120_num = p['res2_num']

        relay = "YES" if gp120_num in RELAY_RESIDUES else ""
        print(f"{i+1:>4} {gp120_str:>12} {cd4_str:>12} {p['vdw']:>8.3f} {p['elec']:>8.3f} {p['polar']:>8.3f} {p['total']:>8.3f} {relay:>6}")

    relay_pairs = [p for p in interface if
                   (p['res1_num'] in RELAY_RESIDUES or p['res2_num'] in RELAY_RESIDUES)]
    relay_sorted = sorted(relay_pairs, key=lambda p: p['total'])

    print(f"\n--- Relay Pocket Pairs (RPGGGDMR, res 372-379) ---")
    print(f"Found {len(relay_pairs)} relay-involved interface pairs")
    print(f"{'gp120 Res':>12} {'CD4 Res':>12} {'VdW':>8} {'Elec':>8} {'Polar':>8} {'Total':>8}")
    print("-" * 60)
    for p in relay_sorted[:20]:
        if p['res1_num'] <= GP120_MAX:
            gp120_str = f"{p['res1']}{p['res1_num']}"
            cd4_str = f"{p['res2']}{p['res2_num']}"
        else:
            gp120_str = f"{p['res2']}{p['res2_num']}"
            cd4_str = f"{p['res1']}{p['res1_num']}"
        print(f"{gp120_str:>12} {cd4_str:>12} {p['vdw']:>8.3f} {p['elec']:>8.3f} {p['polar']:>8.3f} {p['total']:>8.3f}")

    relay_total = sum(p['total'] for p in relay_pairs)
    interface_total = sum(p['total'] for p in interface)
    print(f"\nRelay pocket pairwise sum: {relay_total:.3f} kcal/mol")
    print(f"Total interface pairwise sum: {interface_total:.3f} kcal/mol")
    if interface_total != 0:
        print(f"Relay fraction: {abs(relay_total/interface_total)*100:.1f}%")

    results = {
        'method': 'idecomp=3 pairwise DELTAS',
        'total_pairs_parsed': len(pairs),
        'interface_pairs': len(interface),
        'relay_pairs_count': len(relay_pairs),
        'interface_total_dG': round(interface_total, 4),
        'relay_total_dG': round(relay_total, 4),
        'top20_interface': [],
        'relay_pairs': [],
    }

    for i, p in enumerate(interface_sorted[:20]):
        if p['res1_num'] <= GP120_MAX:
            gp120_res, gp120_num = p['res1'], p['res1_num']
            cd4_res, cd4_num = p['res2'], p['res2_num']
        else:
            gp120_res, gp120_num = p['res2'], p['res2_num']
            cd4_res, cd4_num = p['res1'], p['res1_num']
        results['top20_interface'].append({
            'rank': i+1,
            'gp120_res': gp120_res, 'gp120_num': gp120_num,
            'cd4_res': cd4_res, 'cd4_num': cd4_num,
            'dG_total': p['total'], 'dG_vdw': p['vdw'],
            'dG_elec': p['elec'], 'dG_polar': p['polar'],
            'is_relay': gp120_num in RELAY_RESIDUES,
        })

    for p in relay_sorted:
        if p['res1_num'] <= GP120_MAX:
            gp120_res, gp120_num = p['res1'], p['res1_num']
            cd4_res, cd4_num = p['res2'], p['res2_num']
        else:
            gp120_res, gp120_num = p['res2'], p['res2_num']
            cd4_res, cd4_num = p['res1'], p['res1_num']
        results['relay_pairs'].append({
            'gp120_res': gp120_res, 'gp120_num': gp120_num,
            'cd4_res': cd4_res, 'cd4_num': cd4_num,
            'dG_total': p['total'], 'dG_vdw': p['vdw'],
            'dG_elec': p['elec'], 'dG_polar': p['polar'],
        })

    out_json = DECOMP_FILE.replace('.dat', '_results.json')
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_json}")


if __name__ == '__main__':
    main()
