#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
ADMET Prediction for Top Virtual Screening Hits
Computes drug-likeness, Lipinski, PAINS, and physicochemical properties
using RDKit descriptors (no external web service needed).
"""

import json, os, sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, FilterCatalog, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalogParams

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
OUT_DIR = os.path.join(SUPP_DIR, "results", "09-drug-design")

# Load virtual screening results
batch1 = json.load(open(os.path.join(SUPP_DIR, "results", "09-drug-design", "vscreen_batch1_results.json")))
batch2 = json.load(open(os.path.join(SUPP_DIR, "results", "09-drug-design", "vscreen_batch2_results.json")))
ref_compounds = json.load(open(os.path.join(SUPP_DIR, "results", "09-drug-design", "reference_compounds.json")))

all_compounds = batch1["results"] + batch2["results"]

# Also add reference compounds
for key, ref in ref_compounds.items():
    if "smiles" in ref:
        all_compounds.append({
            "name": ref["name"],
            "smiles": ref["smiles"],
            "tier": "REF",
            "class": ref.get("mechanism", "reference"),
            "relay_best_kcal_mol": None,
            "delta_kcal_mol": None,
        })

# Setup PAINS filter
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains_catalog = FilterCatalog.FilterCatalog(params)

print("=== ADMET Prediction for HIV-1 Relay Pocket Virtual Screening Hits ===")
print(f"Total compounds: {len(all_compounds)}")
print()

results = []

for comp in all_compounds:
    name = comp["name"]
    smiles = comp.get("smiles", "")
    if not smiles:
        continue

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  SKIP {name}: invalid SMILES")
        continue

    # Physicochemical properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Lipinski.NumRotatableBonds(mol)
    rings = Lipinski.RingCount(mol)
    aromatic_rings = Descriptors.NumAromaticRings(mol)
    heavy_atoms = Lipinski.HeavyAtomCount(mol)
    fraction_csp3 = Descriptors.FractionCSP3(mol)

    # Lipinski Rule of Five
    lipinski_violations = 0
    if mw > 500: lipinski_violations += 1
    if logp > 5: lipinski_violations += 1
    if hbd > 5: lipinski_violations += 1
    if hba > 10: lipinski_violations += 1
    lipinski_pass = lipinski_violations <= 1

    # Veber rules (oral bioavailability)
    veber_pass = tpsa <= 140 and rotatable <= 10

    # Lead-likeness (fragment-like)
    leadlike = 250 <= mw <= 350 and logp <= 3.5 and rotatable <= 7

    # Drug-likeness (Ghose)
    ghose = 160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 40 <= heavy_atoms * 12 / mw * 100 <= 130

    # PAINS filter
    pains_match = pains_catalog.HasMatch(mol)
    pains_alerts = []
    if pains_match:
        entry = pains_catalog.GetFirstMatch(mol)
        if entry:
            pains_alerts.append(entry.GetDescription())

    # Synthetic accessibility (Ertl SA score approximation via heavy atoms & rings)
    # True SA score requires sascorer module; approximate here
    sa_approx = min(10, max(1, heavy_atoms / 5 + rings * 0.5))

    # BBB permeability estimate (simple: MW < 400, TPSA < 90, logP 1-3)
    bbb_likely = mw < 400 and tpsa < 90 and 1 < logp < 3

    # hERG risk estimate (logP > 3.7 and basic nitrogen -> higher risk)
    has_basic_n = "[NH]" in smiles or "[nH]" in smiles or "N(" in smiles
    herg_risk = "moderate" if (logp > 3.7 and has_basic_n) else "low"

    entry = {
        "name": name,
        "smiles": smiles,
        "tier": comp.get("tier", "?"),
        "class": comp.get("class", "unknown"),
        "relay_affinity": comp.get("relay_best_kcal_mol"),
        "delta": comp.get("delta_kcal_mol"),
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotatable,
        "rings": rings,
        "aromatic_rings": aromatic_rings,
        "heavy_atoms": heavy_atoms,
        "fraction_csp3": round(fraction_csp3, 3),
        "lipinski_violations": lipinski_violations,
        "lipinski_pass": lipinski_pass,
        "veber_oral": veber_pass,
        "leadlike": leadlike,
        "pains_clean": not pains_match,
        "pains_alerts": pains_alerts,
        "bbb_permeable": bbb_likely,
        "herg_risk": herg_risk,
        "sa_approx": round(sa_approx, 1),
    }
    results.append(entry)

# Sort by tier then affinity
tier_order = {"A": 0, "B": 1, "C": 2, "D": 3, "REF": 4}
results.sort(key=lambda x: (tier_order.get(x["tier"], 5), x.get("relay_affinity") or 0))

# Print summary table
print(f"{'Name':<30} {'MW':>6} {'logP':>5} {'HBD':>3} {'HBA':>3} {'TPSA':>6} {'Rot':>3} {'Lip':>3} {'Veb':>3} {'PAINS':>5} {'hERG':>6} {'Tier':>4}")
print("-" * 105)

pass_count = {"lipinski": 0, "veber": 0, "pains": 0, "all_clean": 0}
for r in results:
    lip = "PASS" if r["lipinski_pass"] else f"FAIL({r['lipinski_violations']})"
    veb = "YES" if r["veber_oral"] else "NO"
    pains = "CLEAN" if r["pains_clean"] else "ALERT"
    print(f"{r['name']:<30} {r['mw']:>6.1f} {r['logp']:>5.2f} {r['hbd']:>3} {r['hba']:>3} {r['tpsa']:>6.1f} {r['rotatable_bonds']:>3} {lip:>7} {veb:>3} {pains:>5} {r['herg_risk']:>6} {r['tier']:>4}")

    if r["lipinski_pass"]: pass_count["lipinski"] += 1
    if r["veber_oral"]: pass_count["veber"] += 1
    if r["pains_clean"]: pass_count["pains"] += 1
    if r["lipinski_pass"] and r["veber_oral"] and r["pains_clean"]:
        pass_count["all_clean"] += 1

print()
print(f"=== Summary ===")
print(f"Total compounds analyzed: {len(results)}")
print(f"Lipinski compliant:       {pass_count['lipinski']} ({pass_count['lipinski']/len(results)*100:.0f}%)")
print(f"Veber oral bioavail:      {pass_count['veber']} ({pass_count['veber']/len(results)*100:.0f}%)")
print(f"PAINS clean:              {pass_count['pains']} ({pass_count['pains']/len(results)*100:.0f}%)")
print(f"All filters passed:       {pass_count['all_clean']} ({pass_count['all_clean']/len(results)*100:.0f}%)")

# Top 10 drug-like selective hits
print()
print("=== Top Drug-Like Selective Hits (Lipinski + Veber + PAINS clean + delta < -1) ===")
top_hits = [r for r in results if r["lipinski_pass"] and r["veber_oral"] and r["pains_clean"]
            and r.get("delta") is not None and r["delta"] < -1.0 and r["tier"] in ("A", "B")]
top_hits.sort(key=lambda x: x.get("delta", 0))
for i, r in enumerate(top_hits[:15]):
    print(f"  {i+1}. {r['name']} | {r['relay_affinity']:.3f} kcal/mol | delta={r['delta']:.3f} | MW={r['mw']:.0f} | logP={r['logp']:.2f} | TPSA={r['tpsa']:.0f}")

# Save full results
output = {
    "analysis": "ADMET / Drug-Likeness Prediction",
    "method": "RDKit descriptors + Lipinski/Veber/PAINS filters",
    "total_compounds": len(results),
    "summary": pass_count,
    "compounds": results,
}

outfile = os.path.join(OUT_DIR, "admet_results.json")
with open(outfile, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to: {outfile}")
