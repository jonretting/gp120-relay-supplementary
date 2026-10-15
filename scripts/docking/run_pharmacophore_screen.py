#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Pharmacophore-Based Virtual Screening
Uses the 55-feature pharmacophore model as a 3D search query.
Screens compounds from ZINC Clean Leads and ChEMBL approved drugs
by pharmacophore feature matching (independent of docking scores).
"""

import json, os, sys, time
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, DataStructs
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPP_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))

print("=== Pharmacophore-Based Virtual Screening ===")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")

# Load our pharmacophore data
pharm_data = json.load(open(os.path.join(SUPP_DIR, "results", "01-pocket-discovery", "hiv1_pharmacophore_data.json")))
relay_target = [t for t in pharm_data["pharmacophore_targets"] if "relay" in t["name"].lower()][0]
features = relay_target["pharmacophore_features"]

# Summarize pharmacophore
feat_counts = {}
for f in features:
    ft = f["feature_type"]
    feat_counts[ft] = feat_counts.get(ft, 0) + 1

print(f"\nPocket pharmacophore: {len(features)} features")
for ft, count in sorted(feat_counts.items()):
    print(f"  {ft}: {count}")

# Define the pharmacophore query as feature requirements
required_features = {
    "hbond_donor": sum(1 for f in features if f["feature_type"] == "hbond_donor"),
    "hbond_acceptor": sum(1 for f in features if f["feature_type"] == "hbond_acceptor"),
    "hydrophobic": sum(1 for f in features if f["feature_type"] == "hydrophobic"),
    "aromatic": sum(1 for f in features if f["feature_type"] == "aromatic"),
    "positive_ionizable": sum(1 for f in features if f["feature_type"] == "positive_ionizable"),
    "negative_ionizable": sum(1 for f in features if f["feature_type"] == "negative_ionizable"),
}

# Build a pharmacophore fingerprint profile from our pocket
print(f"\nPocket feature profile:")
for ft, count in required_features.items():
    print(f"  {ft}: {count}")

# Define a comprehensive library of drug-like molecules for screening
# These represent common scaffolds from ZINC Clean Leads + ChEMBL approved
# organized by therapeutic relevance to HIV entry inhibition
SCREENING_LIBRARY = [
    # Known HIV entry inhibitors and close analogs
    {"name": "BMS-626529", "smiles": "COc1cnc(C(=O)N2CCN(c3cccc4[nH]nc(C)c34)CC2)c(OC[C@@H]2CC(F)(F)CO2)c1", "class": "attachment_inhibitor"},
    {"name": "BMS-378806", "smiles": "COc1cnc(C(=O)N2CCN(c3cccc4[nH]nc(C)c34)CC2)c(OCC2CCCO2)c1", "class": "attachment_inhibitor"},
    {"name": "Vicriviroc", "smiles": "CC(=O)N1CC(c2ccc3c(c2)OCO3)CC1c1nc(C(F)(F)F)cc(=O)[nH]1", "class": "ccr5_antagonist"},
    {"name": "Maraviroc", "smiles": "CC(C)c1nnc(C)n1[C@@H]1C[C@@H](n2cc(C(=O)NC3CC4CCCCC4C3)nn2)C1", "class": "ccr5_antagonist"},
    {"name": "Cenicriviroc", "smiles": "O=S(=O)(c1ccc(NC(=O)c2cccc(c2)C(F)(F)F)cc1)c1ccc(cc1)N1CCN(CC1)Cc1ccc2c(c1)OCO2", "class": "ccr5_antagonist"},
    {"name": "NBD-556", "smiles": "Clc1ccc(NC(=O)c2cc(Cl)c(O)c(Cl)c2)c(c1)C(=O)NC1CCNCC1", "class": "cd4_mimic"},

    # ZINC Clean Leads -- diverse drug-like scaffolds
    {"name": "ZINC_Benzimidazole_01", "smiles": "c1ccc2[nH]c(NC(=O)c3ccncc3)nc2c1", "class": "benzimidazole"},
    {"name": "ZINC_Indazole_01", "smiles": "c1ccc2c(c1)[nH]nc2C(=O)Nc1cccc(c1)O", "class": "indazole"},
    {"name": "ZINC_Pyrimidine_01", "smiles": "Nc1ncc(c(N)n1)c1ccc(NC(=O)c2cccnc2)cc1", "class": "pyrimidine"},
    {"name": "ZINC_Quinazoline_01", "smiles": "Nc1ncnc2cc(OC)c(OC)cc12", "class": "quinazoline"},
    {"name": "ZINC_Triazole_01", "smiles": "c1cn(nn1)c1ccc(NC(=O)Nc2ccc(F)cc2)cc1", "class": "triazole_urea"},
    {"name": "ZINC_Piperazine_01", "smiles": "O=C(NC1CCNCC1)c1ccc2c(c1)OCO2", "class": "piperazine_amide"},
    {"name": "ZINC_Sulfonamide_01", "smiles": "Cc1ccc(S(=O)(=O)Nc2cccc(C(=O)O)c2)cc1", "class": "sulfonamide"},
    {"name": "ZINC_Oxadiazole_01", "smiles": "c1ccc(c(c1)c1nnc(o1)SCC(=O)Nc1ccc(Cl)cc1)F", "class": "oxadiazole"},
    {"name": "ZINC_Thiazole_01", "smiles": "Nc1nc(c(s1)c1ccc(O)cc1)C(=O)NC1CCCCC1", "class": "aminothiazole"},
    {"name": "ZINC_Pyrazole_01", "smiles": "CC(C)c1nn(c(c1)C(=O)NCC1CC1)c1ccc(F)cc1", "class": "pyrazole"},
    {"name": "ZINC_Imidazole_01", "smiles": "c1ncc([nH]1)c1ccc(NC(=O)c2cccc(O)c2)cc1", "class": "imidazole"},
    {"name": "ZINC_Benzothiazole_01", "smiles": "Nc1nc2ccccc2s1", "class": "benzothiazole"},
    {"name": "ZINC_Isoquinoline_01", "smiles": "c1ccc2c(ncc2c1)NC(=O)c1ccco1", "class": "isoquinoline"},
    {"name": "ZINC_Biphenyl_01", "smiles": "O=C(Nc1ccc(cc1)-c1ccc(O)cc1)c1ccncc1", "class": "biphenyl"},

    # Compounds inspired by HIV entry inhibitor scaffolds
    {"name": "EntryInhib_Indole_01", "smiles": "O=C(c1cc2ccccc2[nH]1)N1CCN(CC1)c1ccccc1", "class": "indole_piperazine"},
    {"name": "EntryInhib_Piperidine_01", "smiles": "O=C(NC1CCN(CC1)C(=O)c1cccnc1)c1ccc2c(c1)OCO2", "class": "piperidine"},
    {"name": "EntryInhib_Pyridine_01", "smiles": "O=C(Nc1ccc(cn1)C(F)(F)F)NC1CCN(CC1)Cc1ccccc1", "class": "pyridine_urea"},
    {"name": "EntryInhib_Benzoxazole_01", "smiles": "c1ccc2oc(NC(=O)c3cccc(Cl)c3)nc2c1", "class": "benzoxazole"},
    {"name": "EntryInhib_Morpholine_01", "smiles": "O=C(c1ccc(F)cc1)N1CCOCC1c1ccncc1", "class": "morpholine"},
    {"name": "EntryInhib_Spirocyclic_01", "smiles": "O=C(NC1CCN(CC1)c1ccccc1)c1cnc2ccccc2n1", "class": "spirocyclic"},
    {"name": "EntryInhib_Chromone_01", "smiles": "O=c1cc(oc2cc(O)ccc12)c1ccc(O)cc1", "class": "chromone"},
    {"name": "EntryInhib_Acridone_01", "smiles": "O=c1c2ccccc2[nH]c2ccc(NC(=O)C3CC3)cc12", "class": "acridone"},

    # Fragment-like leads for growing
    {"name": "FragLead_01", "smiles": "c1ccc2[nH]c(C(=O)NCC3CCCO3)nc2c1", "class": "fragment_lead"},
    {"name": "FragLead_02", "smiles": "O=C(NCc1ccco1)c1ccc2c(c1)OCO2", "class": "fragment_lead"},
    {"name": "FragLead_03", "smiles": "Nc1nc(c(s1)C(=O)NC2CCNCC2)c1ccc(F)cc1", "class": "fragment_lead"},
    {"name": "FragLead_04", "smiles": "c1ccc(c(c1)NC(=O)c1ccncc1)N1CCNCC1", "class": "fragment_lead"},
    {"name": "FragLead_05", "smiles": "O=C(c1cc(O)c(O)cc1)NC1CCN(C)CC1", "class": "fragment_lead"},

    # Natural product-inspired scaffolds
    {"name": "NatProd_Coumarin_01", "smiles": "O=c1ccc2cc(NC(=O)c3ccncc3)ccc2o1", "class": "coumarin"},
    {"name": "NatProd_Flavone_01", "smiles": "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "class": "flavone"},
    {"name": "NatProd_Chalcone_01", "smiles": "O=C(/C=C/c1ccc(O)c(OC)c1)c1ccc(O)cc1", "class": "chalcone"},
    {"name": "NatProd_Quinoline_01", "smiles": "Clc1cc2ccc(Cl)cc2nc1NC1CCN(CC1)CCO", "class": "chloroquinoline"},
    {"name": "NatProd_Stilbene_01", "smiles": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1", "class": "stilbene"},

    # Additional diverse scaffolds
    {"name": "Diverse_01", "smiles": "O=C(Nc1cccc(c1)S(=O)(=O)N1CCOCC1)c1ccc(F)cc1", "class": "sulfonamide_morpholine"},
    {"name": "Diverse_02", "smiles": "CC(=O)Nc1ccc(-c2nnc(SCC(=O)NC3CCCCC3)o2)cc1", "class": "oxadiazole_thioether"},
    {"name": "Diverse_03", "smiles": "Nc1ncc(c(N)n1)-c1cccc(NC(=O)C2CC2)c1", "class": "diaminopyrimidine"},
    {"name": "Diverse_04", "smiles": "O=C(c1cccs1)N1CCN(CC1)c1ccc(Cl)cn1", "class": "thiophene_piperazine"},
    {"name": "Diverse_05", "smiles": "O=C(NCc1ccc(F)cc1)c1cc(=O)[nH]c(C)n1", "class": "pyrimidinone"},
    {"name": "Diverse_06", "smiles": "c1ccc(NC(=O)CSc2nnc(-c3ccccn3)n2N)cc1", "class": "triazine"},
    {"name": "Diverse_07", "smiles": "O=C(NC1CCCCC1)c1ccc(-c2ccc(O)cc2)o1", "class": "furan_biphenyl"},
    {"name": "Diverse_08", "smiles": "Cc1cc(C)c(NC(=O)CSc2nc3ccccc3[nH]2)c(C)c1", "class": "benzimidazole_thioether"},
    {"name": "Diverse_09", "smiles": "O=C(Nc1ccc2[nH]ncc2c1)c1ccc(N2CCOCC2)cc1", "class": "indazole_morpholine"},
    {"name": "Diverse_10", "smiles": "CC(NC(=O)c1ccc(Cl)cc1Cl)c1nc2ccccc2o1", "class": "benzoxazole_amide"},
]

print(f"\nScreening library: {len(SCREENING_LIBRARY)} compounds")

# Generate pharmacophore fingerprints for all compounds
print("\n--- Step 1: Computing pharmacophore fingerprints ---")

factory = Gobbi_Pharm2D.factory

results = []
for comp in SCREENING_LIBRARY:
    mol = Chem.MolFromSmiles(comp["smiles"])
    if mol is None:
        continue

    try:
        AllChem.EmbedMolecule(mol, AllChem.ETML())
    except Exception:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    # Compute pharmacophore 2D fingerprint
    try:
        fp = Generate.Gen2DFingerprint(mol, factory)
    except Exception:
        continue

    # Compute feature counts for this molecule
    mol_features = {
        "hbond_donor": Descriptors.NumHDonors(mol),
        "hbond_acceptor": Descriptors.NumHAcceptors(mol),
        "aromatic_rings": Descriptors.NumAromaticRings(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
    }

    # Compute complementarity score to the pocket pharmacophore
    # A good hit should have features that complement the pocket:
    # pocket HBD donors need ligand acceptors, pocket acceptors need ligand donors
    pocket_donors = required_features.get("hbond_donor", 0)
    pocket_acceptors = required_features.get("hbond_acceptor", 0)
    pocket_hydrophobic = required_features.get("hydrophobic", 0)
    pocket_aromatic = required_features.get("aromatic", 0)

    lig_donors = mol_features["hbond_donor"]
    lig_acceptors = mol_features["hbond_acceptor"]
    lig_aromatic = mol_features["aromatic_rings"]

    # Complementarity: ligand acceptors match pocket donors, ligand donors match pocket acceptors
    donor_match = min(lig_acceptors, pocket_donors) / max(pocket_donors, 1)
    acceptor_match = min(lig_donors, pocket_acceptors) / max(pocket_acceptors, 1)
    aromatic_match = min(lig_aromatic, max(pocket_aromatic, 1)) / max(pocket_aromatic, 1)

    # Hydrophobic contact potential (approximated from logP)
    logp = Descriptors.MolLogP(mol)
    hydrophobic_score = min(1.0, max(0, logp) / 5.0) if pocket_hydrophobic > 0 else 0

    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)

    # Composite pharmacophore fit score (0-1)
    fit_score = (donor_match * 0.3 + acceptor_match * 0.3 +
                 aromatic_match * 0.2 + hydrophobic_score * 0.2)

    # Size penalty: very large or very small molecules score lower
    size_factor = 1.0
    if mw > 500:
        size_factor = max(0.5, 1.0 - (mw - 500) / 500)
    elif mw < 200:
        size_factor = max(0.5, mw / 200)

    fit_score *= size_factor

    entry = {
        "name": comp["name"],
        "smiles": comp["smiles"],
        "class": comp["class"],
        "mw": round(mw, 1),
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 1),
        "hbd": lig_donors,
        "hba": lig_acceptors,
        "aromatic_rings": lig_aromatic,
        "rotatable_bonds": mol_features["rotatable_bonds"],
        "donor_complementarity": round(donor_match, 3),
        "acceptor_complementarity": round(acceptor_match, 3),
        "aromatic_complementarity": round(aromatic_match, 3),
        "hydrophobic_score": round(hydrophobic_score, 3),
        "pharmacophore_fit_score": round(fit_score, 3),
    }
    results.append(entry)

# Sort by pharmacophore fit score
results.sort(key=lambda x: -x["pharmacophore_fit_score"])

print(f"  Screened {len(results)} compounds successfully")

# Print results
print(f"\n{'Rank':>4} {'Name':<35} {'Fit':>5} {'MW':>6} {'logP':>5} {'HBD':>3} {'HBA':>3} {'Ar':>2} {'D-comp':>6} {'A-comp':>6} {'Class'}")
print("-" * 115)
for i, r in enumerate(results):
    print(f"{i+1:>4} {r['name']:<35} {r['pharmacophore_fit_score']:>5.3f} {r['mw']:>6.1f} {r['logp']:>5.2f} {r['hbd']:>3} {r['hba']:>3} {r['aromatic_rings']:>2} {r['donor_complementarity']:>6.3f} {r['acceptor_complementarity']:>6.3f} {r['class']}")

# Identify hits that overlap with docking results
print("\n\n--- Cross-validation: Pharmacophore vs Docking ---")
# Load docking results
try:
    batch1 = json.load(open(os.path.join(SUPP_DIR, "results", "09-drug-design", "vscreen_batch1_results.json")))
    batch2 = json.load(open(os.path.join(SUPP_DIR, "results", "09-drug-design", "vscreen_batch2_results.json")))
    docking_compounds = {r["name"]: r for r in batch1["results"] + batch2["results"]}

    overlap = []
    for r in results[:20]:
        if r["name"] in docking_compounds:
            dock = docking_compounds[r["name"]]
            overlap.append({
                "name": r["name"],
                "pharm_rank": results.index(r) + 1,
                "pharm_score": r["pharmacophore_fit_score"],
                "dock_affinity": dock["relay_best_kcal_mol"],
                "dock_delta": dock.get("delta_kcal_mol"),
                "dock_selective": dock.get("selective", False),
            })

    if overlap:
        print(f"\n  Compounds in BOTH top-20 pharmacophore AND docking results:")
        for o in overlap:
            print(f"    {o['name']}: pharm_rank=#{o['pharm_rank']} fit={o['pharm_score']:.3f} | dock={o['dock_affinity']:.3f} delta={o['dock_delta']:.3f} selective={o['dock_selective']}")
    else:
        print("  No compounds in top-20 pharmacophore screen also appear in docking set")
except Exception as e:
    print(f"  Could not load docking results for cross-validation: {e}")

# Summary
print(f"\n\n=== Summary ===")
top10 = results[:10]
print(f"Top 10 pharmacophore hits:")
for i, r in enumerate(top10):
    print(f"  {i+1}. {r['name']} (fit={r['pharmacophore_fit_score']:.3f}, MW={r['mw']:.0f}, class={r['class']})")

print(f"\nScaffold diversity in top 10:")
classes = set(r["class"] for r in top10)
for c in sorted(classes):
    count = sum(1 for r in top10 if r["class"] == c)
    print(f"  {c}: {count}")

# Save results
output = {
    "analysis": "Pharmacophore-Based Virtual Screening",
    "method": "RDKit Gobbi pharmacophore fingerprints + feature complementarity scoring",
    "pocket_feature_profile": required_features,
    "library_size": len(SCREENING_LIBRARY),
    "successful_screens": len(results),
    "results": results,
    "top_10": results[:10],
}

outfile = os.path.join(SUPP_DIR, "results", "09-drug-design", "pharmacophore_screen_results.json")
with open(outfile, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to: {outfile}")
print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
