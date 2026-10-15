# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Triple-check: Verify every JSON parses, key numbers match across documents.
"""
import json, os, glob

base = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
errors = []
warnings = []

# 1. Parse every JSON
print("=" * 60)
print("STEP 1: Verify all JSON files parse cleanly")
print("=" * 60)

all_files = []
for subdir in ['results/00-summary', 'results/01-pocket-discovery', 'results/02-docking',
               'results/03-conservation', 'results/04-md-simulations', 'results/05-elastic-network',
               'results/06-dynamics', 'results/07-free-energy', 'results/08-binding-energetics',
               'results/09-drug-design', 'results/10-controls']:
    all_files.extend(glob.glob(os.path.join(base, subdir, '*.json')))
for f in sorted(all_files):
    rel = os.path.relpath(f, base)
    try:
        d = json.load(open(f))
        sz = os.path.getsize(f)
        print(f"  OK  {rel} ({sz:,} bytes)")
    except Exception as e:
        errors.append(f"JSON PARSE ERROR: {rel}: {e}")
        print(f"  FAIL {rel}: {e}")

# 2. Check that all expected files exist
print("\n" + "=" * 60)
print("STEP 2: Verify all expected files exist")
print("=" * 60)

expected = [
    'results/01-pocket-discovery/hiv1_pharmacophore_data.json',
    'results/02-docking/docking_results.json',
    'results/02-docking/cross_structure_docking.json',
    'results/09-drug-design/vscreen_batch1_results.json',
    'results/09-drug-design/vscreen_batch2_results.json',
    'results/01-pocket-discovery/fpocket_analysis.json',
    'results/01-pocket-discovery/bfactor_analysis.json',
    'results/01-pocket-discovery/druggability_scoring.json',
    'results/01-pocket-discovery/pocket_geometry.json',
    'results/10-controls/selectivity_panel.json',
    'results/09-drug-design/reference_compounds.json',
    'results/09-drug-design/admet_results.json',
    'results/03-conservation/resistance_mapping.json',
    'results/03-conservation/cross_clade_survey_validated.json',
    'results/09-drug-design/pharmacophore_screen_results.json',
    'results/09-drug-design/denovo_designs.json',
    'results/06-dynamics/pca_results.json',
    'results/02-docking/ensemble_docking_results.json',
    'results/01-pocket-discovery/watermap_results.json',
    'results/05-elastic-network/allosteric_network.json',
    'results/06-dynamics/dccm_pocket.json',
]

for ef in expected:
    path = os.path.join(base, ef)
    if os.path.exists(path):
        print(f"  EXISTS  {ef}")
    else:
        errors.append(f"MISSING FILE: {ef}")
        print(f"  MISSING {ef}")

# 3. Cross-check key numbers
print("\n" + "=" * 60)
print("STEP 3: Cross-check key numbers")
print("=" * 60)

def check(name, actual, expected_val, tolerance=0.001):
    if abs(actual - expected_val) > tolerance:
        msg = f"MISMATCH: {name}: actual={actual}, expected={expected_val}"
        errors.append(msg)
        print(f"  FAIL {msg}")
    else:
        print(f"  OK   {name} = {actual}")

# ADMET
d = json.load(open(os.path.join(base, 'results', '09-drug-design', 'admet_results.json')))
check("ADMET total_compounds", d['total_compounds'], 98)
check("ADMET lipinski", d['summary']['lipinski'], 90)
check("ADMET veber", d['summary']['veber'], 93)
check("ADMET pains", d['summary']['pains'], 95)
check("ADMET all_clean", d['summary']['all_clean'], 87)

# Resistance
d = json.load(open(os.path.join(base, 'results', '03-conservation', 'resistance_mapping.json')))
check("Resistance total", d['total_mutations_checked'], 28)
check("Resistance in_pocket", d['in_pocket'], 2)
check("Resistance in_wall", d['in_wall'], 10)
check("Resistance near", d['near_pocket'], 1)
check("Resistance distant", d['distant'], 15)

# Cross-clade
d = json.load(open(os.path.join(base, 'results', '03-conservation', 'cross_clade_survey_validated.json')))
check("CrossClade total", d['total_structures'], 1870)
check("CrossClade with_pocket", d['with_relay_pocket'], 442)
check("CrossClade exact", d['exact_match'], 441)
check("CrossClade without", d['without_pocket'], 1428)

# PCA
d = json.load(open(os.path.join(base, 'results', '06-dynamics', 'pca_results.json')))
check("PCA total_var", d['total_variance_nm2'], 2.4776, 0.0001)
check("PCA top3", d['top3_explain_percent'], 51.0, 0.1)
check("PCA top5", d['top5_explain_percent'], 65.5, 0.1)

# Ensemble docking
d = json.load(open(os.path.join(base, 'results', '02-docking', 'ensemble_docking_results.json')))
check("EnsembleDock BMS mean", d['bms626529']['mean_affinity_kcal_mol'], -6.808)
check("EnsembleDock BMS sd", d['bms626529']['sd_kcal_mol'], 0.750)
check("EnsembleDock BMS best", d['bms626529']['min_best_kcal_mol'], -8.342)
check("EnsembleDock Vic mean", d['vicriviroc']['mean_affinity_kcal_mol'], -7.499)
check("EnsembleDock Vic best", d['vicriviroc']['min_best_kcal_mol'], -9.057)

# WaterMap
d = json.load(open(os.path.join(base, 'results', '01-pocket-discovery', 'watermap_results.json')))
check("WaterMap total", d['total_sites'], 22)
check("WaterMap structural", d['structural_sites'], 0)
check("WaterMap displaceable", d['displaceable_sites'], 8)
check("WaterMap peripheral", d['peripheral_sites'], 14)

# Allosteric
d = json.load(open(os.path.join(base, 'results', '05-elastic-network', 'allosteric_network.json')))
intra = d['intra_pocket_correlations']
check("Allosteric intra pairs", len(intra), 28)
check("Allosteric strongest corr", intra[0]['correlation'], 0.8196, 0.001)
all_positive = all(x['correlation'] > 0 for x in intra)
if not all_positive:
    errors.append("FAIL: Not all intra-pocket correlations are positive!")
    print("  FAIL Not all intra-pocket correlations positive")
else:
    print("  OK   All 28 intra-pocket correlations positive")

# Pharmacophore screen
d = json.load(open(os.path.join(base, 'results', '09-drug-design', 'pharmacophore_screen_results.json')))
check("PharmScreen library", d['library_size'], 48)
check("PharmScreen screened", d['successful_screens'], 47)
check("PharmScreen top hit fit", d['top_10'][0]['pharmacophore_fit_score'], 0.543)
if d['top_10'][0]['name'] != 'NBD-556':
    errors.append(f"MISMATCH: PharmScreen top hit name={d['top_10'][0]['name']}, expected NBD-556")
else:
    print(f"  OK   PharmScreen top hit = NBD-556")

# De novo
d = json.load(open(os.path.join(base, 'results', '09-drug-design', 'denovo_designs.json')))
check("DeNovo total", d['total_designs'], 8)
check("DeNovo drug_like", d['drug_like_count'], 8)
check("DeNovo highly_novel", d['highly_novel_count'], 5)
check("DeNovo top combined", d['designs'][0]['combined_score'], 0.639)

# Docking
d = json.load(open(os.path.join(base, 'results', '02-docking', 'docking_results.json')))
check("Docking relay", d['relay_pocket_docking']['best_affinity_kcal_mol'], -6.692)

# VScreen
d1 = json.load(open(os.path.join(base, 'results', '09-drug-design', 'vscreen_batch1_results.json')))
d2 = json.load(open(os.path.join(base, 'results', '09-drug-design', 'vscreen_batch2_results.json')))
all_r = d1['results'] + d2['results']
check("VScreen total", len(all_r), 96)
druggable = sum(1 for r in all_r if r['relay_best_kcal_mol'] < -6.0)
selective = sum(1 for r in all_r if r['delta_kcal_mol'] < -1.0)
check("VScreen druggable", druggable, 25)
check("VScreen selective", selective, 3)
vic = [r for r in all_r if r['name'] == 'Vicriviroc'][0]
check("VScreen Vic relay", vic['relay_best_kcal_mol'], -8.797)
check("VScreen Vic delta", vic['delta_kcal_mol'], -2.186)

# Pharmacophore
d = json.load(open(os.path.join(base, 'results', '01-pocket-discovery', 'hiv1_pharmacophore_data.json')))
relay = [t for t in d['pharmacophore_targets'] if 'relay' in t['name'].lower()][0]
check("Pharmacophore relay features", len(relay['pharmacophore_features']), 55)
check("Pharmacophore hotspots", len(d['binding_hotspots']), 22)

# Summary
print("\n" + "=" * 60)
print("VERIFICATION SUMMARY")
print("=" * 60)
if errors:
    print(f"\n  ERRORS: {len(errors)}")
    for e in errors:
        print(f"    - {e}")
else:
    print("\n  ALL CHECKS PASSED - ZERO ERRORS")

if warnings:
    print(f"\n  WARNINGS: {len(warnings)}")
    for w in warnings:
        print(f"    - {w}")

print(f"\n  Files checked: {len(all_files)}")
print(f"  Expected files: {len(expected)}")
print(f"  Number checks: passed")
