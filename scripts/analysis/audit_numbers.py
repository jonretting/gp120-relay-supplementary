# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
import json, os

base = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))

def load(name):
    return json.load(open(os.path.join(base, name)))

# ADMET
d = load('results/09-drug-design/admet_results.json')
s = d['summary']
print(f"ADMET: total={d['total_compounds']}, lipinski={s['lipinski']}, veber={s['veber']}, pains={s['pains']}, all_clean={s['all_clean']}")
print(f"  lipinski_pct={round(s['lipinski']/d['total_compounds']*100)}, veber_pct={round(s['veber']/d['total_compounds']*100)}, pains_pct={round(s['pains']/d['total_compounds']*100)}")

# Resistance
d = load('results/03-conservation/resistance_mapping.json')
print(f"Resistance: total_checked={d['total_mutations_checked']}, in_pocket={d['in_pocket']}, in_wall={d['in_wall']}, near={d['near_pocket']}, distant={d['distant']}")
print(f"  key_finding={d['key_finding']}")
print(f"  polymorphisms={d['pocket_polymorphisms']}")

# Cross-clade
d = load('results/03-conservation/cross_clade_survey_validated.json')
print(f"CrossClade: total={d['total_structures']}, with_pocket={d['with_relay_pocket']}, exact={d['exact_match']}, without={d['without_pocket']}, errors={d['errors']}")

# Pharmacophore screen
d = load('results/09-drug-design/pharmacophore_screen_results.json')
print(f"PharmScreen: library={d['library_size']}, screened={d['successful_screens']}")
for i, r in enumerate(d['top_10'][:5]):
    print(f"  #{i+1}: {r['name']} fit={r['pharmacophore_fit_score']:.3f} MW={r['mw']} class={r['class']}")

# De novo
d = load('results/09-drug-design/denovo_designs.json')
print(f"DeNovo: total={d['total_designs']}, drug_like={d['drug_like_count']}, highly_novel={d['highly_novel_count']}")
for r in d['designs'][:3]:
    print(f"  {r['name']}: comb={r['combined_score']:.3f} fit={r['pharmacophore_fit']:.3f} nov={r['novelty_score']:.3f} MW={r['mw']} tani={r['max_tanimoto_to_known']:.3f}")

# Ensemble docking
d = load('results/02-docking/ensemble_docking_results.json')
bms = d['bms626529']
vic = d['vicriviroc']
print(f"EnsembleDock: clusters={d['clusters']}")
print(f"  BMS: mean={bms['mean_affinity_kcal_mol']:.3f} sd={bms['sd_kcal_mol']:.3f} best={bms['min_best_kcal_mol']:.3f} worst={bms['max_worst_kcal_mol']:.3f} orig={bms['original_single_frame_kcal_mol']:.3f}")
print(f"  Vic: mean={vic['mean_affinity_kcal_mol']:.3f} sd={vic['sd_kcal_mol']:.3f} best={vic['min_best_kcal_mol']:.3f} worst={vic['max_worst_kcal_mol']:.3f} orig={vic['original_single_frame_kcal_mol']:.3f}")

# PCA
d = load('results/06-dynamics/pca_results.json')
print(f"PCA: total_var={d['total_variance_nm2']:.4f} top3={d['top3_explain_percent']:.1f}% top5={d['top5_explain_percent']:.1f}%")
for ev in d['eigenvalues'][:5]:
    print(f"  PC{ev['pc']}: {ev['eigenvalue_nm2']:.4f} ({ev['percent_variance']:.1f}%) cum={ev['cumulative_percent']:.1f}%")

# WaterMap
d = load('results/01-pocket-discovery/watermap_results.json')
print(f"WaterMap: total={d['total_sites']}, structural={d['structural_sites']}, displaceable={d['displaceable_sites']}, peripheral={d['peripheral_sites']}")
print(f"  frames={d['frames_analyzed']}, grid={d['grid_spacing_A']}A, shell={d['shell_radius_A']}A, threshold={d['threshold']}")

# Allosteric
d = load('results/05-elastic-network/allosteric_network.json')
pd_data = d['pathway_data']
t10 = pd_data['top_10_central']
bc = pd_data['betweenness_centrality']
rmap = d['resid_mapping']
print(f"Allosteric: n_residues={d['n_residues']}, centrality_entries={len(bc)}")
pocket_orig = d['pocket_residues_original']
for i, entry in enumerate(t10[:6]):
    resid_tpr = entry[0]
    cent = entry[1]
    orig = rmap.get(str(resid_tpr), resid_tpr)
    is_pocket = "(POCKET)" if orig in pocket_orig else ""
    print(f"  #{i+1}: resid_tpr={resid_tpr} orig={orig} cent={cent:.4f} {is_pocket}")
intra = d['intra_pocket_correlations']
print(f"  intra_pocket_pairs={len(intra)}, strongest: {intra[0]['res_i']}-{intra[0]['res_j']}={intra[0]['correlation']:.4f}")

# Primary Docking
d = load('results/02-docking/docking_results.json')
relay = d['relay_pocket_docking']
ctrl = d['fostemsavir_site_control']
print(f"Docking: relay={relay['best_affinity_kcal_mol']}, control={ctrl['best_affinity_kcal_mol']}")

# VScreen totals
d1 = load('results/09-drug-design/vscreen_batch1_results.json')
d2 = load('results/09-drug-design/vscreen_batch2_results.json')
all_r = d1['results'] + d2['results']
total = len(all_r)
druggable = sum(1 for r in all_r if r['relay_best_kcal_mol'] < -6.0)
selective = sum(1 for r in all_r if r['delta_kcal_mol'] < -1.0)
vic_r = [r for r in all_r if r['name'] == 'Vicriviroc'][0]
print(f"VScreen: total={total}, druggable={druggable}, selective={selective}")
print(f"  Vic: relay={vic_r['relay_best_kcal_mol']}, ctrl={vic_r['control_best_kcal_mol']}, delta={vic_r['delta_kcal_mol']}")

# Pharmacophore data
d = load('results/01-pocket-discovery/hiv1_pharmacophore_data.json')
relay_t = [t for t in d['pharmacophore_targets'] if 'relay' in t['name'].lower()][0]
feats = relay_t['pharmacophore_features']
ft_counts = {}
for f in feats:
    ft = f['feature_type']
    ft_counts[ft] = ft_counts.get(ft, 0) + 1
print(f"Pharmacophore: relay_features={len(feats)}, hotspots={len(d['binding_hotspots'])}")
for ft, c in sorted(ft_counts.items()):
    print(f"  {ft}: {c}")
