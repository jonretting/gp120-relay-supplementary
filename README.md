# Supplementary Data Package

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23389294.svg)](https://doi.org/10.5281/zenodo.23389294)

**Study:** Identification of a Novel Druggable Pocket in the HIV-1 gp120 C5 Domain at the Conformational Relay Junction

**Author:** Jonathan Michael Retting, Independent Researcher, Los Angeles, CA

**Manuscript:** [gp120-relay-preprint.pdf](gp120-relay-preprint.pdf)

**Primary structure:** PDB 6MEO (cryo-EM, 3.9 A, Shaik et al. *Nature* 2019) -- HIV-1 gp120 (chain G) complexed with human CD4 (chain A) and human CCR5 (chain B).

**License:** CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)

---

## Directory Structure

```
hiv1-gp120-cd4-ccr5/
|
|-- README.md                          This file
|-- VERIFICATION.md                    Integrity verification procedures
|-- SHA256SUMS                         Checksums for all data files
|
|-- figures/                           Publication figures (PNG, 300 DPI)
|   |-- main/                          10 main-text figures (Fig 1--10)
|   |-- supplementary/                 7 supplementary figures (S1--S7 Fig)
|
|-- results/                           All computational results
|   |-- 00-summary/                    Consolidated validation summaries
|   |-- 01-pocket-discovery/           Pocket geometry, druggability, detection
|   |-- 02-docking/                    AutoDock Vina docking + ensemble docking
|   |-- 03-conservation/               Cross-clade survey, entropy, proteome specificity
|   |-- 04-md-simulations/             Per-structure MD analysis (6MEO + 4 SOSIP, 500ns each)
|   |-- 05-elastic-network/            ANM, GNM coupling, allosteric network, alanine scanning
|   |-- 06-dynamics/                   DCCM matrices, PCA/eigenvalues (all 5 structures)
|   |-- 07-free-energy/                FEP/TI (9 legs)
|   |-- 08-binding-energetics/         MM-GBSA, MM-PBSA, IE-MMPBSA, decomposition
|   |-- 09-drug-design/                De novo designs, pharmacophore, ADMET, virtual screen
|   |-- 10-controls/                   CCR5 ECL2/ECL3 docking, V3 loop, selectivity panel
|
|-- scripts/                           All analysis and pipeline scripts
|   |-- setup/                         System preparation, topology building
|   |-- simulation/                    MD launch, PCA computation
|   |-- analysis/                      Post-processing, DCCM, PCA, conservation, validation
|   |-- docking/                       Docking workflows, virtual screening
|   |-- free-energy/                   FEP setup/analysis, MMPBSA
|   |-- drug-design/                   De novo design, ADMET prediction
|
|-- topology/                          GROMACS topology, force field, MDP parameters
|   |-- topol.top                      Master topology (AMBER99SB-ILDN)
|   |-- topol_Protein_chain_*.itp      Per-chain topology includes
|   |-- posre_Protein_chain_*.itp      Position restraint files
|   |-- md.tpr                         Production run binary (reference)
|   |-- em.mdp, nvt.mdp, npt.mdp      FEP equilibration parameters
|   |-- prod.mdp                       FEP production parameters
|
|-- structures/                        PDB and PDBQT coordinate files
|   |-- hiv_6meo.pdb                   Primary structure (6MEO ternary complex)
|   |-- receptor_last.pdb              Final MD frame receptor coordinates
|   |-- bms626529.pdbqt                BMS-626529 ligand (docking input)
|   |-- vicriviroc.pdbqt               Vicriviroc ligand (docking input)
|   |-- gp120_gp41_relay_receptor.*    Relay pocket receptor (docking-prepared)
|
|-- trajectories/                      Raw trajectories (separate download, 67 GB)
|   |-- md-500ns/                      6MEO 500ns MD
|   |-- sosip-500ns/                   4 SOSIP 500ns cross-validation runs
|   |-- fep-raw/                       FEP dH/dl for 288 lambda windows
|   |-- gp120_v3_loop_receptor.*       V3 loop control receptor
|   |-- ccr5_ecl2_receptor.*           CCR5 ECL2 receptor
|   |-- ccr5_ecl3_receptor.*           CCR5 ECL3 receptor
```

---

## Figures

### Main-Text Figures (figures/main/)

| File | Paper Label | Description |
|------|-------------|-------------|
| Fig01_entry_cascade.png | Fig 1 | HIV-1 entry cascade and untargeted relay gap |
| Fig02_structural_overlay.png | Fig 2 | Relay pocket structure in the gp120 C5 domain (3D render) |
| Fig03_conservation_entropy.png | Fig 3 | Sequence conservation of the relay pocket (Shannon entropy) |
| Fig04_relay_pocket_geometry.png | Fig 4 | Relay pocket geometry with druggability annotation |
| Fig05_mmgbsa_decomposition.png | Fig 5 | MM-GBSA per-residue energy decomposition |
| Fig06_rmsd_timeline.png | Fig 6 | Multi-structure relay pocket RMSD over 500 ns |
| Fig07_dccm_comparison.png | Fig 7 | State-dependent DCCM comparison across conformational states |
| Fig08_allosteric_silencing.png | Fig 8 | Allosteric silencing by drug binding |
| Fig09_fep_cycle.png | Fig 9 | FEP thermodynamic cycle |

### Supplementary Figures (figures/supplementary/)

| File | Paper Label | Description |
|------|-------------|-------------|
| S1_Fig_bfactor_comparison.png | S1 Fig | B-factor comparison across gp120 domains |
| S2_Fig_cross_structure_overlay.png | S2 Fig | Cross-structure pocket overlay (5 structures, Kabsch alignment) |
| S3_Fig_rmsf_heatmap.png | S3 Fig | Per-residue RMSF heatmap (0th percentile in 3/4 SOSIP) |
| S4_Fig_pca_variance.png | S4 Fig | PCA variance and Schlitter entropy comparison |
| S5_Fig_ensemble_docking.png | S5 Fig | Ensemble docking score comparison across conformational states |
| S7_Fig_alanine_scanning.png | S7 Fig | Computational alanine scanning (GNM coupling loss + contacts) |

All figures are 300 DPI PNG. Main-text figures are generated by scripts in the figure pipeline (`figures/generate_figures.py` and individual `fig*.py` modules); 3D renders (Fig 2, Fig 4, S2 Fig) were composed from manual screenshots using `figures/compose_3d_figures.py`.

---

## Results -- File Inventory

### 00-summary/ -- Consolidated Validation

| File | Paper Section | Description |
|------|---------------|-------------|
| cross_validation_consolidated.json | Table 4 | Cross-structure validation summary (5 structures) |
| external_validations.json | S3.7 | fpocket, P2Rank, ConCavity cross-validation |
| md_validation_consolidated.json | S3.5 | MD simulation quality metrics |

### 01-pocket-discovery/ -- Pocket Characterization

| File | Paper Section | Description |
|------|---------------|-------------|
| pocket_geometry.json | S3.1 | Pocket dimensions (12.8 x 13.4 x 13.3 A), center, volume |
| druggability_scoring.json | S3.2 | Halgren 2009 5/5 PASS, volume 1,588 A3 |
| bfactor_analysis.json | S3.1 | B-factor percentile analysis (38th percentile) |
| fpocket_analysis.json | S3.3 | fpocket cavity detection (pockets 5 + 19) |
| fpocket_info.txt | S3.3 | Raw fpocket output |
| p2rank_predictions.csv | S3.3 | P2Rank pocket predictions |
| p2rank_residues.csv | S3.3 | P2Rank per-residue scores |
| concavity_scores_chainG.txt | S3.3 | ConCavity surface concavity |
| watermap_results.json | S3.4 (Table 3) | Hydration analysis: 0 structural waters, 8/22 displaceable |
| hiv1_pharmacophore_data.json | S3.1 | 55-feature pharmacophore model |
| gp120_gp41_relay.phar | S3.1 | Relay pocket pharmacophore definition |
| hiv1_binding_hotspots_and_pharmacophores.md | S3.1 | Binding hotspot documentation |

### 02-docking/ -- Molecular Docking

| File | Paper Section | Description |
|------|---------------|-------------|
| docking_results.json | S3.4 | BMS-626529 docking to 6MEO (-6.7 kcal/mol) |
| cross_structure_docking.json | S3.6 | Cross-structure docking to 6MET |
| ensemble_docking_results.json | S3.8 | 20-frame ensemble docking on 6MEO 500ns (mean -6.81) |
| 5U7O_500ns_ensemble_docking.json | S3.11 | SOSIP 5U7O ensemble docking (mean -4.29) |
| 6MU6_500ns_ensemble_docking.json | S3.11 | SOSIP 6MU6 ensemble docking (mean -3.90) |
| gp120_gp41_relay_vina.txt | S3.4 | Raw Vina output for relay pocket |
| ensemble_frames/ | S3.8 | 20 cluster PDB + PDBQT receptor conformations from 500ns MD |
| ensemble_poses/ | S3.8 | 40 docked poses (20 BMS-626529 + 20 Vicriviroc, one per cluster) |

### 03-conservation/ -- Sequence and Structural Conservation

| File | Paper Section | Description |
|------|---------------|-------------|
| cross_clade_survey_validated.json | S3.9 | 531 PDB entries, 360/360 exact RPGGGDMR |
| entropy_validation_proper.json | S3.9 | Shannon entropy from 981 NCBI gp120 sequences |
| resistance_mapping.json | S3.17 | 28 resistance mutations mapped to pocket |
| glycan_shield_analysis.json | S4.5 | Glycan proximity analysis |
| proteome_motif_scan.json | S3.16 | Human proteome RPGGGDMR scan |
| proteome_specificity_results.json | S3.16 | 0/20,431 human proteins contain motif |

### 04-md-simulations/ -- Molecular Dynamics Results

| File | Paper Section | Description |
|------|---------------|-------------|
| md_500ns_analysis.json | S3.5 | 6MEO 500ns: pocket RMSD 1.97 A, RMSF, Rg, energy |
| 5U7O_500ns_analysis.json | S3.11, Table 5 | 5U7O SOSIP 500ns: pocket RMSD 1.80 A |
| 5U7O_md_analysis.json | S3.10, Table 4 | 5U7O SOSIP 20ns analysis |
| 6MU6_500ns_analysis.json | S3.11, Table 5 | 6MU6 SOSIP 500ns: pocket RMSD 1.72 A |
| 6MU6_md_analysis.json | S3.10, Table 4 | 6MU6 SOSIP 20ns analysis |
| 5FYJ_500ns_analysis.json | S3.11, Table 5 | 5FYJ G471A SOSIP 500ns: pocket RMSD 1.29 A |
| 5FYJ_md_analysis.json | S3.10, Table 4 | 5FYJ SOSIP 20ns analysis |
| 4TVP_500ns_analysis.json | S3.11, Table 5 | 4TVP SOSIP 500ns: pocket RMSD 1.98 A |
| 4TVP_md_analysis.json | S3.10, Table 4 | 4TVP SOSIP 20ns analysis |

### 05-elastic-network/ -- Elastic Network Models

| File | Paper Section | Description |
|------|---------------|-------------|
| allosteric_network.json | S3.12 | DCCM + betweenness centrality from 500ns |
| anm_results.json | S3.12 | Anisotropic Network Model eigenvalues (ProDy) |
| gnm_coupling_results.json | S3.10, Table 4 | GNM relay-CD4 coupling (0.88), 140/140 positive |
| alanine_scanning_results.json | S3.18 | Per-residue coupling disruption (max -2.2%) |

### 06-dynamics/ -- DCCM and PCA

| File | Paper Section | Description |
|------|---------------|-------------|
| dccm_500ns_analysis.json | S3.12, Table 6 | 6MEO 500ns DCCM: 28/28 positive (mean +0.79) |
| dccm_pocket.json | S3.12 | Nested-format intra-pocket DCCM |
| apo_dccm_ternary.json | S3.12 | Apo state DCCM for comparison |
| drug_bound_dccm_ternary.json | S3.12 | Drug-bound state DCCM |
| 5U7O_500ns_dccm.json | S3.11, Table 6 | 5U7O DCCM: 19/28 positive (anti-correlated breathing) |
| 6MU6_500ns_dccm.json | S3.11, Table 6 | 6MU6 DCCM: 24/28 positive |
| 5FYJ_500ns_dccm.json | S3.11, Table 6 | 5FYJ G471A DCCM: 28/28 positive (gate seized) |
| 4TVP_500ns_dccm.json | S3.11, Table 6 | 4TVP DCCM: 20/28 positive (resting gate) |
| pca_results.json | S3.12 | 6MEO PCA: total variance, eigenvalues, % explained |
| pca_eigenvalues.xvg | S3.12 | GROMACS eigenvalue output |
| pca_2d_pc1_pc2.xvg | S3.12 | PC1 vs PC2 projection |
| 5U7O_500ns_pca.json | S3.11, Table 6 | 5U7O PCA: 6.98 nm2, Schlitter 1473 J/mol/K |
| 6MU6_500ns_pca.json | S3.11, Table 6 | 6MU6 PCA: 4.45 nm2, Schlitter 1346 J/mol/K |
| 5FYJ_500ns_pca.json | S3.11, Table 6 | 5FYJ PCA: 2.72 nm2 (lowest), Schlitter 1246 J/mol/K |
| 4TVP_500ns_pca.json | S3.11, Table 6 | 4TVP PCA: 5.44 nm2, Schlitter 1438 J/mol/K |

### 07-free-energy/ -- FEP/TI

| File | Paper Section | Description |
|------|---------------|-------------|
| fep_multirep_summary.json | S3.8 | FEP multi-replicate summary: DDG = +3.95 kcal/mol |
| ternary_rep[1-3]_analysis.json | S3.8 | Per-replicate ternary leg (gp120+CD4+drug) |
| binary_rep[1-3]_analysis.json | S3.8 | Per-replicate binary leg (gp120+drug, no CD4) |
| solvent_rep[1-3]_analysis.json | S3.8 | Per-replicate solvent leg (drug in water) |
| fep_mdp/ | S3.8 | FEP GROMACS MDP parameters (8 files: em/nvt/npt/prod x protein/solvent) |
| fep_raw_example/ | S3.8 | Representative FEP dH/dl output (1 window per leg for format reference) |

### 08-binding-energetics/ -- MM-GBSA and MM-PBSA

| File | Paper Section | Description |
|------|---------------|-------------|
| mmpbsa_results.json | S3.13 | MM-PBSA binding free energy (reduced system) |
| mmpbsa_ternary_results.json | S3.14 | MM-PBSA on full ternary complex |
| noncompetitive_mmpbsa_results.json | S3.14 | Non-competitive binding proof (CD4 still binds) |
| mmgbsa_ligand_results.json | S3.13 | MM-GBSA with ligand (v1.1: PDB residue mapping added Sep 28 2026) |
| FINAL_DECOMP_MMPBSA.dat | S3.13 | Per-residue energy decomposition |
| FINAL_RESULTS_MMPBSA.dat | S3.13 | Raw MMPBSA output |
| FINAL_RESULTS_TERNARY.dat | S3.14 | Raw ternary complex MMPBSA |
| MMGBSA_IE.dat | S3.13 | IE-MMPBSA interaction energy |
| MMGBSA_LIGAND.dat | S3.13 | Ligand-specific MMGBSA |
| MMGBSA_NONCOMPETITIVE.dat | S3.14 | Non-competitive MMGBSA output |

### 09-drug-design/ -- Lead Discovery

| File | Paper Section | Description |
|------|---------------|-------------|
| denovo_designs.json | S3.19 | De novo lead compound designs |
| pharmacophore_screen_results.json | S3.20 | Pharmacophore-based screening hits |
| admet_results.json | S3.20 | ADMET predictions for 98 compounds |
| reference_compounds.json | S3.4 | Reference drug properties (BMS-626529 etc.) |
| reference_compounds.smi | S3.4 | SMILES for reference compounds |
| vscreen_batch1_results.json | S3.20 | Virtual screening batch 1 (48 compounds) |
| vscreen_batch2_results.json | S3.20 | Virtual screening batch 2 (48 compounds) |

### 10-controls/ -- Negative Controls

| File | Paper Section | Description |
|------|---------------|-------------|
| ccr5_ecl2.phar | S3.4 | CCR5 ECL2 pharmacophore (control) |
| ccr5_ecl2_vina.txt | S3.4 | ECL2 docking scores |
| ccr5_ecl3.phar | S3.4 | CCR5 ECL3 pharmacophore (control) |
| ccr5_ecl3_vina.txt | S3.4 | ECL3 docking scores |
| gp120_v3_loop.phar | S3.4 | V3 loop pharmacophore (control) |
| gp120_v3_loop_vina.txt | S3.4 | V3 loop docking scores |
| selectivity_panel.json | S3.4 | CCR5/CCR2/CXCR4 off-target selectivity |

---

## Scripts -- Organization

All scripts are Python (.py) or Bash (.sh) and are organized by workflow stage:

| Directory | Count | Purpose |
|-----------|-------|---------|
| scripts/setup/ | 20 | System preparation, topology building, ligand parameterization |
| scripts/simulation/ | 7 | MD launch, SOSIP launch, PCA computation |
| scripts/analysis/ | 48 | Post-processing, DCCM, PCA, conservation, validation, parsing |
| scripts/docking/ | 6 | Ensemble docking, virtual screening, pharmacophore screen |
| scripts/free-energy/ | 20 | FEP setup/analysis, all MMPBSA variants |
| scripts/drug-design/ | 2 | De novo design, ADMET prediction |

---

## Reproduction

### Software Requirements

- GROMACS >= 2021 (2024+ recommended)
- AutoDock Vina 1.2.7
- Python 3.10+ with: MDAnalysis, ProDy, Biopython, RDKit, numpy, scipy
- AmberTools (for MMPBSA and ligand parameterization)

### Script Path Convention

All scripts use a `DATA_ROOT` variable that resolves to this package root directory (`hiv1-gp120-cd4-ccr5/`). By default, scripts locate `DATA_ROOT` relative to their own location:

```python
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
```

To run scripts from a different directory, set `DATA_ROOT` as an environment variable or modify the path in the script header.

### Quick Start

```bash
# Verify data integrity
sha256sum -c SHA256SUMS

# Regenerate figures from supplementary data
cd paper/figures/
python generate_figures.py --all

# Run pocket analysis
python scripts/analysis/pocket_analysis.py

# Run cross-clade conservation
python scripts/analysis/compute_entropy_proper.py
```

### Trajectory Data (Separate Download)

Raw simulation trajectories (67 GB) are available as five self-contained packages from the trajectory download portal at **https://hiv1-relay-data.jretting.workers.dev**. To integrate them with this supplementary package, extract each package into a `trajectories/` subdirectory at the package root:

```bash
# From the supplementary package root (hiv1-gp120-cd4-ccr5/):
mkdir -p trajectories
cd trajectories

# Download whichever packages you need, then extract:
tar --zstd -xf hiv1-relay-md-500ns.tar.zst          # -> trajectories/md-500ns/
tar --zstd -xf hiv1-relay-sosip-4x500ns.tar.zst      # -> trajectories/sosip-500ns/
tar --zstd -xf hiv1-relay-fep-288lambda.tar.zst       # -> trajectories/fep-raw/
```

After extraction, the resulting structure is:

```
hiv1-gp120-cd4-ccr5/
├── results/
├── scripts/
├── topology/
├── structures/
└── trajectories/                        <- you create this
    ├── md-500ns/6meo/md_500ns.xtc         6MEO 500ns (5.7 GB)
    ├── sosip-500ns/5u7o/md_500ns.xtc      5U7O 500ns (6.5 GB)
    ├── sosip-500ns/6mu6/md_500ns.xtc      6MU6 500ns (5.6 GB)
    ├── sosip-500ns/5fyj/md_500ns.xtc      5FYJ 500ns, G471A (6.1 GB)
    ├── sosip-500ns/4tvp/md_500ns.xtc      4TVP 500ns (5.6 GB)
    └── fep-raw/{ternary,binary,solvent}/   288 dH/dl XVG files (1.1 GB)
```

**Which scripts need which trajectories:**

| Script | Trajectory package needed |
|--------|--------------------------|
| `scripts/analysis/03_analyze_md.py` | md-500ns |
| `scripts/analysis/analyze_500ns_interim.py` | md-500ns |
| `scripts/analysis/compute_dccm_500ns.py` | md-500ns |
| `scripts/analysis/analyze_sosip_500ns.py` | sosip-4x500ns |
| `scripts/analysis/compute_sosip_dccm_500ns.py` | sosip-4x500ns |
| `scripts/analysis/run_sosip_mechanistic.sh` | sosip-4x500ns |
| `scripts/free-energy/run_fep.sh` | fep-288lambda |

Scripts that analyze pre-computed results (JSON files in `results/`) do not require trajectory downloads.

**Verification:** SHA-256 checksums for all trajectory packages are provided at the download portal and in the `SHA256SUMS` file included in this repository.

### MD Simulation Protocol

The primary simulation uses AMBER99SB-ILDN force field with TIP3P water, 150 mM NaCl, PME electrostatics (1.2 nm cutoff), V-rescale thermostat at 300 K, C-rescale barostat at 1 bar. Production runs use 2 fs timestep with h-bond constraints. Topology files are in `topology/`.

---

## Aggregate Statistics

- **Total simulation time:** ~4,000 ns (500ns x 5 structures + 100ns + 1,440ns FEP)
- **Structures analyzed:** 6 PDB entries (6MEO, 6MET, 5U7O, 6MU6, 5FYJ, 4TVP)
- **Cross-clade survey:** 531 PDB entries, 981 NCBI sequences
- **Virtual screening:** 96 compounds across 2 batches
- **FEP replicates:** 9 legs (3 ternary + 3 binary + 3 solvent), 32 lambda windows each

---

---

**Author:** Jonathan Michael Retting, Independent Researcher, Los Angeles, CA
**License:** CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
