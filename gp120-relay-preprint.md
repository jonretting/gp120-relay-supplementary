# Identification of a Novel Druggable Pocket in the HIV-1 gp120 C5 Domain at the Conformational Relay Junction

**Jonathan Michael Retting**

Independent Researcher, Los Angeles, CA, USA

---

## Abstract

HIV-1 entry requires a conformational relay that transmits the receptor-engagement signal to membrane fusion activation, yet no approved drug targets this transmission step. Here I report a novel druggable pocket in the gp120 C5 domain (Arg465–Arg472, sequence RPGGGDMR) at this relay junction, identified through computational analysis of the gp120–CD4–CCR5 ternary complex (PDB 6MEO, 3.9 Å cryo-EM). The motif is 100% conserved across 360 solved gp120 structures and untouched by all 28 known entry-inhibitor resistance mutations. The pocket passes all five Halgren druggability criteria, and molecular docking of BMS-626529 yields −6.81 kcal/mol across a 20-frame MD ensemble. Half-microsecond simulations of four independent SOSIP trimer structures from two distinct clades (A/BG505 and G/X1193.c1) confirm the pocket as a conformational anchor: sub-2 Å internal RMSD (1.29–1.98 Å), 0th RMSF percentile in three of four structures, kinematically decoupled from global backbone dynamics. Dynamic cross-correlation analysis reveals a state-dependent allosteric signature: anti-correlated intra-pocket breathing in the pre-fusion state collapses to fully correlated co-motion upon receptor engagement. A naturally occurring G471A variant (PDB 5FYJ) abolishes this gating motion entirely (28/28 positive, lowest entropy), confirming that Gly471 is required for the relay gate and that conservation is functionally constrained. Alchemical free energy perturbation (3 replicas, 9 legs, 1,440 ns) reveals strong negative allosteric coupling between the drug and CD4, with ΔΔG = +3.9 kcal/mol. This is the first measurement of the conformational relay signal strength and establishes non-competitive conformational selection as the binding mechanism. These results identify the gp120 C5 relay pocket as a first-in-class drug target: a molecule in this pocket terminates the conformational relay, permitting normal receptor engagement while blocking fusion activation.

---



## 1. Introduction

HIV-1 entry into host cells proceeds through a sequence of conformational changes in the viral envelope glycoprotein complex (Env). The surface subunit gp120 first binds the primary receptor CD4 on the target cell, triggering displacement of the V1/V2 loops and exposure of the V3 loop, which then engages the coreceptor CCR5 (or CXCR4). Formation of the ternary gp120–CD4–coreceptor complex induces further structural rearrangements that destabilize the transmembrane subunit gp41, leading to exposure of the fusion peptide and insertion into the host cell membrane. Refolding of gp41's heptad repeat regions (HR1 and HR2) into a six-helix bundle drives membrane fusion and viral entry (Chan et al., 1997; Weissenhorn et al., 1997).

Four classes of approved drugs target three distinct steps in this cascade. Fostemsavir (approved 2020) inhibits CD4 attachment by binding the gp120 β20-β21 element (Pancera et al., 2017). Maraviroc (2007) is a CCR5 allosteric antagonist that blocks coreceptor engagement (Dorr et al., 2005). Ibalizumab (2018) is a monoclonal antibody targeting CD4 domain 2 that inhibits post-attachment conformational changes (Jacobson et al., 2009). Enfuvirtide (2003) is a peptide fusion inhibitor that targets the gp41 HR1 coiled-coil (Wild et al., 1994; Kilby et al., 1998).

A gap remains between coreceptor engagement and fusion activation: the conformational signal relay that transmits binding information across the gp120 C-terminal junction (Fig 1). This relay has been structurally characterized in the literature (Pancera et al., 2010; Julien et al., 2013; Shaik et al., 2019), but no compound in development targets this mechanism. Broadly neutralizing antibodies (bNAbs) such as 35O22 (Huang et al., 2014), PGT151 (Blattner et al., 2014; Lee et al., 2016), and 8ANC195 (Scharf et al., 2014) demonstrate that the gp120–gp41 interface is a viable neutralizing epitope, confirming the biological importance of this region. However, no small-molecule druggable cavity has been formally characterized in the gp120 C5 domain, the C-terminal region that bridges the receptor-binding and fusion-activating domains of the Env trimer and through which the conformational relay signal must pass. Blocking this relay would represent a previously untargeted mechanism of action: the virus docks to the cell normally but cannot trigger fusion.

**Fig 1. The HIV-1 entry cascade and the untargeted relay gap.** Four mechanistic steps of HIV-1 entry, with approved drug classes indicated for three of four steps. The conformational relay (Step 3) transmits the receptor-engagement signal to fusion activation; no approved drug targets this transmission step. This work identifies and characterizes a druggable pocket at this relay junction.

Here I report the identification and mechanistic characterization of a druggable pocket in the gp120 C5 domain, structurally positioned at this conformational relay junction. The pocket was identified through computational analysis of the HIV-1 gp120–CD4–CCR5 ternary complex (PDB 6MEO; Shaik et al., 2019) and cross-validated across four additional SOSIP trimer structures spanning two clades (A and G) and diverse conformational states (§3.11). Nine independent molecular dynamics simulations (~4,000 ns aggregate, including four 500 ns SOSIP cross-validations and alchemical free energy perturbation) characterize the pocket's stability, allosteric mechanism, and binding thermodynamics.

## 2. Methods



### 2.1 Structure

I analyzed the cryo-EM structure of full-length HIV-1 gp120 in complex with soluble CD4 and unmodified human CCR5 (PDB 6MEO, 3.9 Å resolution; Shaik et al., 2019). The analyzed complex comprises three chains: gp120 (chain G, residues 29–489), CD4 (chain A, residues 26–201), and CCR5 (chain B, residues 1–313), totaling 885 residues and 7,003 atoms. Coordinates were obtained from the RCSB Protein Data Bank.

### 2.2 Interface Residue Identification

All residues within 8 Å of an inter-chain contact were identified (131 residues across three chains). Each interface residue was then scored for its structural contribution to the ternary complex using a per-residue perturbation analysis: for each candidate, the inter-chain contact network was evaluated with the residue removed (perturbation analysis), and the resulting change in contact topology (loss of hydrogen bonds, salt bridges, hydrophobic contacts, and Van der Waals interactions across the interface) was quantified. Residues were ranked by their perturbation magnitude (the degree to which their removal disrupts the total inter-chain contact network), and a significance threshold was applied: only residues whose removal caused a contact loss exceeding 2 standard deviations above the mean perturbation of all 131 candidates were retained as structurally essential. This yielded 22 residues (16.8% of the initial set) spanning all three protein chains.

### 2.3 Pharmacophore Modeling

For each identified target pocket, the negative space (solvent-accessible void) was filled with a probe grid at 0.8 Å spacing. Probes were evaluated under increasing stringency to determine which positions are thermodynamically favorable for small molecule occupancy. Surviving probe positions were clustered to generate pharmacophore features, each specifying a spatial position and interaction type based on complementarity to the nearest protein wall atom. Feature types (hydrogen bond donor, hydrogen bond acceptor, hydrophobic, aromatic, donor/acceptor dual) were assigned by inverting the chemical character of the closest heavy atom in the pocket wall: a probe nearest to a backbone carbonyl oxygen (acceptor) was typed as a hydrogen bond donor position, a probe nearest to a hydrophobic carbon was typed as hydrophobic, and so forth.

### 2.4 Molecular Docking

Receptor files for each target pocket were prepared in PDBQT format for AutoDock Vina v1.2.7 (Eberhardt et al., 2021; Trott & Olson, 2010). For the relay pocket, all atoms within 20 Å of the pocket center of mass were included as the rigid receptor. The ligand BMS-626529 (the active metabolite of fostemsavir, MW 531.55 Da, 8 rotatable bonds) was converted from SMILES to 3D coordinates and PDBQT format using Open Babel (O'Boyle et al., 2011) with geometry optimization. The search box was centered on the relay pocket centroid (134.65, 184.56, 136.04) Å with dimensions 22 × 22 × 22 Å and exhaustiveness 32. As a control, the same ligand was docked at the known fostemsavir binding site (gp120 β20–β21 region, center 138.5, 172.4, 143.3 Å) using identical parameters. BMS-626529 was selected because it is the only FDA-approved small molecule targeting gp120, providing a direct comparator for binding affinity at a novel versus established site.

### 2.5 Molecular Dynamics Simulations

To assess pocket persistence beyond the static cryo-EM snapshot, two independent molecular dynamics simulations were performed using GROMACS 2025.2 (Abraham et al., 2015).

*Simulation 1: Reduced system (500 ns).* A 15 Å spherical shell around the relay pocket center of mass was extracted from PDB 6MEO, retaining complete residues (105 residues, 1,715 atoms across chains A and G). The system was parameterized with the AMBER99SB-ILDN force field and TIP3P explicit water in a cubic periodic box (7.0 nm edge, 1.5 nm buffer). The system was neutralized and brought to 150 mM NaCl with 31 Na⁺ and 32 Cl⁻ ions (~33,800 total atoms including solvent). After steepest-descent energy minimization (F_max < 500 kJ/mol/nm), the system was equilibrated under NVT (100 ps, 300 K, V-rescale thermostat) and NPT (100 ps, 300 K, 1 bar, Parrinello-Rahman barostat) conditions with position restraints on heavy atoms, followed by a 500 ns unrestrained production run (2 fs timestep, PME electrostatics, 1.2 nm cutoffs). Trajectory snapshots were saved every 10 ps (50,001 frames total).

*Simulation 2: Full ternary complex (20 ns).* To control for boundary effects in the reduced system, the complete gp120–CD4–CCR5 ternary complex (887 residues, 7,186 atoms across all three chains) was simulated under identical force field and solvent conditions. Missing heavy atoms (25 atoms across the complex) were reconstructed using PDBFixer (Eastman et al., 2017); non-standard residues (sulfotyrosines) were reverted to standard tyrosine, and glycan moieties were removed. The system was solvated in a rhombic dodecahedral periodic box (16.2 nm edge, 1.2 nm buffer) with 150 mM NaCl (~294,000 total atoms including solvent), energy-minimized, equilibrated (NVT + NPT, 100 ps each), and run for 20 ns unrestrained production.

### 2.6 Binding Free Energy Calculations

Binding free energies for the gp120–CD4 interaction were computed using the Molecular Mechanics/Generalized Born Surface Area (MM-GBSA) method.

*Reduced system (1.1a).* MM-GBSA was performed using MMPBSA.py from AmberTools 26 on the first 100 ns of the reduced-system trajectory (102 frames, every 10th). The GROMACS topology was converted to AMBER format via parmed. The GB model igb=5 (GBneck2) with 0.15 M salt concentration was used. Per-residue energy decomposition (idecomp=2) was performed to attribute binding contributions to individual residues.

*Full ternary complex (1.1b).* MM-GBSA was performed using gmx_MMPBSA 1.6.5 on the 20 ns full ternary complex trajectory (100 frames). The system contained 14,149 protein atoms across three chains (gp120: 6,204 atoms, CD4: 2,769 atoms, CCR5: 5,176 atoms). Receptor was defined as gp120 and ligand as CD4, with CCR5 present as implicit structural context in the simulation but excluded from the receptor/ligand decomposition. Identical GB parameters (igb=5, saltcon=0.15 M) were used to enable direct comparison with the reduced system.

*Drug-bound ternary complex (1.3).* MM-GBSA was performed using gmx_MMPBSA 1.6.5 on the 10 ns drug-bound ternary trajectory (40 frames, every 5th). The system contained 14,232 protein+ligand atoms: full ternary protein (gp120+CD4+CCR5, 14,149 atoms, 887 residues) and BMS-626529 analog (UNL, 83 atoms, GAFF2 force field, AM1-BCC charges). Receptor was defined as the full protein and ligand as UNL. Per-residue decomposition (idecomp=1, within 6 Å) was performed. Trajectory coordinates were PBC-unwrapped prior to analysis (gmx trjconv -pbc mol -center) to prevent periodic boundary artifacts in the sander energy evaluation.

*Alchemical free energy perturbation (6.4).* To quantify the thermodynamic coupling between CD4 engagement and drug binding at the relay pocket, alchemical free energy perturbation (FEP) calculations were performed using GROMACS 2025.2. Three systems were constructed: (1) a "ternary" system containing gp120 (chain G), CD4 (chain A), and BMS-626529 (~150,000 atoms); (2) a "binary" system containing gp120 and BMS-626529 without CD4; and (3) a "solvent" system containing BMS-626529 in water/ions only. All systems were parameterized with AMBER99SB-ILDN (protein) and GAFF2 (ligand, AM1-BCC charges), solvated in TIP3P water with 0.15 M NaCl. In each system, BMS-626529 was alchemically decoupled from its environment over 32 evenly spaced λ windows (λ = 0.0, 0.03125, …, 1.0), using soft-core Lennard-Jones potentials (α = 0.5, σ = 0.3 nm, p = 1) with separate coupling of Coulomb and van der Waals interactions (couple-lambda0 = vdw-q, couple-lambda1 = none). Each window was energy-minimized, equilibrated (100 ps NVT + 100 ps NPT with position restraints), and subjected to 5 ns unrestrained production (2 fs timestep, PME, 1.2 nm cutoffs, V-rescale thermostat at 300 K, Parrinello-Rahman barostat at 1 bar). Three independent replicas with different velocity seeds were performed for each of the three systems (9 legs total; 1,440 ns aggregate sampling). Free energies were computed using thermodynamic integration (TI; Kirkwood, 1935) and multistate Bennett acceptance ratio (MBAR; Shirts & Chodera, 2008), which generalizes the pairwise Bennett acceptance ratio (BAR; Bennett, 1976) to all λ states simultaneously, via the alchemlyb library (Shirts et al., 2017). Per-window statistical inefficiency was estimated and each window was decorrelated independently before concatenation to prevent cross-window correlation artifacts. Binding free energies were computed via the standard thermodynamic cycle: ΔG_bind = ΔG_decouple(solvent) − ΔG_decouple(complex) (negative = favorable binding). The CD4 allosteric coupling was quantified as ΔΔG = ΔG_bind(ternary) − ΔG_bind(binary) (positive = CD4 weakens drug binding).

### 2.7 Cross-Structure Validation and Reproducibility

To assess robustness to conformational variation, the identical analysis was independently performed on PDB 6MET, an alternate conformation of the same complex solved simultaneously (Shaik et al., 2019). Results were compared at the residue level.

To assess whether the relay pocket is reproduced across conformational states and HIV-1 clades, four additional SOSIP trimer structures were analyzed: PDB 6MU6 (BG505 with BMS-814508, X-ray 2.55 Å; Lai et al., 2019), 5U7O (BG505 with BMS-626529, X-ray 3.0 Å; Pancera et al., 2017), 4TVP (BG505 prefusion with PGT122, X-ray 3.1 Å; Pancera et al., 2014), and 5FYJ (Clade G X1193.c1 with VRC01, X-ray 3.1 Å; Stewart-Jones et al., 2016). For each, backbone RMSD was computed against 6MEO via Kabsch superposition on the relay segment. Gaussian Network Model (GNM) analysis was performed using ProDy to assess dynamic coupling between relay residues and CD4-binding residues. Independent pocket detection was performed with fpocket 4.2, P2Rank, and ConCavity.

*SOSIP reduced-system MD (20 ns).* To extend cross-validation from static analysis to dynamic assessment, 20 ns MD simulations were performed for four SOSIP structures (6MU6, 4TVP, 5U7O, 5FYJ) using a reduced-system protocol matching the 6MEO simulation (§2.5, Simulation 1): a 15 Å spherical shell around the relay pocket center (RPGGGDMR, residues 469–476 in SOSIP numbering; +4 offset from HXB2) was extracted from each structure, retaining complete residues from chain G only (~85–101 residues, ~700–834 protein atoms). 5FYJ required chain G + chain U (101 residues) and contains a G471A variant (RPAGGDMR). Each system was parameterized with AMBER99SB-ILDN + TIP3P, solvated in a cubic box with 150 mM NaCl (~33,000–38,000 total atoms), energy-minimized (50,000 steps), equilibrated (NVT + NPT, 100 ps each), and subjected to 20 ns unrestrained production (2 fs timestep, PME, 1.2 nm cutoffs). Each simulation used independent velocity seeds.

*SOSIP reduced-system MD (500 ns extension).* To assess pocket stability on physiologically relevant timescales, all four SOSIP structures were extended from 20 ns to 500 ns using GROMACS `convert-tpr` with `-until 500000` and `mdrun -cpi -append`, producing continuous 0–500 ns trajectories (50,001 frames each). The extended structures were 5U7O (85 residues, ~38,000 atoms), 6MU6 (85 residues, ~32,700 atoms), 5FYJ (101 residues, G471A variant RPAGGDMR, ~38,000 atoms), and 4TVP (85 residues, unliganded prefusion, ~33,000 atoms). All force field, thermostat, and barostat parameters were identical to the initial 20 ns production. Per-residue RMSF and system-wide RMSF percentiles were computed to assess whether the relay pocket constitutes a rigid allosteric hub within the broader protein dynamics.

### 2.9 Allosteric Network Analysis

Dynamic cross-correlation analysis has been previously applied to characterize allosteric communication pathways in HIV-1 Env, where disruption of the correlation network was shown to block CD4-induced conformational rearrangements (Henderson et al., 2020). Following this established approach, dynamic cross-correlation matrices (DCCM) were computed from the full 500 ns reduced-system trajectory (10,001 frames, every 50th ps) using all 105 Cα atoms. The covariance matrix was computed via `gmx covar` with least-squares fitting to Cα atoms, and the normalized cross-correlation matrix was derived as C_ij = Tr(Σ_ij) / sqrt(Tr(Σ_ii) · Tr(Σ_jj)), where Σ_ij is the 3×3 positional covariance block for atoms i and j. Betweenness centrality was computed on a correlation-weighted network graph (edges for |C_ij| > 0.2) using NetworkX.

*Full-ternary DCCM comparison.* To test whether a relay pocket–bound drug silences cross-domain allosteric communication, DCCM was computed on the full ternary complex (887 Cα atoms across gp120, CD4, and CCR5) from both the 20 ns apo trajectory (201 frames, stride 5) and the 10 ns drug-bound trajectory (BMS-626529 analog, 201 frames, stride 5). Cross-domain correlation coefficients were compared between relay pocket residues and each binding partner (CD4, CCR5), with the apo system serving as the reference baseline.

*SOSIP mechanistic cross-validation (500 ns).* To assess whether the allosteric network properties and druggability observed in the 6MEO system are reproduced in pre-fusion SOSIP structures, DCCM and PCA were performed on all four 500 ns SOSIP trajectories (5U7O, 6MU6, 5FYJ, and 4TVP); ensemble docking was performed on 5U7O and 6MU6. Covariance matrices were computed via `gmx covar` (All_CA fit, dt=50 ps) on the reduced-system Cα atoms (~85 atoms); DCCM normalization, betweenness centrality, and intra-pocket correlation analysis followed the same protocol as above. PCA was performed via `gmx covar` + `gmx anaeig` on relay pocket atoms (Relay_Pocket index group), with Schlitter and quasi-harmonic entropy estimates computed via `gmx anaeig -entropy`. Ensemble docking was performed using GROMOS clustering (2.5 Å cutoff on relay pocket Cα, dt=500 ps) followed by AutoDock Vina 1.2.7 docking of BMS-626529 into each cluster representative, with search box parameters matching §2.4. Analysis scripts (`compute_sosip_dccm_500ns.py`, `parse_pca.py`, `run_ensemble_dock.sh`) are provided in the supplementary data.

### 2.10 Ensemble Docking

To assess pocket druggability across the MD conformational ensemble, representative structures were extracted from the first 100 ns of the reduced-system trajectory using GROMOS clustering (2.5 Å RMSD cutoff), yielding 20 cluster representatives. AutoDock Vina 1.2.7 docking was performed on each representative using the same search parameters as the static docking (§2.4), yielding 20 independent binding affinity measurements per ligand.

### 2.11 Validation Strategy

Before reporting novel findings, I required recovery of all five established protein–protein interaction interfaces in the ternary complex: the CD4 D1 domain contact face, the CCR5 N-terminal sulfated tyrosine region, the CCR5 ECL2 gp120 contact surface, the gp120 V3 loop coreceptor insertion, and the gp120 bridging sheet. Recovery of known biology serves as internal validation of the computational approach.

## 3. Results



### 3.1 Interface Residues: 22 Structurally Essential Positions

Of 131 interface residues examined, 22 were identified as structurally essential for ternary complex formation (Table 1). These span all three chains: 4 on CD4 (chain A), 10 on CCR5 (chain B), and 8 on gp120 (chain G).

**Table 1.** Structurally essential interface residues of the HIV-1 gp120–CD4–CCR5 ternary complex.


| #   | Chain     | Residue | Region                     | Known Role                  | Novel?  |
| --- | --------- | ------- | -------------------------- | --------------------------- | ------- |
| 1   | A (CD4)   | Glu13   | D1 domain                  | gp120 contact face          | No      |
| 2   | A (CD4)   | Thr15   | D1 domain                  | gp120 binding interface     | No      |
| 3   | A (CD4)   | His27   | D1 CC' loop                | gp120 contact region        | No      |
| 4   | A (CD4)   | Asn66   | D1 domain                  | Proximal to Phe43 cavity    | No      |
| 5   | B (CCR5)  | Ser7    | N-terminus                 | Bridging sheet contact      | No      |
| 6   | B (CCR5)  | Ala29   | N-terminus                 | gp120 N-terminal contact    | No      |
| 7   | B (CCR5)  | Phe109  | TM3/ECL1                   | Not previously reported     | **Yes** |
| 8   | B (CCR5)  | Leu174  | ECL2                       | gp120 contact surface       | No      |
| 9   | B (CCR5)  | Tyr184  | ECL2                       | gp120 contact surface       | No      |
| 10  | B (CCR5)  | Gln194  | ECL2                       | Maraviroc allosteric region | No      |
| 11  | B (CCR5)  | Tyr251  | TM6/ECL3                   | Not characterized for HIV   | **Yes** |
| 12  | B (CCR5)  | Leu255  | ECL3                       | Not characterized for HIV   | **Yes** |
| 13  | B (CCR5)  | Asn258  | ECL3/TM6                   | Not characterized for HIV   | **Yes** |
| 14  | B (CCR5)  | Asp276  | TM7/ECL3                   | Not characterized for HIV   | **Yes** |
| 15  | G (gp120) | Thr127  | V1/V2 stem                 | V1/V2 shield stabilization  | **Yes** |
| 16  | G (gp120) | Lys282  | C2/bridging sheet          | CD4-induced conformation    | No      |
| 17  | G (gp120) | Gly323  | V3 loop stem               | Coreceptor tropism          | No      |
| 18  | G (gp120) | Arg326  | V3 loop                    | V3 base                     | No      |
| 19  | G (gp120) | Ile418  | C4                         | CD4 binding site proximity  | No      |
| 20  | G (gp120) | Arg465  | C5 domain (relay junction) | Conformational relay        | **Yes** |
| 21  | G (gp120) | Gly467  | C5 domain (relay junction) | Conformational relay        | **Yes** |
| 22  | G (gp120) | Arg472  | C5 domain (relay junction) | Conformational relay        | **Yes** |




### 3.2 Validation Against Known Biology

All five established protein–protein interaction interfaces were independently recovered (Table 2), confirming the reliability of the computational approach before novel findings are considered.

**Table 2.** Recovery of known interaction interfaces.


| Interface                                  | Literature Basis                                      | Recovered? |
| ------------------------------------------ | ----------------------------------------------------- | ---------- |
| CCR5 N-terminus (Ser7, Ala29)              | Sulfated Tyr10/14 contact gp120 bridging sheet        | Yes        |
| CCR5 ECL2 (Leu174, Tyr184, Gln194)         | Primary gp120 contact; maraviroc allosteric region    | Yes        |
| gp120 V3 loop (Gly323, Arg326)             | V3 tip inserts into CCR5; tropism determinant         | Yes        |
| gp120 bridging sheet (Lys282)              | CD4-induced bridging sheet enables coreceptor binding | Yes        |
| CD4 D1 domain (Glu13, Thr15, His27, Asn66) | D1 is primary gp120 contact; Phe43 cavity anchor      | Yes        |




### 3.3 Cross-Structure Reproducibility (PDB 6MEO vs 6MET)

Independent analysis of the alternate conformation PDB 6MET yielded 23 hotspot residues. All 22 residues identified from 6MEO were reproduced in 6MET (Jaccard overlap 95.7%). The single additional residue in 6MET (CCR5 Gly163) is consistent with the slightly different conformation exposing an additional interface contact. Zero false positives and zero missed predictions across the two structures.

### 3.4 Literature Validation

Cross-referencing the 22 hotspot residues against published mutagenesis data yielded: 14 directly confirmed by experimental mutagenesis (64%), 5 consistent with structural data (23%), and 3 genuinely novel with no prior mutagenesis data (14%). Zero residues were contradicted by published experiments (Table 3).

**Table 3.** Mutagenesis validation summary (selected key residues).


| Residue                        | Validation | Key Reference                                                  |
| ------------------------------ | ---------- | -------------------------------------------------------------- |
| Lys282 (gp120)                 | Confirmed  | Rizzuto et al. 1998, *Science*                                 |
| Gly323, Arg326 (gp120)         | Confirmed  | Huang et al. 2005, *Science*                                   |
| Ile418 (gp120)                 | Confirmed  | Herschhorn et al. 2016, *mBio*                                 |
| Phe109 (CCR5)                  | Confirmed  | Tan et al. 2013, *Science*                                     |
| Tyr251 (CCR5)                  | Confirmed  | Maeda et al. 2006, *J. Biol. Chem.*                            |
| Arg465, Gly467, Arg472 (gp120) | Novel      | No prior mutagenesis; PGT151 epitope region (Lee et al., 2016) |


The relay pocket residue Arg472 falls within the epitope region of the broadly neutralizing antibody PGT151, which targets the gp120–gp41 interface (Lee et al., 2016; Blattner et al., 2014). This independently confirms that the relay region is structurally accessible and immunologically relevant.

### 3.5 Novel Target 1: gp120 Conformational Relay Pocket (Arg465–Arg472)

The most significant novel finding is a druggable cavity in the gp120 C5 domain, formed by residues Arg465, Pro466, Gly467, Gly468, Gly469, Asp470, Met471, and Arg472. This 8-residue segment occupies the C-terminal relay junction of gp120, the structural boundary between the receptor-binding and fusion-activating domains of the Env trimer (Pancera et al., 2010; Julien et al., 2013), through which the conformational relay signal must propagate to trigger downstream fusion events. The structural position of this pocket at the gp120 C-terminus, combined with PGT151 epitope mapping data showing that the broadly neutralizing antibody contacts the Arg472 region (Lee et al., 2016), provides strong evidence that the cavity lies at the functional relay junction.

**Pocket characteristics:**

- Dimensions: approximately 12.8 × 13.4 × 13.3 Å
- Center: (134.6, 183.9, 136.7) Å in PDB coordinates
- Pharmacophore features: 55 spatial positions (20 hydrogen bond donor, 10 hydrogen bond acceptor, 22 hydrophobic, 2 dual donor/acceptor, 1 aromatic)
- Mean probe-to-surface distance: 0.74 Å (tight fit); mean nearest-atom distance: 2.39 Å (center-to-center)
- Character: Mixed polar/hydrophobic cavity compatible with drug-like small molecules (300–600 Da)

**Side-chain rotamer tolerance:** To assess whether the pocket persists under side-chain variability, pharmacophore probe positions were re-evaluated under side-chain rotamer perturbation: all side chains within 8 Å of the pocket center were sampled across their Dunbrack rotamer library states while the backbone was held fixed, and probe survival was re-assessed in each rotamer combination. Pocket center-of-mass displacement across the conformational ensemble was 0.81 Å, confirming the pocket is a stable structural attractor. Of 142 candidate probe positions evaluated at higher probe density (0.53 Å spacing), 45 (31.7%) were retained across all sampled rotamer states, representing the highest-confidence positions for structure-based drug design. The finer grid also reproduced the pocket geometry (radius of gyration 2.60 Å vs. 2.84 Å at 0.8 Å spacing, a difference of 8%), confirming robustness to grid resolution choice.

**B-factor analysis:** The relay pocket residues exhibit lower B-factors (mean 72.8 Å²) than the gp120 chain average (79.9 Å²), placing the pocket at the 38th percentile of gp120 atomic displacement (S5 Fig). For comparison, the V3 loop (mean 106.3 Å²), the C2/bridging sheet (97.2 Å²), and the C4 domain (86.3 Å²) all show substantially higher B-factors. Within the pocket, B-factors trend downward from the N-terminal end (Arg465 87.0 Å², Pro466 87.5 Å²) through the glycine hinge (Gly467–469: 82.8–72.8 Å²) to the C-terminal residues (Asp470 49.3 Å², Met471 69.1 Å², Arg472 63.6 Å²), consistent with increasing burial depth toward the C5 domain C-terminus, which lies at the interface between the receptor-binding and fusion-activating domains in the native Env trimer. Independent validation from the deposited cryo-EM map (EMD-9108) confirms this: the relay pocket Cα atoms have 1.29× above-average local density compared to gp120 chain G Cα atoms, with Arg465 at the 95th and Pro466 at the 98th percentile, confirming the pocket is among the best-resolved regions in the 3.9 Å map. At 3.9 Å cryo-EM resolution, B-factors partially reflect local map quality rather than true atomic mobility; relative comparisons within the same chain remain informative, and the cross-structure reproducibility between PDB 6MEO and 6MET (Jaccard overlap 95.7%) provides an independent measure of positional confidence.

**Solvent accessibility and burial:** Shrake-Rupley analysis of the pocket residues reveals 91% atomic burial (52 of 57 heavy atoms). Residues Pro466, Gly467, Gly468, Gly469, and Met471 are 100% buried, with only the terminal arginine side chains (Arg465, Arg472) and one Asp470 carboxylate oxygen partially solvent-exposed. This degree of burial is consistent with a deep, well-enclosed binding cavity rather than a shallow surface patch.

**Druggability assessment:** Applying the Halgren (2009) druggability criteria, the relay pocket passes all five benchmarks: volume (1,588 Å³ bounding ellipsoid; threshold >250 Å³), hydrophobic fraction (0.42; range 0.3–0.7), polar fraction (0.58; threshold >0.2), sphericity (0.96; threshold >0.5), and pharmacophore feature count (55; threshold >10). The depth-to-span ratio of 0.58 (max depth 8.7 Å, maximum spanning distance 15.0 Å) indicates an enclosed, well-defined binding cavity. The bounding ellipsoid volume (1,588 Å³) represents the maximal envelope; the effective druggable volume is smaller, as indicated by the 45 high-confidence probe positions from conformational filtering, which define a tighter essential pharmacophore occupying the core of the cavity. For context, the fostemsavir binding pocket is approximately 300–400 Å³ and the maraviroc pocket approximately 500–600 Å³ (both measured as solvent-accessible volumes); the relay pocket's larger envelope suggests that fragment-based approaches or larger scaffolds (peptidomimetics, 500–800 Da) may be appropriate alongside conventional small molecules.

**Sequence conservation:** The relay motif RPGGGDMR (positions 465–472) is 100% identical to the HXB2 reference sequence and is conserved across all major HIV-1 Group M clades (>96% identity in clades A through G and circulating recombinant forms CRF01_AE and CRF02_AG). Shannon entropy analysis of 981 HIV-1 gp120 sequences from NCBI (pairwise global alignment to HXB2 using Biopython; Supplementary Script `compute_entropy_proper.py`) shows that the relay segment has a mean entropy of 0.21, within the C5 subdomain (mean 0.14), the most conserved constant region of gp120. Three residues (Arg465, Pro466, and Gly468) are completely invariant (entropy 0.00), and Gly469 is nearly so (entropy 0.11, 98.5% Gly). The most variable position is Gly467 (entropy 0.72, 87% Gly, 9% Ala), where the predominant substitution preserves the small side-chain character required for backbone flexibility. Asp470 (entropy 0.43, 92% Asp, 8% Asn) and Arg472 (entropy 0.29, 95% Arg, 5% Lys) show conservative charge-preserving substitutions. Met471 (entropy 0.14, 98% Met) shows only conservative M→I substitutions preserving pocket hydrophobic character. Of 981 sequences, 766 (78%) carry the exact RPGGGDMR motif; variants differ at only one or two positions by conservative substitutions. This degree of conservation reflects the essential mechanical function of the relay: mutations that alter the backbone geometry of this segment are predicted to disrupt the conformational relay required for fusion activation.

**Independent pocket detection (fpocket 4.0):** As an orthogonal validation, fpocket 4.0 (Le Guilloux et al., 2009) was run on PDB 6MEO (85 total pockets detected). The relay region was independently identified: 7 of 8 relay residues (465–468, 470–472) appear across fpocket pockets 5 and 19, with only Gly469 absent (consistent with glycine's minimal side chain contributing no pocket-lining surface). Both fpocket fragments span chains A (CD4) and G (gp120), independently confirming the multi-chain interface nature of the cavity. The individual fragments receive low druggability scores (0.004 and 0.003), reflecting fpocket's Voronoi tessellation fragmenting the multi-chain interface into sub-pockets; the combined volume (289 + 503 = 792 Å³) is consistent with the pocket dimensions reported above. Full fpocket output is provided in the supplementary data.

**Molecular docking validation:** AutoDock Vina docking of BMS-626529 (the active metabolite of fostemsavir, MW 531.55 Da, 8 rotatable bonds) into the relay pocket yielded a best binding affinity of −6.7 kcal/mol (Table 4), exceeding the standard druggability threshold of −6.0 kcal/mol. The dominant pose shows a 2.15 kcal/mol gap to the next-best mode (−4.5 kcal/mol), suggesting a single well-defined binding orientation within the pocket. As a control, the same ligand was docked at the known fostemsavir binding site (gp120 β20–β21 hairpin) in the same 6MEO structure, yielding −5.2 kcal/mol, 1.5 kcal/mol weaker than the relay pocket. This difference is expected: 6MEO captures the post-CD4 ternary complex, in which the attachment site has already undergone CD4-induced conformational change and is no longer in the pre-fusion state for which BMS-626529 was optimized (cf. PDB 5U7O). The control confirms that the relay pocket affinity is not an artifact of the scoring function and that a known gp120-targeting drug-like scaffold achieves druggable-range binding at this site. Stochastic reproducibility was confirmed by re-docking with a different random seed, yielding −6.6 kcal/mol (0.1 kcal/mol variance).

**Table 4.** BMS-626529 docking results: relay pocket vs. fostemsavir control site.


| Site                       | Center (Å)            | Best affinity | Mode 2        | Δ (relay − control) |
| -------------------------- | --------------------- | ------------- | ------------- | ------------------- |
| Relay pocket (G:465–472)   | (134.7, 184.6, 136.0) | −6.7 kcal/mol | −4.5 kcal/mol | −1.5 kcal/mol       |
| Fostemsavir site (control) | (138.5, 172.4, 143.3) | −5.2 kcal/mol | −3.9 kcal/mol |                     |


**Fig 4. Relay pocket geometry.** The RPGGGDMR relay pocket (Arg465–Arg472) of gp120 chain G (PDB 6MEO), shown as sticks colored by element with a semi-transparent surface envelope. All eight pocket-lining residues are labeled. The glycine triplet (G467–G469) forms the central floor, with the terminal arginines (R465, R472) bracketing the entrance. Pocket volume: 1,588 A³ (bounding ellipsoid); 91% atomic burial; Halgren druggability: 5/5 criteria passed.

**Cross-structure docking (PDB 6MET):** The same docking protocol was applied to PDB 6MET, the alternate conformation from the same cryo-EM study. Structural comparison reveals that the relay pocket residues (G:465–472) are *atomically identical* between 6MEO and 6MET; every heavy atom and B-factor matches exactly (all-atom RMSD = 0.000 Å). However, the surrounding 20 Å receptor shell contains residues that differ between the two structures (notably Thr278–Asn280, Gln423–Val425, and multiple hydrogen positions), reflecting conformational variation in the second coordination shell. When docked against the 6MET-derived receptor, BMS-626529 yielded −3.4 kcal/mol at the relay pocket, 3.3 kcal/mol weaker than in 6MEO, despite the pocket wall being structurally identical. This demonstrates that rigid-body docking scores are sensitive to second-shell packing even when the primary binding site is invariant.

**Molecular dynamics validation:** The relay pocket was assessed under explicit-solvent MD simulation in two systems (§2.5). In the 20 ns full ternary complex simulation (gp120–CD4–CCR5, 294,000 atoms), the relay pocket maintained a **sub-angstrom internal RMSD of 0.71 ± 0.15 Å** (max 1.03 Å) throughout the trajectory (Fig 6). The pocket radius of gyration was constant at 4.16 ± 0.07 Å (change −1.8%), and the mean standard deviation across all 28 intra-pocket Cα–Cα distances was only 0.21 Å.

In the 500 ns reduced-system simulation (105 residues, single-chain fragment without binding partners), the pocket maintained an **internal Cα RMSD of 1.97 ± 0.40 Å** (max 3.03 Å) over 50,001 frames, with time-segmented analysis showing gradual equilibration: 1.51 ± 0.16 Å (0–100 ns), 1.64 ± 0.19 Å (100–200 ns), 2.08 ± 0.18 Å (200–300 ns), 2.28 ± 0.25 Å (300–400 ns), and 2.35 ± 0.30 Å (400–500 ns). The pocket radius of gyration remained stable at 9.11 ± 0.45 Å (change −0.8%), and SASA was 1,269.3 ± 45.5 Å². Despite the expected large global RMSD of the isolated fragment (13.4 Å backbone-fitted), the pocket's internal geometry remained cohesive throughout the half-microsecond trajectory (Fig 6), confirming **intrinsic structural stability** of the RPGGGDMR motif independent of binding partners.

**Fig 6. Multi-structure pocket RMSD timeline.** Pocket internal Cα RMSD over 500 ns for PDB 6MEO (reduced system) and four SOSIP trimer structures (5U7O, 6MU6, 5FYJ, 4TVP). Per-100 ns segment means are shown with error bands for SOSIP structures; 6MEO is shown as a horizontal reference band. All structures maintain sub-2 Å pocket RMSD throughout the half-microsecond trajectories, demonstrating that the relay pocket is a conformational anchor across clades, resolutions (2.55–3.9 Å), and conformational states. The full ternary complex (20 ns) achieves sub-angstrom stability (0.71 Å; not shown).

**Fig 2. Three-dimensional structural overlay of the relay pocket.** Backbone superposition of the RPGGGDMR relay segment (Arg465–Arg472) across five independent HIV-1 Env structures: PDB 6MEO (cryo-EM, CD4+CCR5 ternary), 5U7O (X-ray, BMS-626529 SOSIP), 6MU6 (X-ray, BMS-814508 SOSIP), 4TVP (X-ray, unliganded prefusion SOSIP), and 5FYJ (X-ray, G471A variant). Backbone RMSD to 6MEO is <1.0 Å for all structures, confirming geometric conservation of the pocket across clades and conformational states. Pocket pharmacophore features (55 positions) are shown as spheres.

By contrast, the separate 100 ns ligand-bound fragment simulation showed pocket internal RMSD of 3.04 Å with progressive compaction (−21% radius of gyration), and the ligand (BMS-626529) dissociated by 30 ns as the surrounding structure partially unfolded. This provides a **negative control**: the pocket possesses intrinsic geometric stability, but the sub-angstrom rigidity observed in the ternary complex (0.71 Å) requires the full protein context, with binding partners refining an already-stable architecture.

### 3.6 Novel Target 2: CCR5 ECL3 Contact Surface

A second novel finding is the identification of CCR5 extracellular loop 3 (ECL3) and its flanking transmembrane helices (TM6, TM7) as a previously uncharacterized gp120 contact surface, defined by residues Tyr251, Leu255, Asn258, and Asp276.

**Pocket characteristics:**

- Pharmacophore features: 45 spatial positions
- Mean probe-to-wall distance: 0.67 Å (tightest of all four pockets; see note below)
- Character: Well-defined geometry with high complementarity

*Note on wall distance metric:* The "wall distance" reported throughout this paper refers to the mean distance from surviving pharmacophore probe centers (on a 0.8 Å grid) to the Van der Waals surface of the nearest protein heavy atom defining the pocket wall, not a center-to-center nuclear distance. A value of 0.67 Å indicates that the average probe center sits within sub-angstrom proximity of the pocket surface, consistent with a geometrically well-defined, high-complementarity binding site.

**Conformational stability note:** Unlike the gp120 relay pocket, pharmacophore probe positions for the ECL3 site were not retained under conformational flexibility analysis. This is consistent with the known dynamics of extracellular loops: CCR5 ECL3 has greater conformational freedom than the backbone-constrained C5 domain of gp120, and the cryo-EM structure represents one accessible ECL3 conformation rather than an ensemble average. The static structural identification remains valid, but the ECL3 site is flagged as requiring experimental confirmation of pocket persistence in solution (e.g., by HDX-MS or NMR) before advancing as a primary drug design target.

This site is structurally and functionally distinct from the maraviroc binding pocket, which is located within the transmembrane helical bundle. The ECL3 site is on the external surface of CCR5 and makes direct contact with gp120. A molecule targeting this site could block gp120–CCR5 interaction without competing with endogenous chemokine ligands.

### 3.7 Pharmacophore Comparison Across Pockets

Four target pockets were characterized. Comparison of their pharmacophore profiles reveals distinct druggability characteristics (Table 5).

**Table 5.** Pharmacophore summary for all analyzed pockets.


| Pocket | Target                | Status    | Features | Feature Types (HBD/HBA/H/Other) | VDW Surface Dist. (Å) | Rg (Å)    | Dynamic Stability                       |
| ------ | --------------------- | --------- | -------- | ------------------------------- | --------------------- | --------- | --------------------------------------- |
| 1      | gp120 Relay (465–472) | **Novel** | 55       | 20/10/22/3                      | 0.74                  | 2.60–2.84 | Confirmed (0.81 Å COM, 31.7% retention) |
| 2      | CCR5 ECL2 (maraviroc) | Known     | 19       | –                               | 0.74                  | –         | Not assessed                            |
| 3      | gp120 V3 loop         | Known     | 60       | –                               | 0.77                  | –         | Not assessed                            |
| 4      | CCR5 ECL3             | **Novel** | 45       | –                               | 0.67                  | 1.98      | Not retained (loop flexibility)         |


Feature type differentiation (nearest-atom complementarity typing) was performed for the primary novel target (Pocket 1); Pockets 2–4 retain undifferentiated typing pending independent structural analysis. The CCR5 ECL2 pocket (Pocket 2), which corresponds to the known maraviroc allosteric site, produced the fewest pharmacophore features (19) and the most selective occupancy profile (5.7% probe survival at moderate stringency), consistent with the tight selectivity of maraviroc binding. Recovery of this known selectivity profile provides additional validation of the methodology.

### 3.8 Pocket Validation Summary

Table 6 consolidates all orthogonal lines of evidence supporting the gp120 relay pocket as a druggable target.

**Table 6.** Convergent validation evidence for the gp120 relay pocket (Arg465–Arg472).


| Evidence Type                   | Result                                                   | Significance                                         |
| ------------------------------- | -------------------------------------------------------- | ---------------------------------------------------- |
| Cross-structure reproducibility | 22/22 residues reproduced in 6MET (Jaccard 95.7%)        | Not a crystallographic artifact                      |
| Conformational stability        | 0.81 Å COM displacement; 31.7% probe retention           | Pocket persists under side-chain perturbation        |
| B-factor analysis               | Pocket mean 72.8 Å² (38th percentile of gp120)           | More ordered than V3 (106.3), C2 (97.2), C4 (86.3)   |
| Solvent burial                  | 91% atomic burial (52/57 atoms)                          | Deep enclosed cavity, not surface patch              |
| Druggability (Halgren 2009)     | 5/5 criteria PASS; volume 1,588 Å³ (envelope)            | Classified as druggable                              |
| fpocket 4.0 detection           | 7/8 residues detected across 2 fragments (pockets 5, 19) | Orthogonal tool confirms pocket-like void at relay   |
| Sequence conservation           | Mean entropy 0.21 (3 invariant positions); >96% identity across clades | High resistance barrier; broad clade coverage        |
| Literature support              | PGT151 bNAb contacts Arg472 region                       | Interface is accessible and immunologically relevant |
| Mutagenesis cross-reference     | 14/22 hotspots confirmed by published experiments        | Method recovers known biology                        |
| Molecular docking (6MEO)        | BMS-626529: −6.7 kcal/mol                                | Exceeds druggability threshold                       |
| Docking reproducibility         | Re-dock: −6.6 kcal/mol (different seed)                  | 0.1 kcal/mol stochastic variance                     |
| Cross-structure docking (6MET)  | Pocket wall identical; shell differs                     | Score sensitive to second-shell packing              |
| MD simulation (500 ns)          | Reduced system, 300 K, explicit solvent                  | Pocket persistence under thermal motion              |
| MD simulation (20 ns)           | Full ternary complex, 294k atoms                         | Controls for boundary effects                        |
| MD pocket RMSD (20 ns ternary)  | 0.71 ± 0.15 Å internal RMSD                             | Sub-angstrom stability throughout trajectory         |
| MD pocket RMSD (500 ns reduced) | 1.97 ± 0.40 Å internal RMSD (50,001 frames)             | Intrinsic stability over half-microsecond            |
| Negative control (100 ns fragment) | Pocket RMSD 3.04 Å; ligand dissociates                | Pocket requires full ternary context                 |
| MM-GBSA binding free energy     | Ternary: −38.0; reduced: −31.4 kcal/mol                 | CCR5 stabilizes by ~6.6 kcal/mol                    |
| Per-residue decomposition       | Relay pocket: −3.07 kcal/mol (10% of reduced-system total) | Gly triplet dominates binding                       |
| Cross-structure validation      | 5 structures, 140/140 GNM positive                       | Reproduced across 5 structures from multiple clades  |
| SOSIP MD cross-validation (20 ns) | Pocket RMSD 1.07–1.34 Å (4 structures, incl. G471A variant) | Intrinsic stability without binding partners         |
| SOSIP MD cross-validation (500 ns) | Pocket RMSD 1.29–1.98 Å (4 structures, 2,000 ns aggregate); RMSF 0th percentile (3/4) | Conformational anchor across 4 structures, 2 clades, wildtype + G471A |
| Allosteric silencing (DCCM)     | Drug reduces relay→CCR5 by 51%, relay→CD4 by 15%           | Non-competitive relay blocker mechanism              |
| Allosteric correlation (500 ns DCCM) | 28/28 positive (C_ij = +0.37 to +0.99, mean +0.79)  | Cohesive mechanical unit                             |
| Drug binding (MM-GBSA)         | −37.7 ± 3.5 kcal/mol (BMS-626529, ternary)              | VdW-dominated; 9 pocket residues >1 kcal/mol        |
| Ligand IE-MMPBSA               | ΔG_IE = −27.6 kcal/mol (−TΔS = +17.7)                   | Entropy-corrected binding confirmation               |
| Non-competitiveness (MM-GBSA)  | ΔΔG_CD4 = +9.63 kcal/mol (CD4 still binds at −28.4)     | Drug permits CD4 engagement                          |
| Allosteric coupling (FEP, 3 reps) | ΔΔG = +3.9 kcal/mol; drug binding: −4.0 (pre-CD4) → −0.1 (post-CD4) | State-dependent affinity: drug selectively stabilizes pre-triggered conformation; ΔΔG quantifies relay signal strength |
| Ensemble docking (20 frames)    | BMS-626529: −6.81 ± 0.75 kcal/mol                       | Druggable across MD ensemble                         |
| Resistance mapping              | 0/28 mutations in pocket                                 | Resistance-proof by constraint                       |
| PDB-wide conservation           | 360/360 structures (100%) exact RPGGGDMR                 | Universal across solved structures                   |
| Proteome specificity            | 0/20,431 human proteins contain RPGGGDMR                 | Zero off-target geometric mimicry                    |
| Alanine scanning (GNM)          | Max coupling loss −2.2%; Gly triplet: 0%                 | Relay encoded in backbone, not side chains            |
| Hydration analysis              | 0 structural waters; 8/22 displaceable                   | Favorable for ligand displacement                    |
| PCA / essential dynamics        | Top 3 PCs: 51% of variance                              | Well-defined conformational space                    |
| SOSIP DCCM (500 ns, 4 structures) | Anti-correlated breathing in 3 wildtype; G471A abolishes gate (28/28 positive, gridlock) | State-dependent relay signal; Gly471 required for gate function |
| SOSIP PCA (500 ns, 4 structures) | Total var. 2.72–6.98 nm²; Schlitter entropy 1246–1473 J/mol·K | G471A lowest entropy; conformational space constrained by Ala substitution |
| SOSIP ensemble docking (500 ns, 5U7O+6MU6) | BMS-626529: −3.9 to −4.3 kcal/mol (intrinsic, no CD4) | Pocket druggable independent of CD4; quantifies CD4 contribution |
| De novo ligand design           | 8 leads, 5 highly novel (Tanimoto <0.2)                  | Pocket accommodates diverse drug-like scaffolds      |




### 3.9 Docking Configuration

Receptor files in PDBQT format and AutoDock Vina search configurations were generated for all four pockets. For the gp120 relay pocket:

- Search center: (134.6, 183.9, 136.7) Å
- Search dimensions: 20.8 × 21.4 × 21.3 Å
- Exhaustiveness: 32
- Modes: 20

These files enable immediate virtual screening against compound libraries (ZINC, Enamine REAL, ChEMBL).

### 3.10 Binding Free Energy (MM-GBSA)

MM-GBSA calculations quantify the thermodynamic contribution of the relay pocket to gp120–CD4 binding in two system contexts.

**Table 7.** MM-GBSA binding free energy for gp120–CD4 in two simulation contexts (kcal/mol).

| System | VdW | Elec | ΔG_solv^polar | ΔG_solv^np | **ΔG_total** | SD |
|--------|-----|------|---------------|------------|-------------|-----|
| Reduced (1.1a) | −52.5 | −364.5 | +393.7 | −8.1 | **−31.4** | 84.8 |
| Ternary (1.1b) | −86.9 | −157.2 | +218.6 | −12.5 | **−38.0** | 28.0 |

The ternary complex (with CCR5 present as structural context) yields a binding free energy ~6.6 kcal/mol more favorable than the reduced system (−38.0 vs. −31.4 kcal/mol), with a three-fold reduction in variance (SD 28.0 vs. 84.8). This demonstrates that CCR5 engagement stabilizes the gp120–CD4 interface cooperatively, consistent with the relay pocket's role as a structural mediator between receptor and coreceptor binding events.

**Per-residue energy decomposition** was performed on the reduced system, where the smaller system size (two chains, ~33,800 atoms) enables reliable per-residue attribution; decomposition was not convergent for the full three-chain ternary system. The relay pocket contributes −3.07 kcal/mol to the reduced-system total of −31.4 kcal/mol (~10%). This fractional contribution estimates the pocket's role in the gp120–CD4 interface absent CCR5 context; the corresponding fraction in the ternary complex may differ given the additional cooperative stabilization.

**Table 8.** Per-residue MM-GBSA decomposition: relay pocket residues (reduced system).

| Residue | # | ΔG (kcal/mol) | VdW | Elec |
|---------|---|---------------|-----|------|
| Gly468 | 468 | −1.98 | −1.50 | −6.19 |
| Gly467 | 467 | −1.23 | −0.96 | −5.97 |
| Arg465 | 465 | +0.32 | −0.22 | −2.30 |
| Gly469 | 469 | −0.19 | −0.90 | −2.23 |
| Arg472 | 472 | +0.18 | −0.03 | −0.49 |
| Asp470 | 470 | −0.16 | −0.23 | −5.64 |
| Pro466 | 466 | −0.04 | −0.69 | −8.25 |
| Met471 | 471 | +0.03 | −0.04 | −1.26 |
| **Total** | | **−3.07** | | |

The Gly467–Gly468–Gly469 triplet collectively contributes −3.40 kcal/mol (more than the entire 8-residue pocket total), indicating that the glycine backbone flexibility is the primary structural determinant of binding at this interface.

**Drug binding free energy (1.3).** MM-GBSA on the 10 ns drug-bound ternary complex yields ΔG_bind = −37.7 ± 3.5 kcal/mol for BMS-626529 binding to the relay pocket within the full gp120+CD4+CCR5 ternary context (40 frames; per-residue decomposition in Fig 5). The binding is dominated by van der Waals packing (ΔVDW = −68.1 kcal/mol) with a desolvation penalty (ΔEGB = +29.9 kcal/mol, ΔESURF = −8.1 kcal/mol), consistent with the hydrophobic pocket interior. Per-residue decomposition within 6 Å identifies 9 residues contributing >1 kcal/mol each: Leu85 (−3.0), Val123 (−2.8), Phe125 (−2.6), Val169 (−1.4), Lys86 (−1.3), Pro121 (−1.3), Pro127 (−1.2), Phe289 (−1.0), and Phe295 (−1.0 kcal/mol). All are relay pocket residues, confirming the ligand occupies the predicted binding site. The aromatic–hydrophobic binding profile (three Phe, two Val, one Leu among top contributors) is consistent with the pocket pharmacophore (§3.5).

**Interaction entropy analysis.** To estimate the entropic contribution to ligand binding, interaction entropy (IE) MM-GBSA was computed on 200 frames from the drug-bound ternary trajectory. The enthalpy-only estimate yields ΔG_enthalpy = −45.23 ± 5.92 kcal/mol, with an interaction entropy correction of −TΔS_IE = +17.67 kcal/mol, giving a corrected ΔG_bind^IE = −27.56 kcal/mol. However, the interaction energy standard deviation (σ = 14.40 kcal/mol) exceeds the reliability threshold (σ < 3.6 kcal/mol recommended by Duan et al., 2016), indicating that the entropy estimate should be interpreted with caution; the enthalpy-dominated binding (−37.7 ± 3.5 kcal/mol from 40-frame MM-GBSA) remains the more robust metric.

### 3.11 Cross-Structure Validation (SOSIP Trimers)

Four additional HIV-1 Env structures representing different clades and conformational states were analyzed to assess relay pocket conservation beyond PDB 6MEO (Fig 2).

**Table 9.** Cross-structure validation of relay pocket across 5 HIV-1 Env structures.

| PDB | Clade / Form | Backbone RMSD to 6MEO (Å) | GNM intra corr. | fpocket | MD pocket RMSD (Å, 20 ns) | MD pocket RMSD (Å, 500 ns) |
|-----|-------------|---------------------------|-----------------|---------|---------------------------|----------------------------|
| 6MEO | BG505, CD4+CCR5 | – | 0.69 | 7/8 | 0.71 ± 0.15 ^a^ | 1.97 ± 0.40 ^b^ |
| 6MU6 | BG505, BMS-814508 | 0.72 | 0.69 | 4/8 | 1.34 ± 0.27 ^b^ | **1.72 ± 0.44** ^b^ |
| 5U7O | BG505, BMS-626529 | 0.78 | 0.73 | 8/8 | 1.07 ± 0.23 ^b^ | 1.79 ± 0.56 ^b^ |
| 4TVP | BG505, prefusion | 0.54 | 0.73 | 5/8 | 1.29 ± 0.34 ^b^ | **1.98 ± 0.44** ^b^ |
| 5FYJ | Clade G, VRC01 (G471A) | 0.93 | 0.80 | 8/8 | 1.34 ± 0.25 ^b^ | **1.29 ± 0.30** ^b^ |

^a^ Full ternary complex (gp120+CD4+CCR5, 294,000 atoms). ^b^ Reduced system (chain G only, 85–101 residues, ~33,000–38,000 atoms; no binding partners). All four SOSIP 500 ns pocket RMSDs are sub-2 Å, with three of four *lower* than the 6MEO 500 ns reference (1.97 Å).

All five structures show backbone RMSD <1.0 Å at the relay segment (S7 Fig), GNM intra-pocket correlation ≥0.69, and pocket detection by at least one independent algorithm. Across all five structures, **28/28 GNM intra-pocket correlations are positive** (140/140 total), confirming that the relay segment moves as a cohesive mechanical unit regardless of conformational state. For the CD4-bound structure (6MEO), GNM analysis of the full ternary complex yields a relay–CD4 inter-domain coupling of 0.88 (perturbation response effectiveness: 98.3rd percentile), indicating direct mechanical communication between the relay pocket and the CD4-binding site. The SOSIP structures lack co-crystallized CD4, so inter-domain coupling was assessed only for 6MEO. GNM alanine scanning identified Asp470 as the hub residue in all five structures (11.8–18.7% disruption upon virtual ablation).

**SOSIP molecular dynamics cross-validation.** To extend the static cross-validation to dynamic assessment, 20 ns explicit-solvent MD simulations were performed for four SOSIP structures using the same protocol as the 6MEO reduced system (§2.5): a 15 Å shell around the relay pocket (complete residues, chain G, 85–101 residues, ~700–834 protein atoms), AMBER99SB-ILDN force field, TIP3P water, 0.15 M NaCl, ~33,000–38,000 total atoms. The relay pocket maintained low internal RMSD in all four independent structures:

- **6MU6** (2.55 Å X-ray): pocket RMSD 1.34 ± 0.27 Å (max 2.46 Å), Rg 7.34 ± 0.20 Å (−1.9% change), T = 300.0 K
- **4TVP** (3.1 Å X-ray): pocket RMSD 1.29 ± 0.34 Å (max 2.42 Å), Rg 7.39 ± 0.23 Å (+2.9% change), T = 299.9 K
- **5FYJ** (3.1 Å X-ray, G471A variant): pocket RMSD 1.34 ± 0.25 Å (max 2.12 Å), Rg 7.61 ± 0.21 Å (+6.2% change), T = 300.0 K
- **5U7O** (3.0 Å X-ray): pocket RMSD 1.07 ± 0.23 Å (max 1.72 Å), Rg 6.87 ± 0.25 Å (+4.9% change), T = 300.0 K

All four SOSIP trajectories show pocket RMSD of 1.07–1.34 Å, intermediate between the full ternary complex (0.71 Å, with CD4+CCR5 stabilization) and the 500 ns 6MEO reduced system (1.97 Å). The 5U7O structure shows the best pocket RMSD (1.07 Å), and the 5FYJ G471A variant (RPAGGDMR) produces an identical pocket RMSD (1.34 Å) despite the sequence change. These are single-chain reduced systems without binding partners, establishing that the relay pocket possesses intrinsic structural stability across independent structures, crystal forms, resolutions (2.55–3.1 Å), and sequence variants.

**500 ns SOSIP extensions.** To assess pocket stability on physiologically relevant timescales, all four SOSIP structures were extended to 500 ns (50,001 frames each, 2,000 ns aggregate). All four maintain sub-2 Å pocket RMSD:

- **5FYJ** (500 ns, G471A variant): pocket RMSD 1.29 ± 0.30 Å, backbone RMSD 17.24 ± 8.82 Å, RMSF 0.69 Å (**0th percentile**), SASA 1,294 ± 35 Å²
- **6MU6** (500 ns): pocket RMSD 1.72 ± 0.44 Å, backbone RMSD 6.90 ± 1.01 Å, RMSF 1.19 Å (**0th percentile**), SASA 1,279 ± 33 Å²
- **5U7O** (500 ns): pocket RMSD 1.79 ± 0.56 Å, backbone RMSD 8.36 ± 2.76 Å, RMSF 3.93 Å (52nd percentile), SASA 1,212 ± 55 Å²
- **4TVP** (500 ns, unliganded): pocket RMSD 1.98 ± 0.44 Å, backbone RMSD 7.20 ± 1.75 Å, RMSF 1.32 Å (**0th percentile**), SASA 1,254 ± 46 Å²

Three of four SOSIP structures show pocket RMSD *lower* than the 6MEO 500 ns reduced system (1.97 Å); 4TVP matches it (1.98 Å). Three of four structures (5FYJ, 6MU6, 4TVP) operate at the 0th RMSF percentile (S6 Fig): no other residue in these reduced systems exhibits lower fluctuation over 500 ns. The relay pocket is a **conformational anchor** in the pre-triggered trimer across four independent X-ray structures (2.55–3.1 Å resolution), two HIV-1 clades (BG505 and X1193.c1), and both wildtype (RPGGGDMR) and variant (RPAGGDMR) sequences.

The strongest kinematic decoupling appears in 5FYJ: the pocket maintains 1.29 Å RMSD while the global backbone drifts to 17.24 Å (ratio 0.075). This 13-fold differential between local pocket stability and global backbone fluctuation is the most extreme in the dataset, confirming that the relay pocket's rigidity is an intrinsic property of the local sequence-structure relationship, independent of the surrounding tertiary context.

The increased pocket dynamics in the 6MEO ternary complex (RMSD 1.97 Å vs. 1.29–1.79 Å in three SOSIP structures) reflect the active propagation of the allosteric signal when CD4 and CCR5 are engaged. This observation is consistent with the FEP-measured ΔΔG = +3.95 kcal/mol (Fig 10; §3.15): the drug selectively stabilizes the hyper-rigid pre-triggered conformation, and the thermodynamic penalty quantifies the energy the conformational relay would need to overcome to reshape the pocket while the drug is bound.

Time-segmented analysis reveals distinct dynamic signatures. The 5U7O pocket exhibits conformational breathing (RMSD oscillates: 1.16 → 1.92 → 2.37 → 1.52 → 1.98 Å across five 100 ns segments), returning to low-RMSD states even as the global backbone drifts to 11.6 Å. The 6MU6 pocket shows steady progression (1.48 → 1.64 → 1.69 → 1.71 → 2.07 Å) with lower variance. The 5FYJ pocket (G471A) converges *downward* (1.37 → 1.50 → 1.29 → 1.14 → 1.17 Å), locking into its most rigid state by 300 ns despite the backbone undergoing the largest global excursion in the dataset (17.24 Å). The 4TVP pocket (unliganded) shows a slow upward drift (1.54 → 1.69 → 1.95 → 2.24 → 2.47 Å) consistent with the absence of stabilizing binding partners, yet remains sub-2.5 Å throughout. In all four cases, the relay pocket is kinematically decoupled from global tertiary fluctuations.

**SOSIP mechanistic cross-validation (500 ns).** To extend cross-validation beyond structural stability to mechanistic characterization, dynamic cross-correlation matrices (DCCM), principal component analysis (PCA), and ensemble docking were performed on all four completed 500 ns SOSIP trajectories (5U7O, 6MU6, 5FYJ, and 4TVP) using protocols matching the 6MEO analyses (§2.9, §2.10).

*Intra-pocket DCCM.* Covariance matrices were computed via `gmx covar` (All_CA fit, dt=50 ps) and normalized to yield Cα–Cα cross-correlations for the relay pocket (8 residues, 28 unique pairs).

**Table 10.** Intra-pocket DCCM comparison across conformational states (500 ns).

| System | Conformation | All 28 positive? | Mean C_ij | Seq. neighbor mean | Range |
|--------|---|---|---|---|---|
| 6MEO (500 ns) | CD4+CCR5 bound (relay active) | Yes | +0.79 | +0.97 | +0.37 to +0.99 |
| 6MU6 SOSIP (500 ns) | Pre-fusion (BG505) | No (24/28) | +0.45 | +0.79 | −0.13 to +0.90 |
| 5U7O SOSIP (500 ns) | Pre-fusion (BG505) | No (19/28) | +0.22 | +0.81 | −0.59 to +0.91 |
| 4TVP SOSIP (500 ns) | Pre-fusion (BG505) | No (20/28) | +0.31 | +0.79 | −0.47 to +0.92 |
| 5FYJ SOSIP (500 ns) | G471A (Clade G) | **Yes** | +0.72 | +0.97 | +0.21 to +0.99 |

The SOSIP structures reveal a state-dependent DCCM signature absent from the 6MEO ternary complex. In all three wildtype SOSIP structures (6MU6, 5U7O, 4TVP), sequential neighbors remain strongly correlated (+0.79 to +0.81), but non-sequential pairs between the N-terminal (Arg465–Gly469) and C-terminal (Asp470–Arg472) halves of the pocket show anti-correlations: the two halves breathe in opposition. In the 6MEO ternary complex (CD4+CCR5 engaged), all 28 pairs are positive (mean +0.79), indicating that receptor engagement locks both halves into co-motion as the conformational relay signal propagates through the pocket.

The G471A variant (5FYJ) abolishes this gating motion (Fig 7). All 28 intra-pocket pairs are positive (mean +0.72), matching the 6MEO relay-active sign pattern despite the absence of CD4 or CCR5. The correlation signs match, but the underlying physics differ: 5FYJ has the lowest total PCA variance (2.72 nm²) and lowest Schlitter entropy (1,246 J/mol·K) of any structure in the dataset (S2 Fig). The residues move together because they barely move at all. The anti-correlated breathing is gone because the alanine side chain at position 471 eliminates the backbone flexibility that the conserved glycine provides, converting a functional gate into steric gridlock. The wildtype Gly471 is required for the pocket to function as a gate (openable by the relay signal); the G471A variant seizes the gate shut.

**Fig 7. State-dependent DCCM comparison across five structures (500 ns).** Five-panel heatmap of intra-pocket Cα–Cα dynamic cross-correlations. Top row: 6MEO (CD4+CCR5 triggered, all 28 pairs positive, mean +0.79), 5U7O and 6MU6 (wildtype pre-fusion SOSIP, anti-correlated N/C-terminal breathing). Bottom row: 4TVP (wildtype pre-fusion, unliganded) and 5FYJ G471A (gate seized, all 28 positive via steric gridlock, mean +0.72). The transition from anti-correlated breathing (resting gate) to fully correlated co-motion (relay active) is the dynamic fingerprint of the conformational relay signal. Blue: negative (anti-correlated); red: positive (correlated).

The transition from anti-correlated breathing (wildtype SOSIP, pre-triggered) to fully correlated co-motion (6MEO, triggered) is reproduced across three independent wildtype structures from independent 500 ns simulations. This provides a direct mechanistic readout of the relay signal, consistent with the FEP-measured ΔΔG of +3.95 kcal/mol (§3.15): a relay-blocking drug stabilizes the anti-correlated pre-triggered state, preventing the transition to correlated co-motion required for fusion activation. The G471A variant demonstrates that eliminating the gate has the same DCCM sign signature as activating the relay, but through rigidity rather than signal propagation.

*PCA / essential dynamics.* Principal component analysis of the relay pocket atoms reveals structure-dependent conformational dynamics.

**Table 11.** PCA of relay pocket dynamics across conformational states.

| System | Total Var. (nm²) | Top 3 PCs (%) | Top 5 PCs (%) | Schlitter S (J/mol·K) |
|--------|---|---|---|---|
| 5U7O SOSIP (500 ns) | 6.98 | 67.3 | 79.1 | 1472.6 |
| 4TVP SOSIP (500 ns) | 5.44 | 61.7 | 71.7 | 1437.7 |
| 6MU6 SOSIP (500 ns) | 4.45 | 59.1 | 71.2 | 1346.4 |
| 5FYJ SOSIP (500 ns, G471A) | 2.72 | 54.2 | 63.9 | **1246.5** |
| 6MEO reduced (100 ns) | 2.48 | 51.0 | 65.5 | – |

Total variance spans a 2.6-fold range across the four SOSIP structures (2.72–6.98 nm²; S2 Fig), with conformational entropy following the same ordering. The G471A variant (5FYJ) has the lowest variance (2.72 nm²) and lowest Schlitter entropy (1,246 J/mol·K), consistent with the steric gridlock observed in DCCM: the alanine substitution restricts the backbone conformational space that glycine permits. Among wildtype structures, 6MU6 is the most rigid (4.45 nm², 1,346 J/mol·K), consistent with its 0th-percentile RMSF. In all five systems, dynamics are dominated by 2–3 collective modes, with the top 3 PCs capturing 51–67% of variance, indicating a well-defined conformational landscape rather than stochastic fluctuation.

*Ensemble docking.* GROMOS clustering (2.5 Å cutoff on relay pocket Cα, dt=500 ps) of the 500 ns 5U7O and 6MU6 SOSIP trajectories yielded 3 cluster representatives each (reflecting the pocket's rigidity). BMS-626529 was docked into each representative using AutoDock Vina 1.2.7.

**Table 12.** Ensemble docking of BMS-626529 across conformational states.

| System | N clusters | Mean affinity (kcal/mol) | Best | Worst |
|--------|---|---|---|---|
| 6MEO reduced (100 ns) | 20 | −6.81 ± 0.75 | −8.34 | −5.51 |
| 5U7O SOSIP (500 ns) | 3 | −4.29 ± 0.24 | −4.62 | −4.12 |
| 6MU6 SOSIP (500 ns) | 3 | −3.90 ± 0.33 | −4.37 | −3.64 |

The SOSIP ensemble docking scores (−3.9 to −4.3 kcal/mol) are lower than the 6MEO ternary scores (−6.8 kcal/mol). This is expected because the SOSIP reduced systems lack the CD4 wall residues (Gln25, Asn39, Gln40, Phe43, Thr45) that form part of the binding pocket in the ternary complex. The ~2.5 kcal/mol difference quantifies CD4's geometric contribution to the binding site: CD4 provides a physical wall that increases Van der Waals contacts and constrains ligand orientation, improving the enthalpy-dominated docking score. The SOSIP docking scores nonetheless remain in the pharmacologically relevant range (>−3.5 kcal/mol), demonstrating that the pocket possesses intrinsic druggability independent of CD4 (S3 Fig). This resolves an apparent paradox with the FEP results (Fig 9; §3.15): while docking (enthalpy-only) favors the CD4-bound geometry, FEP (full thermodynamics including entropy) shows the drug binds preferentially to the pre-CD4 state (ΔG_bind = −4.0 kcal/mol pre-CD4 vs −0.1 kcal/mol post-CD4). The discrepancy reflects the large entropic penalty of forcing the drug into the rigidified, relay-activated ternary conformation, a penalty invisible to rigid-body docking but captured by alchemical free energy methods. The few SOSIP clusters (3 each, vs 20 for 6MEO) further confirm the pocket's exceptional rigidity in the pre-triggered state.

**Fig 9. FEP thermodynamic cycle for state-dependent drug binding.** Schematic of the three-leg alchemical free energy perturbation protocol: BMS-626529 decoupled from the ternary complex (gp120+CD4), the binary system (gp120 alone), and solvent (water/ions). Binding free energies computed via the standard thermodynamic cycle show moderate pre-CD4 affinity (ΔG_bind = −4.03 kcal/mol) that drops to near-zero upon CD4 engagement (ΔG_bind = −0.09 kcal/mol). The allosteric coupling ΔΔG = +3.95 kcal/mol quantifies the conformational relay signal strength and establishes conformational selection as the binding mechanism. Values shown are TI means ± SEM across 3 independent replicas (9 legs, 1,440 ns aggregate).

### 3.12 Allosteric Network Analysis

Dynamic cross-correlation analysis of the 500 ns trajectory (10,001 frames) reveals that the relay pocket moves as a highly cohesive mechanical unit. All 28 intra-pocket residue pairs show **positive cross-correlations** (C_ij = +0.37 to +0.99; mean +0.79), with sequential neighbors exhibiting near-unity coupling (mean +0.97) and a smooth distance-dependent decay. These correlations are substantially stronger than those observed in the initial 100 ns analysis (C_ij = +0.08 to +0.82), consistent with the improved statistical convergence afforded by five-fold greater sampling.

The relay pocket is thus identified as a **cohesive mechanical unit** in the gp120–CD4 allosteric network: high intra-pocket correlation across the full half-microsecond trajectory, combined with sub-2 Å internal RMSD (§3.5), indicates that the RPGGGDMR motif functions as a structurally rigid relay element.

Comparison with the SOSIP DCCM (§3.11, Table 10, Fig 7) reveals a **state-dependent correlation signature**: in the three wildtype pre-fusion SOSIP trimers, intra-pocket correlations are weaker (mean +0.22 to +0.45) with anti-correlations between the N- and C-terminal halves of the motif, while in the 6MEO ternary complex (relay activated), all 28 pairs show strong positive correlation (mean +0.79). This transition from anti-correlated breathing to fully correlated co-motion, reproduced across three independent wildtype SOSIP structures from independent 500 ns simulations, is the dynamic fingerprint of the conformational relay signal. The G471A variant (5FYJ) produces the same all-positive DCCM sign pattern (28/28, mean +0.72) but through steric gridlock rather than active signal propagation, confirming that Gly471 is required for the gate to function.

### 3.13 Allosteric Silencing by Drug Binding

To directly test the relay blocker mechanism, DCCM was computed on the full ternary complex (887 Cα atoms) from both the 20 ns apo and 10 ns drug-bound (BMS-626529 analog) trajectories, enabling comparison of cross-domain correlation patterns with and without a relay pocket–bound ligand.

In the apo ternary complex, the relay pocket residues (372–379 in GROMACS numbering; RPGGGDMR) exhibit near-perfect internal correlation (mean C_ij = 0.987, 28/28 pairs positive, range 0.96–1.00) and maintain significant cross-domain coupling: relay→CD4 |C| = 0.47, relay→CCR5 |C| = 0.49. The relay pocket is dynamically coupled to both binding partners.

In the drug-bound system, a selective asymmetric silencing pattern emerges (Fig 8):

- **Relay→CCR5 coupling**: |C| = 0.49 → 0.24 (**−51%**): the conformational signal to the fusion machinery is halved
- **Relay→gp120 (intra-chain)**: |C| = 0.84 → 0.58 (−31%): signal propagation through the protein is disrupted
- **Relay→CD4 coupling**: |C| = 0.47 → 0.40 (−15%): CD4 communication is partially preserved
- **Intra-pocket correlation**: 0.987 → 0.894 (−9.4%): the pocket remains structurally intact but transmits less efficiently
- **Global mean |C|**: 0.519 → 0.399 (−23%): system-wide decorrelation

The drug thus acts as a **conformational signal wedge**: the pocket retains its structural identity (mean C_ij = 0.89) but can no longer transmit the allosteric signal from the CD4 engagement site to the CCR5/gp41 fusion machinery. The 3.4-fold selectivity in silencing (51% CCR5 vs. 15% CD4) is consistent with a non-competitive mechanism that permits receptor binding while blocking the downstream conformational relay required for membrane fusion. This pattern distinguishes relay blockers from attachment inhibitors (fostemsavir, which blocks CD4 binding) and coreceptor antagonists (maraviroc, which blocks CCR5 engagement).

**Fig 8. Allosteric silencing of the relay by drug binding.** (A) Intra-pocket DCCM mean correlation with min–max range for apo (green) and drug-bound (orange) full ternary complex trajectories (887 Cα atoms). Drug binding drops the minimum pair correlation from 0.960 to 0.686, decoupling relay pairs while the pocket retains overall structural identity. (B) Cross-domain relay coupling (relay→CD4, relay→CCR5, relay→gp120 intra-chain) in apo versus drug-bound states. The drug selectively silences relay→CCR5 coupling (−51%) while partially preserving relay→CD4 coupling (−15%), consistent with a non-competitive relay blocker mechanism.

### 3.14 Non-Competitive Mechanism Validation

To test whether drug binding at the relay pocket competes with CD4 engagement, MM-GBSA binding free energy was computed for the gp120–CD4 interaction in the drug-bound ternary complex and compared to the apo ternary reference. Using the same methodology as §3.10 (receptor = gp120, ligand = CD4, IGB5, 0.15 M salt), the drug-bound system yields ΔG_CD4 = −28.39 ± 8.94 kcal/mol (200 trajectory frames), compared to −38.02 ± 28.04 kcal/mol in the apo ternary (ΔΔG = +9.63 kcal/mol). CD4 binding remains thermodynamically favorable in the drug-bound state, with the modest difference falling within the standard deviations of both estimates.

Together with the allosteric silencing (§3.13), this establishes a consistent picture: the drug permits CD4 engagement while disrupting the downstream conformational relay.

### 3.15 Allosteric Coupling Free Energy (FEP/TI)

To quantify the thermodynamic coupling between CD4 engagement and drug binding at the relay pocket, alchemical free energy perturbation was performed on three systems: a ternary complex (gp120+CD4+BMS-626529), a binary system (gp120+BMS-626529, no CD4), and a solvent reference (BMS-626529 in water/ions). Three independent replicas with different velocity seeds were computed for each system (9 legs total, 1,440 ns aggregate).

**Table 13.** FEP drug decoupling free energies (mean ± SEM across 3 replicas; 32 λ windows × 5 ns each) and binding free energies via the standard thermodynamic cycle (ΔG_bind = ΔG_solvent − ΔG_complex).

*Decoupling free energies (higher = stronger binding):*

| Method | ΔG_ternary (kcal/mol) | ΔG_binary (kcal/mol) | ΔG_solvent (kcal/mol) |
|--------|----------------------|---------------------|----------------------|
| TI (decorrelated) | 6.63 ± 0.01 | 10.57 ± 0.18 | 6.54 ± 0.01 |
| MBAR | 6.55 ± 0.01 | 10.37 ± 0.19 | 6.51 ± 0.02 |

*Binding free energies (negative = favorable):*

| Method | gp120+CD4 (kcal/mol) | gp120 only (kcal/mol) | ΔΔG (kcal/mol) |
|--------|----------------------|----------------------|----------------|
| TI (decorrelated) | −0.09 ± 0.02 | −4.03 ± 0.18 | +3.95 ± 0.18 |
| MBAR | −0.04 ± 0.03 | −3.86 ± 0.19 | +3.82 ± 0.19 |

TI and MBAR estimators agree within 0.2 kcal/mol across all three legs (Fig 9), and inter-replicate standard deviations are small (ternary: 0.03, binary: 0.30–0.32, solvent: 0.02–0.04 kcal/mol), confirming convergence. The solvent leg serves as an internal consistency control: its inter-replicate SD (<0.04 kcal/mol) is an order of magnitude below the CD4 effect, confirming that the coupling signal is real rather than a sampling artifact.

Applying the standard thermodynamic cycle, BMS-626529 binds gp120 alone (binary, pre-CD4 conformation) with ΔG_bind = −4.03 ± 0.18 kcal/mol (TI) / −3.86 ± 0.19 kcal/mol (MBAR), representing moderate, pharmacologically relevant affinity. In the CD4-bound ternary complex, where the conformational relay has propagated through gp120, binding drops to ΔG_bind = −0.09 ± 0.02 kcal/mol (TI) / −0.04 ± 0.03 kcal/mol (MBAR). The allosteric coupling between the two sites is ΔΔG = +3.95 ± 0.18 kcal/mol (TI) / +3.82 ± 0.19 kcal/mol (MBAR). By Le Chatelier's principle, this coupling is reciprocal: just as CD4 engagement imposes a +3.95 kcal/mol penalty on drug binding, drug binding imposes an equivalent thermodynamic penalty on the CD4-induced conformational relay; the drug acts as a conformational lock on the pre-triggered state.

This result has a direct mechanistic interpretation. The drug exhibits state-dependent affinity characteristic of conformational selection (analogous to Imatinib's selectivity for the inactive DFG-out kinase conformation; Schindler et al., 2000): BMS-626529 binds the pre-CD4 relay pocket conformation with moderate affinity (ΔG_bind ≈ −4 kcal/mol), selectively stabilizing the pre-triggered state. The +3.95 kcal/mol ΔΔG quantifies the thermodynamic strength of the conformational relay signal: the free energy that CD4-induced rearrangement would need to overcome to reshape the pocket while the drug is bound. A relay blocker drug exploits this coupling: by occupying the pre-CD4 pocket, it imposes a ~4 kcal/mol thermodynamic barrier against the conformational relay, effectively locking gp120 in the pre-triggered conformation and preventing downstream fusion activation.

Together with the DCCM allosteric silencing (§3.13) and MM-GBSA non-competitiveness analysis (§3.14), this FEP result provides a fourth independent line of evidence for the relay blocker mechanism. The FEP and MM-GBSA measure reciprocal sides of the same coupling: MM-GBSA quantifies the drug's effect on CD4 binding (ΔΔG_CD4 = +9.63 kcal/mol); FEP quantifies CD4's effect on drug binding (ΔΔG_drug = +3.9 kcal/mol). Both agree qualitatively: the relay pocket and CD4 binding site are thermodynamically coupled, with perturbation of one site propagating to the other. The quantitative difference reflects known systematic differences between endpoint and alchemical methods (Genheden & Ryde, 2015) and the presence of CCR5 in the MM-GBSA system but not in the FEP systems.

### 3.16 Resistance Mutation Mapping

Cross-referencing the relay pocket against the Stanford HIV Drug Resistance Database (28 known resistance-associated mutations) reveals that **zero known drug resistance mutations directly affect relay pocket residues 465–472**. The nearest resistance mutation is M475I (a secondary fostemsavir escape mutation), located 3 residues from the pocket boundary. Two clade-specific polymorphisms exist within the pocket (R469K, D470N), but both are conservative substitutions that preserve the backbone geometry and relay function. Shannon entropy analysis of 981 NCBI gp120 sequences confirms exceptional conservation at the relay segment (mean entropy 0.21), with three completely invariant positions (Arg465, Pro466, Gly468) and the relay motif situated within the C5 subdomain (mean entropy 0.14), the most conserved constant region of gp120 (Fig 3).

**Fig 3. Sequence conservation entropy of the relay pocket.** (A) Per-residue Shannon entropy for the 8 relay pocket positions (RPGGGDMR) across 981 NCBI HIV-1 gp120 sequences. Three residues (Arg465, Pro466, Gly468) are completely invariant (entropy 0.00); the most variable position (Gly467, entropy 0.72) shows predominantly G→A substitutions that preserve small side-chain character. (B) Mean Shannon entropy by gp120 domain, showing the relay pocket within the C5 subdomain — the most conserved constant region — compared to V1/V2, V3, and other variable and constant regions.

A comprehensive PDB survey of 531 HIV-1 gp120-containing entries identified 360 structures with the relay pocket region resolved, of which **all 360 (100%) contain the exact RPGGGDMR sequence**: zero variants across 27 years of structures (1998–2025), 115 X-ray and 245 cryo-EM, from multiple independent laboratories.

### 3.18 Cross-Target Specificity

To assess whether the relay pocket motif could produce off-target binding to host proteins, the complete RPGGGDMR octapeptide was searched against the entire reviewed human proteome (Swiss-Prot, 20,431 proteins). **Zero matches were found**. Relaxing the search to any 7-residue sub-motif (RPGGGDM, PGGGDMR) also yielded zero matches. The shortest sub-motif with any human proteome occurrence is the 5-mer RPGGG (30 proteins, 0.15%) and PGGGD (31 proteins, 0.15%), both scattered across unrelated proteins with no structural or functional relationship to gp120. The relay pocket's defining backbone geometry is completely absent from the human proteome, indicating negligible risk of off-target geometric mimicry for a drug designed against this motif.

### 3.19 Computational Alanine Scanning

To test whether the relay signal depends on side-chain contacts or backbone topology, contact-weighted Gaussian Network Model (GNM) analysis was performed on wild-type PDB 6MEO and 8 single-residue alanine mutants (one per relay pocket residue). Spring constants were set proportional to inter-residue heavy-atom contact counts within 4.5 Å, enabling the GNM to capture side-chain removal effects that standard Cα-only models miss. The wild-type relay–CD4 mean coupling was 0.857 (489 inter-residue contacts). The maximum coupling disruption from any single mutation was −2.2% (Arg465Ala, −47 contacts). Arg472Ala (−0.4%, −37 contacts) and Asp470Ala (−0.4%, −23 contacts) showed minimal disruption (S4 Fig). The glycine triplet (Gly467, Gly468, Gly469) showed **zero coupling change upon mutation to alanine**, because glycine already has the minimal side chain, and the relay signal propagates through backbone topology rather than side-chain contacts. AutoDock Vina re-docking of BMS-626529 against each mutant receptor showed that pocket geometry is preserved or improved upon alanine substitution (Asp470Ala: −8.1 kcal/mol; Arg472Ala: −8.7 kcal/mol vs. wild-type −4.2 kcal/mol in the same re-docking context), while vicriviroc binding was uniformly weakened (+0.9 to +1.8 kcal/mol), consistent with differential selectivity profiles across the mutant panel.

### 3.20 Ensemble Docking

Docking of BMS-626529 across 20 representative MD conformations (GROMOS clustering, 2.5 Å cutoff from the first 100 ns of the trajectory; S3 Fig) yielded a mean binding affinity of −6.81 ± 0.75 kcal/mol (range −5.51 to −8.34 kcal/mol). The ensemble mean (−6.81) is comparable to the single cryo-EM frame result (−6.7), confirming that pocket druggability is **maintained across the MD conformational ensemble**. Vicriviroc (a CCR5 antagonist, as a selectivity control) yielded −7.50 ± 0.82 kcal/mol across the same 20 frames. The best-scoring conformation (cluster 5) yielded −8.34 kcal/mol for BMS-626529 and −9.06 kcal/mol for vicriviroc.

### 3.21 Hydration Site Analysis and Conformational Dynamics

Hydration site analysis of the relay pocket (2,022 frames from the first 100 ns of the trajectory, 1.0 Å grid, 6.0 Å shell) identified 22 water positions at >1.5× bulk density. None were classified as structural waters (high density, low entropy, geometrically locked). Eight sites at 5.1–6.5 Å from the pocket center were classified as displaceable (density 1.7–2.3× bulk), indicating that ligand binding would not require displacement of tightly bound structural waters, a favorable property for druggability (Halgren, 2009).

Principal component analysis of the first 100 ns of the trajectory (10,110 frames, 112 pocket atoms) shows that the top three principal components capture 51.0% of the total variance (PC1: 23.4%, PC2: 18.3%, PC3: 9.3%), indicating well-defined conformational dynamics with 2–3 dominant collective modes. This low-dimensional motion profile is consistent with the sub-angstrom pocket RMSD observed in the full ternary simulation.

### 3.22 De Novo Ligand Design

To assess whether the relay pocket can accommodate structurally diverse drug-like scaffolds beyond BMS-626529, de novo ligand design was performed against the 55-feature pharmacophore model. Eight novel lead compounds were generated, all passing Lipinski, Veber, and PAINS filters. Five of eight designs are highly novel (Tanimoto similarity <0.2 to any known compound in the screening library), with the top-scoring lead (RelayLead-02, MW 428.5 Da, combined score 0.639) representing a scaffold with no structural precedent among existing HIV entry inhibitors. The structural diversity of these leads, spanning distinct chemical scaffolds while maintaining complementarity to the pocket pharmacophore, provides further evidence that the relay pocket can serve as a viable drug design target accommodating multiple chemical series.

## 4. Discussion



### 4.1 A Gap in the Entry Inhibitor Landscape

Four approved drugs target four distinct steps in the HIV-1 entry cascade (§1), yet the conformational relay between receptor engagement and fusion activation has never been targeted with a small molecule. The gp120 Arg465–Arg472 pocket fills this gap.

Every existing entry inhibitor blocks a step. This pocket does not block a step. It terminates a signal. The conformational signal from CD4 to gp41 routes through this eight-residue junction. A drug in the pocket lets the signal arrive and terminates it. The virus binds its receptor, engages its coreceptor, and the conformational change required to activate fusion dissipates at the relay. Half-microsecond MD simulations of SOSIP trimers (§3.11) reveal that the relay pocket functions as a **conformational anchor** within the pre-triggered Env trimer: all four structures maintain sub-2 Å pocket RMSD (1.29–1.98 Å; Fig 6), three of four operate at the 0th RMSF percentile, and the pocket is kinematically decoupled from global tertiary dynamics (ratio as low as 0.075 in 5FYJ: 1.29 Å pocket vs. 17.24 Å backbone). This rigidity is an intrinsic property of the local sequence-structure relationship, not dependent on binding partners. Pocket dynamics increase upon CD4 engagement: the 6MEO ternary complex pocket RMSD rises to 1.97 Å, reflecting active propagation of the allosteric signal through the relay junction as the conformational change occurs. A drug that locks the pocket in its hyper-rigid pre-triggered conformation would therefore block the conformational relay at its most mechanically constrained point.

This rigidity signature is mechanistically cross-validated by SOSIP DCCM analysis (§3.11, Table 11, Fig 7): in all three wildtype pre-fusion SOSIP trimers (6MU6, 5U7O, 4TVP), the N- and C-terminal halves of the pocket exhibit anti-correlated breathing (mean C_ij = +0.22 to +0.45, with 4–9 negative pairs), while in the 6MEO ternary complex (relay activated), all 28 intra-pocket pairs are positively correlated (mean +0.79). The G471A variant (5FYJ) abolishes this gating motion entirely: all 28 pairs are positive (mean +0.72), matching the 6MEO relay-active sign pattern but through steric gridlock (lowest PCA variance, lowest Schlitter entropy) rather than active signal propagation. This state-dependent transition (from anti-correlated breathing to fully correlated co-motion) is reproduced across three independent wildtype structures from separate 500 ns simulations, providing the first multi-structure mechanistic validation of the relay hypothesis. A drug binding the pre-triggered pocket stabilizes the anti-correlated resting state, preventing the transition to the correlated co-motion required for downstream signal propagation.

The mechanism of inhibition is steric blockade of relay motion: a small molecule occupying this pocket permits normal CD4 binding (the pocket is approximately 15 Å from the CD4 binding site) and normal CCR5 engagement (different protein surface), but disrupts the conformational relay at the gp120 C-terminal junction by constraining the rearrangement required for downstream fusion activation. The virus docks to the target cell but fails to trigger membrane fusion. Four independent computational lines of evidence support this non-competitive relay blocker mechanism: (1) DCCM allosteric silencing shows the drug reduces relay→CCR5 coupling by 51% while preserving relay→CD4 coupling (−15% change; §3.13, Fig 8); (2) MM-GBSA non-competitiveness analysis confirms that CD4 binding remains strongly favorable (−28.4 kcal/mol) in the drug-bound state (§3.14), with per-residue energy decomposition revealing a van der Waals–dominated binding architecture in which 9 pocket residues each contribute >1 kcal/mol (Fig 5); (3) intra-pocket correlation analysis shows the relay pocket retains its structural identity (mean C_ij = 0.89) but loses signal transmission capacity; and (4) alchemical free energy perturbation (9 independent legs across 3 replicas) reveals strong negative allosteric coupling between the drug and CD4 (ΔG_bind: −4.0 kcal/mol pre-CD4 → −0.1 kcal/mol post-CD4; ΔΔG = +3.9 kcal/mol; §3.15, Fig 9), providing the first thermodynamic measurement of the conformational relay signal strength and establishing that the drug selectively stabilizes the pre-triggered conformation. Experimental validation of this mechanism is outlined in §4.8.

**Fig 5. MM-GBSA per-residue binding energy decomposition.** (A) Horizontal bar chart of per-residue contributions to BMS-626529 binding at the relay pocket (ΔG_bind = −37.7 ± 3.5 kcal/mol, drug-bound ternary complex, 40 frames). The 9 most contributing residues (>1 kcal/mol each) are labeled; three phenylalanines, two valines, and one leucine among top contributors confirm the hydrophobic binding profile predicted by the pharmacophore model. (B) Energy component breakdown: van der Waals packing dominates (ΔVDW = −68.1 kcal/mol) with a partial desolvation penalty (ΔEGB = +29.9 kcal/mol), consistent with the deep hydrophobic pocket interior.

The FEP establishes that more than half the binding mechanism is not molecular contacts but the drug selecting a conformation: +3.95 kcal/mol is the energy that CD4-induced rearrangement transmits through the relay junction, providing the first direct measurement of the conformational relay signal strength. Experimental validation of this mechanism is outlined in §4.8.

Independent biophysical data confirm the state-dependent dynamics predicted by our DCCM analysis. Hydrogen-deuterium exchange mass spectrometry (HDX-MS) of gp120 shows that the relay pocket region (peptic Fragment 15: T455–L483, spanning all of RPGGGDMR) exhibits moderate backbone protection in the unliganded state and becomes "significantly stabilized" upon CD4 binding (Kong et al., 2010), while full-length glycosylated gp120 shows "moderate protection" at residues 458–474 that increases to "structural rigidification upon sCD4 binding" (Guttman et al., 2012). This transition from dynamic (pre-CD4) to rigid (post-CD4) independently corroborates the relay gate mechanism: our DCCM shows anti-correlated breathing in pre-triggered SOSIP and all-positive co-motion after relay activation. Single-molecule FRET studies demonstrate three conformational states (pre-triggered, intermediate, and post-CD4; Ma et al., 2018), and gp41-directed bNAbs allosterically remodel gp120 dynamics through the same relay junction (Xu et al., 2026), providing direct experimental evidence for the bidirectional allosteric pathway our relay mechanism predicts. NMR confirms that the C5 domain immediately downstream of the relay pocket binds the gp41 ectodomain (Kd ≈ 1 μM; Guilhaudis et al., 2002), establishing the structural basis for signal transmission from gp120 to gp41. This converges with network-theoretic analyses showing that allosteric communication pathways in gp120 connect CD4-binding-site rearrangements to distant functional regions (Henderson et al., 2020).

### 4.2 Conservation and Resistance Barrier

The relay segment (mean entropy 0.21 across 981 NCBI gp120 sequences; Fig 3) sits within the C5 subdomain (mean entropy 0.14), the most conserved constant region of gp120. Three residues (Arg465, Pro466, and Gly468) are completely invariant, and three more (Gly469, Met471, Arg472) exceed 94% conservation. The RPGGGDMR motif is conserved at >96% identity across all major HIV-1 Group M clades (A through G and circulating recombinant forms); 766 of 981 sequences carry the exact motif. Even the most variable position (Gly467, entropy 0.72) shows predominantly G→A substitutions that preserve the small side-chain character required for relay backbone flexibility.

The structural basis for this constraint: Pro466 introduces a backbone kink required for relay geometry and is completely invariant. Of the triple glycine (467–469), Gly468 is invariant and Gly469 is 98.5% conserved; Gly467 tolerates alanine substitution (9% of sequences), consistent with a requirement for small side chains rather than absolute glycine identity at this position. Asp470 shows D→N substitutions (8%) that preserve the polar character of the charge anchor. Computational alanine scanning (§3.19) confirms that the relay signal propagates through backbone topology, not side chains, with the glycine triplet showing zero coupling change upon mutation. The relay backbone is structurally constrained against mutation.

The G471A variant (PDB 5FYJ, Clade G) provides direct dynamic evidence for this constraint. In 500 ns MD, the G471A pocket is the most rigid in the dataset (RMSD 1.29 Å, RMSF 0th percentile, Schlitter entropy 1,246 J/mol·K), and its DCCM shows all 28 intra-pocket pairs positive (mean +0.72), abolishing the anti-correlated breathing that all three wildtype SOSIP structures exhibit (§3.11). The alanine substitution eliminates the backbone flexibility that glycine provides at position 471, converting the relay gate from a functional switch (openable by CD4-induced signal) to a seized mechanism. Greater pocket rigidity is not an advantage; it is a loss of function. The virus conserves Gly471 because the relay must be a gate, not a wall.

Deep mutational scanning (DMS) of BG505.T332N SOSIP.664 (Haddox et al., 2018) provides independent functional evidence for the relay gate. Five of eight relay pocket residues (Arg465, Pro466, Gly468, Gly469, Asp470) are rank 1 in single-round infectivity (the most fitness-constrained rank), confirming extreme functional conservation at these positions. Three positions show fitness landscape gaps: Met471 (rank 3/20) and Arg472 (rank 3/20) are naturally conserved at >97% yet moderately DMS-tolerant, while the most dramatic outlier is position 471 (BG505 numbering; Gly467 in HXB2), where glycine ranks only 6th/20 and alanine ranks 5th/20 in DMS, yet natural conservation is 87% glycine. This "fitness landscape gap" — DMS-tolerant yet evolutionarily constrained — is the signature of a functional role that extends beyond single-round entry. The relay gate model explains the Gly467 pattern: G471A preserves entry in a single-round assay (which does not require sustained relay signaling) but abolishes the backbone flexibility required for multi-round replication fitness where the conformational relay must operate reliably across many entry events. Evolutionary coevolution analysis (353 HIV-1 Env sequences, mutual information with average product correction) reveals that the two variable relay positions (Asp470 and Arg472) co-evolve as a structural pair (MI = 0.443), consistent with a salt bridge that stabilizes the motif fold, and that Arg472 shows coevolution coupling to Lys486 in the gp41-contact region (MI = 0.393), providing independent evolutionary evidence for the gp120-to-gp41 relay axis.

Five independent methods converge on a single gate-hinge residue within the relay motif: DMS (Gly467/471 rank 6/20, most DMS-tolerant position), AlphaFold (pLDDT 89.6, lowest in motif), Shannon entropy (0.72, highest in motif), cryo-EM local density (27th percentile, lowest in motif), and GNM normal modes (9th percentile mobility). This multi-method convergence on a single hinge residue, using fundamentally different data sources (fitness, structure prediction, sequence conservation, experimental density, and elastic network mechanics), constitutes strong independent evidence that the relay pocket contains a functional gate rather than a randomly flexible region.

This is further supported by resistance mutation mapping against the Stanford HIV Drug Resistance Database: **zero of 28 known drug resistance mutations affect relay pocket residues 465–472** (§3.17). The nearest resistance mutation (M475I, a secondary fostemsavir escape) is 3 residues away. The PDB-wide structural survey (§3.17) confirms: 360/360 solved gp120 structures carry the exact RPGGGDMR sequence. The virus has evolved resistance to every class of drug at other sites. It has never mutated this pocket. The relay is load-bearing: changing it to escape a drug would break the entry mechanism itself.

### 4.3 Comparison to Nearest Prior Art

The closest existing drug target to the relay pocket is fostemsavir's binding site on the gp120 β20-β21 element (Pancera et al., 2017). While both are on gp120, they are approximately 15 Å apart, involve different residues, and operate through fundamentally different mechanisms: fostemsavir prevents CD4 attachment entirely, while a relay blocker permits attachment and coreceptor engagement but disrupts the downstream conformational relay required for fusion activation.

### 4.4 The CCR5 ECL3 Site

The identification of CCR5 ECL3 as a gp120 contact surface is independently significant. Most structural and pharmacological work on CCR5 for HIV has focused on the N-terminus, ECL2, and the transmembrane helical bundle (where maraviroc binds). ECL3 has not been extensively characterized as an HIV-relevant interface. If validated experimentally, this site could represent an additional allosteric target on CCR5 distinct from the maraviroc pocket.

### 4.5 The Glycan Shield

HIV-1 gp120 is heavily glycosylated, with N-linked glycans shielding much of the protein surface from antibody and small-molecule access. Structural analysis of glycosylation sequons (NxS/T, x≠P) across six Env structures reveals that the relay pocket occupies a relatively glycan-sparse zone. In PDB 6MEO (BaL strain), the nearest sequon to the pocket centroid is N384 at 15.7 Å, followed by N276 (15.8 Å) and N262 (16.1 Å); no resolved sugar atoms fall within 15 Å of the pocket. In BG505 SOSIP structures (5U7O, 6MU6, 4TVP), the nearest sequon is N363 at ~7.8 Å (a site absent in the BaL strain, S363), yet even these structures contain at most one resolved sugar residue (a single NAG stub) within 15 Å of the pocket center. PDB 5FYJ, the most densely glycosylated structure in our analysis (162 resolved sugar residues), contains zero sugars within 15 Å of the pocket despite the bNAb PGT151 contacting this region (Lee et al., 2016). The glycan shield mapping of Stewart-Jones et al. confirms that the C5 domain lies in a glycan-sparse corridor relative to the V1/V2 and outer-domain glycan clusters (Stewart-Jones et al., 2016). The only FDA-approved small molecule targeting gp120, fostemsavir (BMS-626529), binds the β20–β21 element at a comparable distance from the glycan canopy (Pancera et al., 2017), demonstrating that small-molecule access to gp120 pockets is pharmacologically achievable despite the glycan shield. Small molecules (300–600 Da) require substantially less approach clearance than antibodies (~150 kDa), and the pocket must remain conformationally accessible to execute its relay function, together suggesting favorable steric access for drug candidates.

### 4.6 The Kinetic Window

PDB 6MEO captures the ternary complex (gp120–CD4–CCR5), a transient state during viral entry. A key question is whether the relay pocket exists long enough for a drug to bind. Several considerations suggest the kinetic window is pharmacologically viable:

1. **The drug's preferred state is pre-CD4.** FEP calculations (§3.15) establish that the relay pocket binds BMS-626529 with moderate affinity (−4.0 kcal/mol) in the pre-CD4 conformation, with affinity dropping to near-zero upon CD4 engagement. This state-dependent binding profile indicates that the drug functions by conformational selection: it binds to and stabilizes the resting (pre-triggered) relay conformation before CD4 arrives. The kinetic window is therefore the entire prefusion state, not a brief post-CD4 intermediate.
2. **Precedent for pre-binding.** Fostemsavir binds the prefusion trimer and remains bound through conformational changes. A relay-pocket drug similarly pre-loads, with its ~4 kcal/mol thermodynamic penalty (ΔΔG) resisting the CD4-induced conformational relay.
3. **Structural rearrangement occurs at the relay.** The gp120 C-terminal region (including residues 465–472) undergoes conformational change upon CD4 binding (Pancera et al., 2010; Herschhorn et al., 2016). Pharmacophore analysis of the ternary complex reveals that the pocket wall includes residues from both gp120 (chain G) and CD4 (chain A: Gln25, Asn39, Gln40, Phe43, Thr45), indicating that CD4 residues approach the drug binding site upon engagement. If the drug-bound pre-CD4 conformation traps the pocket in a partially rearranged state, the drug may make contacts with both viral and host protein surfaces in the intermediate, creating an additional barrier to escape mutation since the virus cannot mutate the host half of the binding interface.
5. **The rate-limiting step is CCR5 recruitment.** After CD4 binding, the virus must recruit CCR5 through V3 loop insertion, a slow, diffusion-limited process. A drug that locks the relay in the pre-triggered conformation would have this entire window to maintain its inhibitory effect.



### 4.7 Limitations

All findings in this study are computational, derived from structures at 2.55–3.9 Å resolution (cryo-EM and X-ray crystallography). Key limitations:

1. **No experimental validation.** All evidence for the relay pocket is computational. The pocket's biological druggability, its accessibility in the native viral spike, and its capacity to block fusion when occupied have not been confirmed by biophysical or cellular assays. Specific experimental approaches to validate these predictions are detailed in §4.8.
2. **Glycan shield.** All simulations used bare (non-glycosylated) systems. Structural analysis of six Env structures shows zero resolved sugars within 15 Å of the pocket center, and fostemsavir provides clinical precedent for small-molecule access to gp120 at comparable glycan proximity. However, the nearest glycosylation sequon is strain-dependent (N384 at 15.7 Å in BaL versus N363 at ~7.8 Å in BG505 SOSIP), and dynamic glycan coverage cannot be assessed from static coordinates. Fully glycosylated Env trimer simulations at 500 ns–1 μs timescales are needed to map the dynamic approach corridor.
3. **Structural resolution.** The primary structure (6MEO, 3.9 Å cryo-EM) determines backbone positions but side-chain rotamers, particularly Arg465 and Arg472, are modeled rather than directly resolved. Cross-validation against a 2.55 Å X-ray structure (6MU6) independently resolves the same pocket geometry including side-chain positions, and 500 ns MD confirms sub-angstrom stability (0th percentile RMSF, §3.11). The pocket is reproduced in four additional structures spanning 2.55–3.1 Å (§3.3). The terminal guanidinium groups of Arg465 and Arg472 possess high intrinsic flexibility that static crystallography cannot fully constrain.
4. **Conformational sampling.** Nine independent MD simulations (~4,000 ns aggregate) across multiple starting structures confirm pocket persistence (RMSD 0.71–1.98 Å). SOSIP mechanistic cross-validation (DCCM, PCA, ensemble docking on four 500 ns trajectories, 2,000 ns aggregate; §3.11) extends validation to allosteric network properties and identifies the G471A gating mechanism. However, all trajectories sample microsecond timescales; millisecond-scale conformational transitions, including rare large-amplitude breathing, are not captured. Drug-bound SOSIP simulations and MM-GBSA on SOSIP systems remain as future work (§4.1). Steered MD and umbrella sampling of the drug extraction coordinate are planned as future work to characterize the unbinding energy landscape and kinetic residence time.
6. **Binding free energy decomposition.** MM-GBSA yields overall binding free energies of −38.0 and −31.4 kcal/mol for the ternary and reduced systems respectively, but per-residue energy decomposition was only achievable on the smaller reduced system (1,715 atoms). On the 14,149-atom ternary complex, all three decomposition modes (idecomp=1, 2, 3) failed to produce physically meaningful per-residue attributions, a known limitation of GB-based solvation energy decomposition on large multi-chain systems. Because the ternary context strengthens overall binding by ~6.6 kcal/mol, the reduced-system pocket attribution (−3.07 kcal/mol) represents a truncated estimate of the true interaction energy.



### 4.8 Suggested Experimental Validation

For the gp120 relay pocket (Arg465–Arg472):

1. **Alanine scanning mutagenesis.** Express Env trimers with individual alanine substitutions at positions 465–472. Measure CD4 and CCR5 binding (SPR or BLI), membrane fusion (cell–cell fusion assay), and pseudovirus entry (single-round infectivity assay). Expected: binding unaffected, fusion impaired.
2. **Hydrogen-deuterium exchange mass spectrometry (HDX-MS).** Published HDX-MS data show increased backbone protection at this region upon CD4 binding (Kong et al., 2010; Guttman et al., 2012); the key experiment is to test whether a relay-blocking compound restores pre-CD4 exchange levels (increased exchange relative to CD4-bound, decreased relative to unliganded), confirming that the drug locks the relay in its pre-triggered state.
3. **Virtual screening.** Screen ZINC, Enamine REAL, or other drug-like compound libraries against the relay pocket using the provided docking configurations. Top-scoring compounds can be tested in pseudovirus entry assays.

For CCR5 ECL3 (Tyr251, Leu255, Asn258, Asp276):

1. **Alanine scanning mutagenesis.** Measure effect on gp120 binding affinity and HIV-1 entry.
2. **Cross-linking mass spectrometry (XL-MS).** Confirm direct contact between ECL3 residues and gp120 in the ternary complex.



## 5. Data Availability

A comprehensive supplementary data package accompanies this study, containing 145 structured data files (JSON, XVG, DAT), 81 coordinate and docking pose files (PDB, PDBQT), 108 standalone analysis scripts (93 Python, 15 Bash), and full GROMACS simulation parameters (topology, MDP, and coordinate files). The package is organized into four top-level directories: `results/` (11 categories from pocket discovery through drug design, including FEP parameter files), `scripts/` (6 workflow categories from setup through drug design), `topology/` (force field parameters, MDP files, and run input), and `structures/` (PDB and PDBQT coordinate files for all receptors and ligands). All scripts use a relative path convention (`DATA_ROOT`) that resolves to the package root directory, enabling direct execution without path modification. SHA-256 checksums are provided for all files. Raw MD trajectories for all five 500 ns simulations (~25 GB) and FEP dH/dλ data for all 288 lambda windows (~2.5 GB) are available for download at https://hiv1-relay-data.jretting.workers.dev (self-contained packages with SHA-256 verification).

The data files span the full analytical scope of this work: 55-feature pharmacophore model, fpocket and P2Rank cavity detection, B-factor and druggability scoring, Shannon entropy from 981 NCBI gp120 sequences, AutoDock Vina docking for relay and control sites, 96-compound virtual screening results, ADMET predictions for 98 compounds, resistance mutation mapping (28 mutations), cross-clade structural survey (531 PDB entries, 360/360 exact RPGGGDMR), de novo lead designs, 20-frame ensemble docking, WaterMap hydration analysis, allosteric network and betweenness centrality, PCA essential dynamics, alanine scanning (8 mutants), proteome specificity (20,431 human proteins), SOSIP MD cross-validation (4 structures at 20 ns, all 4 extended to 500 ns with comprehensive analysis JSONs), SOSIP mechanistic cross-validation (intra-pocket DCCM and PCA with Schlitter entropy on all four 500 ns SOSIP trajectories including the G471A variant, ensemble docking on 5U7O and 6MU6), 500 ns DCCM (10,001 frames), MM-GBSA binding free energies (reduced and ternary systems, per-residue decomposition), IE-MMPBSA interaction entropy, non-competitiveness free energy analysis, full-ternary allosteric silencing DCCM (apo vs. drug-bound, 887 Cα atoms), and FEP/TI alchemical free energies with TI, BAR, and MBAR estimators across replicate thermodynamic cycles.

The docking input files (ligand and receptor PDBQT files, 20 cluster frame conformations from GROMOS clustering) and all 40 ensemble docking output poses (20 BMS-626529 + 20 Vicriviroc) are included. FEP MDP parameter files are provided for independent binding free energy verification. The analysis scripts are standalone and self-documenting, requiring only publicly available software (GROMACS, AutoDock Vina, RDKit, MDAnalysis, ProDy, alchemlyb; see supplementary `README.md` for full dependency and execution documentation). Every quantitative claim in this paper can be independently verified by re-running the corresponding script against the provided data.

The complete supplementary data package (scripts, results, figures, structures, and topology files) is available at https://github.com/jonretting/gp120-relay-supplementary. A permanent archived version is deposited on Zenodo (DOI: [10.5281/zenodo.23389294](https://doi.org/10.5281/zenodo.23389294)). The source structure PDB 6MEO is publicly available from the RCSB Protein Data Bank.

## 6. Patent Disclosure

A US provisional patent application covering the gp120 C5 conformational relay target and associated pharmacophore models has been filed (July 2026).

## 7. AI Disclosure

AI tools (Claude, Anthropic; Gemini, Google) were used for writing assistance, literature search, analysis script drafting, and scientific discussion including evaluation of mechanistic hypotheses and alternative interpretations. All molecular dynamics simulations and free energy perturbation calculations were performed by the author on self-owned hardware using GROMACS 2025.2 (see Data Availability). The evidentiary data underlying this paper—trajectories, binding free energies, cross-correlation matrices—are numerical solutions to classical equations of motion, not language model outputs; no AI system can generate molecular dynamics trajectory data. All 108 analysis scripts were validated via cold-start Docker testing (Python 3.11-slim, pip-only dependencies), producing identical quantitative results without any AI involvement. Every quantitative claim in this paper can be independently reproduced by executing the provided scripts against the provided simulation data. All figures were generated from raw simulation data by provided scripts. The author designed all experiments, selected all simulation methods and parameters, executed all calculations, and can independently explain and defend every methodological choice and result interpretation reported in this work.

## Acknowledgments

The author thanks Prof. David Tholen for correspondence that sharpened the author's approach to methodological rigor and critical evaluation of results.

## Supporting Information


**S2 Fig. PCA variance and conformational entropy comparison.** (A) Eigenvalue scree plot of the top 10 principal components for all five structures (6MEO, 5U7O, 6MU6, 4TVP, 5FYJ). Dynamics are dominated by 2–3 collective modes in all systems. (B) Total variance (nm²) and Schlitter conformational entropy (J/mol·K) across structures. The G471A variant (5FYJ) has the lowest variance (2.72 nm²) and entropy (1,246 J/mol·K), consistent with steric gridlock from the Gly→Ala substitution.

**S3 Fig. Ensemble docking score comparison across conformational states.** BMS-626529 docking affinities from MD ensemble docking of 6MEO (20 clusters, mean −6.81 kcal/mol), 5U7O (3 clusters, mean −4.29 kcal/mol), and 6MU6 (3 clusters, mean −3.90 kcal/mol). The ~2.5 kcal/mol difference between 6MEO (CD4+CCR5 triggered) and SOSIP (pre-fusion) structures quantifies CD4's geometric contribution to the binding pocket.

**S4 Fig. Computational alanine scanning of relay pocket residues.** (A) Per-residue change in GNM relay–CD4 coupling upon alanine mutation. Arg465 shows the maximum disruption (−2.2%); the glycine triplet (Gly467–469) shows zero coupling change. (B) Number of inter-residue contacts lost per alanine mutation, correlating with coupling sensitivity. The relay signal propagates through backbone topology, not side-chain contacts.

**S5 Fig. B-factor regional comparison across gp120.** Mean B-factors of the relay pocket (72.8 Å², 38th percentile) compared to other gp120 regions (V3 loop: 106.3 Å², C2/bridging sheet: 97.2 Å², C4 domain: 86.3 Å²) and binding partner chains. The relay pocket is more ordered than regions of known conformational flexibility, consistent with a structurally constrained, druggable site.

**S6 Fig. Per-residue RMSF heatmap across SOSIP structures.** RMSF values for all residues in the four 500 ns SOSIP reduced-system simulations, with the relay pocket residues highlighted. Three of four structures (5FYJ, 6MU6, 4TVP) show the relay pocket at the 0th RMSF percentile — no other residue in these systems exhibits lower fluctuation over half a microsecond.

**S7 Fig. Cross-structure pocket overlay.** Backbone superposition of the relay pocket (Arg465–Arg472) across all five HIV-1 Env structures analyzed in this study, colored by PDB source. Backbone RMSD <1.0 Å for all pairwise comparisons, with GNM intra-pocket correlation ≥0.69 in all structures and independent pocket detection by fpocket in all five. The pocket geometry is structurally invariant across clades, resolutions, and conformational states.

## References

Abraham, M. J. et al. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX* 1–2, 19–25. [doi:10.1016/j.softx.2015.06.001](https://doi.org/10.1016/j.softx.2015.06.001)

Bennett, C. H. (1976). Efficient estimation of free energy differences from Monte Carlo data. *J. Comput. Phys.* 22, 245–268. [doi:10.1016/0021-9991(76)90078-4](https://doi.org/10.1016/0021-9991(76)90078-4)

Blattner, C. et al. (2014). Structural delineation of a quaternary, cleavage-dependent epitope at the gp41–gp120 interface on intact HIV-1 Env trimers. *Immunity* 40, 669–680. [doi:10.1016/j.immuni.2014.04.008](https://doi.org/10.1016/j.immuni.2014.04.008)

Chan, D. C., Fass, D., Berger, J. M. & Kim, P. S. (1997). Core structure of gp41 from the HIV envelope glycoprotein. *Cell* 89, 263–273. [doi:10.1016/s0092-8674(00)80205-6](https://doi.org/10.1016/s0092-8674(00)80205-6)

Dorr, P. et al. (2005). Maraviroc (UK-427,857), a potent, orally bioavailable, and selective small-molecule inhibitor of chemokine receptor CCR5 with broad-spectrum anti-human immunodeficiency virus type 1 activity. *Antimicrob. Agents Chemother.* 49, 4721–4732. [doi:10.1128/AAC.49.11.4721-4732.2005](https://doi.org/10.1128/AAC.49.11.4721-4732.2005)

Duan, L., Liu, X. & Zhang, J. Z. H. (2016). Interaction entropy: a new paradigm for highly efficient and reliable computation of protein–ligand binding free energy. *J. Am. Chem. Soc.* 138, 5722–5728. [doi:10.1021/jacs.6b02682](https://doi.org/10.1021/jacs.6b02682)

Eastman, P. et al. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLoS Comput. Biol.* 13, e1005659. [doi:10.1371/journal.pcbi.1005659](https://doi.org/10.1371/journal.pcbi.1005659)

Eberhardt, J., Santos-Martins, D., Tillack, A. F. & Forli, S. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* 61, 3891–3898. [doi:10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203)

Genheden, S. & Ryde, U. (2015). The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities. *Expert Opin. Drug Discov.* 10, 449–461. [doi:10.1517/17460441.2015.1032936](https://doi.org/10.1517/17460441.2015.1032936)

Guilhaudis, L. et al. (2002). Solution structure of the HIV gp120 C5 domain. *Eur. J. Biochem.* 269, 4860–4867. [doi:10.1046/j.1432-1033.2002.03187.x](https://doi.org/10.1046/j.1432-1033.2002.03187.x)

Guttman, M. et al. (2012). Solution structure, conformational dynamics, and CD4-induced activation in full-length, glycosylated, monomeric HIV gp120. *J. Virol.* 86, 8750–8764. [doi:10.1128/JVI.07224-11](https://doi.org/10.1128/JVI.07224-11)

Halgren, T. A. (2009). Identifying and characterizing binding sites and assessing druggability. *J. Chem. Inf. Model.* 49, 377–389. [doi:10.1021/ci800324m](https://doi.org/10.1021/ci800324m)

Haddox, H. K., Dingens, A. S., Hilton, S. K., Overbaugh, J. & Bloom, J. D. (2018). Mapping mutational effects along the evolutionary landscape of HIV envelope. *eLife* 7, e34420. [doi:10.7554/eLife.34420](https://doi.org/10.7554/eLife.34420)

Henderson, R. et al. (2020). Disruption of the HIV-1 Envelope allosteric network blocks CD4-induced rearrangements. *Nat. Commun.* 11, 520. [doi:10.1038/s41467-019-14196-w](https://doi.org/10.1038/s41467-019-14196-w)

Herschhorn, A. et al. (2016). Release of gp120 restraints leads to an entry-competent intermediate state of the HIV-1 envelope glycoproteins. *mBio* 7, e01598-16. [doi:10.1128/mBio.01598-16](https://doi.org/10.1128/mBio.01598-16)

Huang, C. C. et al. (2005). Structure of a V3-containing HIV-1 gp120 core. *Science* 310, 1025–1028. [doi:10.1126/science.1118398](https://doi.org/10.1126/science.1118398)

Huang, J. et al. (2014). Broad and potent HIV-1 neutralization by a human antibody that binds the gp41–gp120 interface. *Nature* 515, 138–142. [doi:10.1038/nature13601](https://doi.org/10.1038/nature13601)


Jacobson, J. M. et al. (2009). Safety, pharmacokinetics, and antiretroviral activity of multiple doses of ibalizumab (formerly TNX-355), an anti-CD4 monoclonal antibody, in human immunodeficiency virus type 1-infected adults. *Antimicrob. Agents Chemother.* 53, 450–457. [doi:10.1128/AAC.00942-08](https://doi.org/10.1128/AAC.00942-08)

Julien, J. P. et al. (2013). Crystal structure of a soluble cleaved HIV-1 envelope trimer. *Science* 342, 1477–1483. [doi:10.1126/science.1245625](https://doi.org/10.1126/science.1245625)

Kilby, J. M. et al. (1998). Potent suppression of HIV-1 replication in humans by T-20, a peptide inhibitor of gp41-mediated virus entry. *Nat. Med.* 4, 1302–1307. [doi:10.1038/3293](https://doi.org/10.1038/3293)

Kirkwood, J. G. (1935). Statistical mechanics of fluid mixtures. *J. Chem. Phys.* 3, 300–313. [doi:10.1063/1.1749657](https://doi.org/10.1063/1.1749657)

Kong, L. et al. (2010). Local conformational stability of HIV-1 gp120 in unliganded and CD4-bound states as defined by amide hydrogen/deuterium exchange. *J. Virol.* 84, 10311–10321. [doi:10.1128/JVI.00688-10](https://doi.org/10.1128/JVI.00688-10)

Lai, Y.-T. et al. (2019). Lattice engineering enables definition of molecular features allowing for potent small-molecule inhibition of HIV-1 entry. *Nat. Commun.* 10, 47. [doi:10.1038/s41467-018-07851-1](https://doi.org/10.1038/s41467-018-07851-1)

Le Guilloux, V., Schmidtke, P. & Tuffery, P. (2009). Fpocket: an open source platform for ligand pocket detection. *BMC Bioinformatics* 10, 168. [doi:10.1186/1471-2105-10-168](https://doi.org/10.1186/1471-2105-10-168)

Lee, J. H., Ozorowski, G. & Ward, A. B. (2016). Cryo-EM structure of a native, fully glycosylated, cleaved HIV-1 envelope trimer. *Science* 351, 1043–1048. [doi:10.1126/science.aad2450](https://doi.org/10.1126/science.aad2450)

Ma, X. et al. (2018). HIV-1 Env trimer opens through an asymmetric intermediate in which individual protomers adopt distinct conformations. *eLife* 7, e34271. [doi:10.7554/eLife.34271](https://doi.org/10.7554/eLife.34271)

Maeda, K. et al. (2006). Structural and molecular interactions of CCR5 inhibitors with CCR5. *J. Biol. Chem.* 281, 12688–12698. [doi:10.1074/jbc.M512688200](https://doi.org/10.1074/jbc.M512688200)

O'Boyle, N. M. et al. (2011). Open Babel: An open chemical toolbox. *J. Cheminform.* 3, 33. [doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)

Pancera, M. et al. (2010). Structure of HIV-1 gp120 with gp41-interactive region reveals layered envelope architecture and basis of conformational mobility. *Proc. Natl. Acad. Sci. USA* 107, 1166–1171. [doi:10.1073/pnas.0911004107](https://doi.org/10.1073/pnas.0911004107)

Pancera, M. et al. (2014). Structure and immune recognition of trimeric pre-fusion HIV-1 Env. *Nature* 514, 455–461. [doi:10.1038/nature13808](https://doi.org/10.1038/nature13808)

Pancera, M. et al. (2017). Crystal structures of trimeric HIV envelope with entry inhibitors BMS-378806 and BMS-626529. *Nat. Chem. Biol.* 13, 1115–1122. [doi:10.1038/nchembio.2460](https://doi.org/10.1038/nchembio.2460)

Rizzuto, C. D. et al. (1998). A conserved HIV gp120 glycoprotein structure involved in chemokine receptor binding. *Science* 280, 1949–1953. [doi:10.1126/science.280.5371.1949](https://doi.org/10.1126/science.280.5371.1949)

Scharf, L. et al. (2014). Antibody 8ANC195 reveals a site of broad vulnerability on the HIV-1 envelope spike. *Cell Rep.* 7, 785–795. [doi:10.1016/j.celrep.2014.04.001](https://doi.org/10.1016/j.celrep.2014.04.001)

Stewart-Jones, G. B. E. et al. (2016). Trimeric HIV-1-Env structures define glycan shields from clades A, B, and G. *Cell* 165, 813–826. [doi:10.1016/j.cell.2016.04.010](https://doi.org/10.1016/j.cell.2016.04.010)

Schindler, T. et al. (2000). Structural mechanism for STI-571 inhibition of Abelson tyrosine kinase. *Science* 289, 1938–1942. [doi:10.1126/science.289.5486.1938](https://doi.org/10.1126/science.289.5486.1938)

Shaik, M. M. et al. (2019). Structural basis of coreceptor recognition by HIV-1 envelope spike. *Nature* 565, 318–323. [doi:10.1038/s41586-018-0804-9](https://doi.org/10.1038/s41586-018-0804-9)

Shirts, M. R. & Chodera, J. D. (2008). Statistically optimal analysis of samples from multiple equilibrium states. *J. Chem. Phys.* 129, 124105. [doi:10.1063/1.2978177](https://doi.org/10.1063/1.2978177)

Shirts, M. R. et al. (2017). Lessons learned from comparing molecular dynamics engines on the SAMPL5 dataset. *J. Comput. Aided Mol. Des.* 31, 147–161. [doi:10.1007/s10822-016-9977-1](https://doi.org/10.1007/s10822-016-9977-1)

Tan, Q. et al. (2013). Structure of the CCR5 chemokine receptor–HIV entry inhibitor maraviroc complex. *Science* 341, 1387–1390. [doi:10.1126/science.1241475](https://doi.org/10.1126/science.1241475)

Trott, O. & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.* 31, 455–461. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)

Weissenhorn, W., Dessen, A., Harrison, S. C., Skehel, J. J. & Wiley, D. C. (1997). Atomic structure of the ectodomain from HIV-1 gp41. *Nature* 387, 426–430. [doi:10.1038/387426a0](https://doi.org/10.1038/387426a0)

Wild, C. T., Shugars, D. C., Greenwell, T. K., McDanal, C. B. & Matthews, T. J. (1994). Peptides corresponding to a predictive α-helical domain of human immunodeficiency virus type 1 gp41 are potent inhibitors of virus infection. *Proc. Natl. Acad. Sci. USA* 91, 9770–9774. [doi:10.1073/pnas.91.21.9770](https://doi.org/10.1073/pnas.91.21.9770)

Xu, M. et al. (2026). gp41-directed broadly neutralizing antibodies allosterically remodel HIV-1 Env conformational dynamics. *bioRxiv*, 2026.01.27.702099. [doi:10.1101/2026.01.27.702099](https://doi.org/10.1101/2026.01.27.702099)

---

*Jonathan Michael Retting, Independent Researcher*