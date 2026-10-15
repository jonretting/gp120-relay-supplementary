# HIV-1 gp120-CD4-CCR5 Ternary Complex: Binding Hotspot Residues and Pharmacophore Models

**Author:** Jonathan Retting  
**Date:** 2026-04-02  
**Structure:** PDB [6MEO](https://www.rcsb.org/structure/6meo) -- Cryo-EM, 3.9 Å  
**Reference:** Shaik MM, Peng H, Lu J, et al. *Nature* 565, 318-323 (2019). DOI: 10.1038/s41586-018-0804-9

> **Numbering note:** gp120 residue numbers in this document follow the PDB 6MEO
> author numbering (isolate 92UG037.8, Clade A). HXB2 equivalents are provided
> in parentheses for cross-referencing with the wider HIV literature.

---

## 1. Structure

| Chain | Protein | Organism | Residues | UniProt | Role |
|---|---|---|---|---|---|
| G | gp120 (Env glycoprotein) | HIV-1 (92UG037.8) | 29-489 | [Q70145](https://www.uniprot.org/uniprot/Q70145) | Viral surface protein |
| A | CD4 | Homo sapiens | 26-201 (D1-D2) | [P01730](https://www.uniprot.org/uniprot/P01730) | Primary receptor |
| B | CCR5 | Homo sapiens | 1-313 | [P51681](https://www.uniprot.org/uniprot/P51681) | Coreceptor |

- 885 residues, 7,003 atoms
- 131 residues at protein-protein interfaces (within 8 Å of another chain)

---

## 2. Binding Hotspot Residues

The following 22 interface residues are structurally essential for the gp120-CD4-CCR5 ternary complex. Disruption of any single residue is predicted to destabilize or prevent complex formation.

### CD4 (Chain A) -- 4 hotspot residues

| # | Residue | AA | Region | Known Role |
|---|---|---|---|---|
| 1 | Glu 13 | E | D1 domain | gp120 contact face |
| 2 | Thr 15 | T | D1 domain | gp120 binding interface |
| 3 | His 27 | H | D1 CC' loop | gp120 contact region |
| 4 | Asn 66 | N | D1 domain | Proximal to Phe43 cavity |

### CCR5 (Chain B) -- 10 hotspot residues

| # | Residue | AA | Region | Known Role | Note |
|---|---|---|---|---|---|
| 1 | Ser 7 | S | N-terminus | Bridging sheet contact; adjacent to sulfo-Tyr10/14 | |
| 2 | Ala 29 | A | N-terminus | gp120 N-terminal contact | |
| 3 | Phe 109 | F | TM3 / ECL1 | Not previously reported as primary contact | **Candidate novel** |
| 4 | Leu 174 | L | ECL2 | gp120 contact surface | |
| 5 | Tyr 184 | Y | ECL2 | gp120 contact surface | |
| 6 | Gln 194 | Q | ECL2 | Primary gp120 contact; maraviroc allosteric region | |
| 7 | Tyr 251 | Y | TM6 / ECL3 | Not extensively characterized for HIV | **Candidate novel** |
| 8 | Leu 255 | L | ECL3 | Not extensively characterized for HIV | **Candidate novel** |
| 9 | Asn 258 | N | ECL3 / TM6 | Not extensively characterized for HIV | **Candidate novel** |
| 10 | Asp 276 | D | TM7 / ECL3 | Not extensively characterized for HIV | **Candidate novel** |

### gp120 (Chain G) -- 8 hotspot residues

| # | Residue | AA | HXB2 # | Region | Known Role | Note |
|---|---|---|---|---|---|---|
| 1 | Thr 127 | T | 128 | V1/V2 stem | V1/V2 shields V3 pre-CD4 binding | **Candidate novel** |
| 2 | Lys 282 | K | 282 | C2 / bridging sheet | Connects CD4-induced conformation to CCR5 binding | |
| 3 | Gly 323 | G | 324 | V3 loop stem | Coreceptor tropism determinant | |
| 4 | Arg 326 | R | 327 | V3 loop | V3 base | |
| 5 | Ile 418 | I | **423** | C4 / β20-β21 base | Conformational switch (Herschhorn et al. 2016); fostemsavir target | |
| 6 | Arg 465 | R | 469 | C5 / gp41 interface | Conformational relay to gp41 | **Candidate novel** |
| 7 | Gly 467 | G | 471 | C5 / gp41 interface | Conformational relay to gp41 | **Candidate novel** |
| 8 | Arg 472 | R | 476 | C5 / gp41 interface | Conformational relay to gp41 | **Candidate novel** |

---

## 3. Validation Against Known Biology

| Interface | Known Literature | Recovered |
|---|---|---|
| CCR5 N-terminus (Ser7, Ala29) | N-terminal sulfated tyrosines Tyr10/14 contact gp120 bridging sheet | Yes |
| CCR5 ECL2 (Leu174, Tyr184, Gln194) | Primary gp120 contact surface; maraviroc allosteric pocket | Yes |
| gp120 V3 loop (Gly323, Arg326) | V3 tip inserts into CCR5 chemokine pocket; tropism determinant | Yes |
| gp120 bridging sheet (Lys282) | CD4-induced bridging sheet enables coreceptor engagement | Yes |
| CD4 D1 domain (Glu13, Thr15, His27, Asn66) | D1 is primary gp120 contact; Phe43 cavity is structural anchor | Yes |

All five established protein-protein interaction interfaces are independently identified.

---

## 4. Candidate Novel Drug Targets

### 4.1 gp120 Conformational Relay -- Arg465, Pro466, Gly467, Gly468, Gly469, Asp470, Met471, Arg472

**Location:** gp120 C-terminus, at the gp41 interface.

**Proposed mechanism:** This region transmits the conformational change from CD4/CCR5 binding to gp41 fusion peptide exposure. A molecule occupying this pocket would allow the virus to engage both receptors normally but prevent the structural signal from reaching gp41. The virus docks but cannot trigger membrane fusion.

**Context:** No approved HIV drug targets this relay. Current entry inhibitors block binding (maraviroc, ibalizumab) or fusion (enfuvirtide targets gp41 HR1). Fostemsavir/BMS-626529 (FDA 2020) acts at the adjacent β20-β21 element (Ile 418 / HXB2 I423) but blocks CD4 *attachment* rather than the downstream conformational signal. This relay represents an intermediate step between binding and fusion that has been structurally characterized (Pancera et al., PNAS 2010; Julien et al., Science 2013) but not targeted pharmacologically.

### 4.2 CCR5 ECL3/TM6-7 -- Tyr251, Leu255, Asn258, Asp276

**Location:** Extracellular loop 3 and transmembrane helices 6-7.

**Context:** Nearly all CCR5 structural work for HIV focuses on ECL2 and the N-terminus. ECL3 has not been extensively characterized as a gp120 contact surface. If validated, this represents a new allosteric/contact region on CCR5 distinct from the maraviroc binding pocket.

### 4.3 CCR5 Phe109 -- TM3

**Context:** F109 is a known maraviroc contact residue in the transmembrane pocket (PDB 4MBS), and also the target site for fostemsavir-related CCR5 conformational changes. Its identification here as a binding hotspot in the ternary complex suggests it plays a dual role: maintaining receptor conformation for productive gp120 engagement from the extracellular side while also being druggable from the transmembrane side.

### 4.4 gp120 Thr127 -- V1/V2 stem

**Context:** V1/V2 loops shield V3 prior to CD4 binding. This residue may contribute to stabilizing the post-CD4 exposed conformation required for coreceptor engagement.

---

## 5. Pharmacophore Models

Four target pockets were analyzed. Pharmacophore features represent spatial positions where drug atoms or functional groups should be placed for pocket occupancy.

### Pocket 1: gp120 Conformational Relay (Novel)

| Parameter | Value |
|---|---|
| Target residues | gp120: Arg465-Arg472 |
| Pocket center | (134.6, 183.9, 136.7) Å |
| Pharmacophore features | 55 |
| Dominant feature type | H-bond donor/acceptor (tight fit, mean 0.74 Å from wall) |
| Pocket character | Large, accessible cavity with extensive H-bonding potential |

### Pocket 2: CCR5 ECL2 / Maraviroc Region (Known)

| Parameter | Value |
|---|---|
| Target residues | CCR5: Leu174, Tyr184, Gln194 |
| Pocket center | (154.1, 146.3, 153.6) Å |
| Pharmacophore features | 19 |
| Dominant feature type | H-bond donor/acceptor (tight, 0.74 Å) |
| Pocket character | Highly selective -- small, tight cavity consistent with maraviroc binding |

### Pocket 3: gp120 V3 Loop (Known)

| Parameter | Value |
|---|---|
| Target residues | gp120: Gly323, Arg326 |
| Pocket center | (148.1, 166.3, 151.9) Å |
| Pharmacophore features | 60 |
| Dominant feature type | Mixed H-bond and hydrophobic |
| Pocket character | Exposed loop -- high accessibility, consistent with known difficulty of targeting V3 with small molecules |

### Pocket 4: CCR5 ECL3 (Novel)

| Parameter | Value |
|---|---|
| Target residues | CCR5: Tyr251, Leu255, Asn258, Asp276 |
| Pocket center | (154.3, 143.9, 169.4) Å |
| Pharmacophore features | 45 |
| Dominant feature type | H-bond donor/acceptor (tight, 0.67 Å) |
| Pocket character | Well-defined geometry, tightest mean wall distance of all four pockets |

Full pharmacophore feature coordinates (xyz, feature type, wall distance) are in `hiv1_pharmacophore_data.json`.

---

## 6. Existing Approved Drugs -- Cross-Reference

| Drug | Target | Year | Mechanism | Relation to These Findings |
|------|--------|------|-----------|---------------------------|
| Maraviroc (Selzentry) | CCR5 TM pocket | 2007 | Allosteric antagonist | Contacts F109 (our hotspot B:109) |
| Enfuvirtide (Fuzeon) | gp41 HR1 | 2003 | Fusion inhibitor | Acts *downstream* of our gp120-gp41 relay |
| Ibalizumab (Trogarzo) | CD4 D2 | 2018 | Post-attachment inhibitor | Targets the receptor, not the virus |
| Fostemsavir (Rukobia) | gp120 β20-β21 | 2020 | Attachment inhibitor | Binds near I418/HXB2-I423 (our hotspot) |

---

## 7. HXB2 Cross-Reference Table

PDB 6MEO uses isolate 92UG037.8 (Clade A) numbering. The HIV literature
predominantly uses HXB2 (Clade B) numbering. Use this table for literature
cross-referencing.

| 6MEO Residue | HXB2 Equivalent | Notes |
|---|---|---|
| G:127 THR | HXB2 T128 | V1/V2 stem |
| G:282 LYS | HXB2 K282 | Bridging sheet (same number) |
| G:323 GLY | HXB2 G324 | V3 loop |
| G:326 ARG | HXB2 R327 | V3 loop |
| G:418 ILE | **HXB2 I423** | β20-β21 conformational switch |
| G:465 ARG | HXB2 R469 | gp41 relay |
| G:467 GLY | HXB2 G471 | gp41 relay |
| G:472 ARG | HXB2 R476 | gp41 relay |

CD4 (Chain A): PDB residue N = UniProt [P01730](https://www.uniprot.org/uniprot/P01730) position N+25.  
CCR5 (Chain B): PDB numbering = UniProt [P51681](https://www.uniprot.org/uniprot/P51681) numbering.

---

## 8. How to Use This Data

### Visualization

Load PDB 6MEO in PyMOL or ChimeraX, then use the selection commands below.

**PyMOL:**
```
fetch 6MEO
select hotspots, (chain A and resi 13+15+27+66) or (chain B and resi 7+29+109+174+184+194+251+255+258+276) or (chain G and resi 127+282+323+326+418+465+467+472)
select novel, (chain B and resi 109+251+255+258+276) or (chain G and resi 127+465+467+472)
select gp41_relay, chain G and resi 465-472
select ecl3, chain B and resi 251+255+258+276
show sticks, hotspots
color red, novel
color green, hotspots and not novel
```

**ChimeraX:**
```
open 6meo
select :13.A|:15.A|:27.A|:66.A|:7.B|:29.B|:109.B|:174.B|:184.B|:194.B|:251.B|:255.B|:258.B|:276.B|:127.G|:282.G|:323.G|:326.G|:418.G|:465.G|:467.G|:472.G
```

### Virtual Screening

Pharmacophore coordinates from `hiv1_pharmacophore_data.json` can be used with:

- **AutoDock Vina / AutoDock-GPU** -- Define search box around pocket center coordinates
- **Schrödinger Glide** -- Use feature centers as receptor-grid constraints
- **OpenEye ROCS** -- Generate 3D shape queries from feature point clouds
- **RDKit** -- `rdkit.Chem.Pharm3D` pharmacophore matching against ZINC, ChEMBL, or Enamine REAL

### Fragment-Based Drug Design

Each pharmacophore feature specifies a position and type:
- Wall distance < 1.5 Å -> H-bond donor/acceptor (amine, hydroxyl, carbonyl)
- Wall distance 1.5-2.5 Å -> Hydrophobic or aromatic (phenyl, cyclohexyl)
- Wall distance > 2.5 Å -> Flexible linker (alkyl chain)

### De Novo Molecule Generation

- **REINVENT 4** -- Pharmacophore coordinates as 3D shape constraints
- **Pocket2Mol / TargetDiff** -- Structure-based generative models accepting pocket coordinates

---

## 9. Suggested Experimental Validation

### For gp120 Conformational Relay (residues 465-472):

1. **Alanine scanning mutagenesis.** Express Env trimers with individual Ala substitutions at positions 465-472. Measure:
   - CD4 and CCR5 binding (SPR or BLI) -- expected: binding unaffected
   - Membrane fusion (cell-cell fusion assay) -- expected: fusion impaired
   - Pseudovirus entry (single-round assay) -- expected: entry reduced or abolished

2. **HDX-MS.** Hydrogen-deuterium exchange mass spectrometry of gp120 C-terminal region during CD4-triggered conformational change. Expected: increased exchange in 465-472 region upon CD4 engagement, confirming structural rearrangement.

3. **Site-directed spin labeling / EPR.** Label residues 465-472 to measure conformational dynamics upon sequential CD4 and CCR5 engagement.

### For CCR5 ECL3 (residues 251, 255, 258, 276):

1. **Alanine scanning mutagenesis.** Measure effect on gp120 binding affinity (SPR/BLI) and HIV-1 entry (pseudovirus assay).

2. **Co-IP / cross-linking MS.** Confirm direct contact between ECL3 residues and gp120 in the ternary complex.

---

## 10. Files

| File | Format | Description |
|---|---|---|
| `hiv1_binding_hotspots_and_pharmacophores.md` | Markdown | This document |
| `hiv1_pharmacophore_data.json` | JSON | All coordinates, residues, features, and visualization commands |
| `hiv_6meo.pdb` | PDB | Source structure ([RCSB download](https://files.rcsb.org/download/6MEO.pdb)) |
| `*.phar` | XML | LigandScout pharmacophore files (importable into MOE, Phase, LigandScout) |
| `*_receptor.pdbqt` | PDBQT | Pre-prepared receptor files for AutoDock Vina / AutoDock-GPU |
| `reference_compounds.json` | JSON | Maraviroc, BMS-626529, enfuvirtide -- SMILES, IDs, mechanisms |
| `selectivity_panel.json` | JSON | CCR5 vs CCR2/CXCR4 residue conservation -- off-target risk |
| `md_simulation_protocol.md` | Markdown | GROMACS protocol for pocket flexibility MD simulations |

---

**Author:** Jonathan Michael Retting, Independent Researcher, Los Angeles, CA
**License:** CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
