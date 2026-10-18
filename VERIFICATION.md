# Independent Verification Guide

## How to Verify These Findings Using Standard Tools

This document describes how to independently verify the binding hotspot
residues and pharmacophore models in this package using publicly available
structural data and standard molecular biology tools.

---

## 1. Cross-Structure Consistency

Download PDB 6MET (alternate conformation of the same complex, same study):

```
curl -O https://files.rcsb.org/download/6MET.pdb
```

Load both 6MEO and 6MET in PyMOL or ChimeraX and verify that the
hotspot residues listed in `results/01-pocket-discovery/hiv1_binding_hotspots_and_pharmacophores.md`
occupy equivalent structural positions in both conformations:

```pymol
fetch 6MEO
fetch 6MET
align 6MET, 6MEO
show sticks, 6MEO and chain G and resi 323+326+418+465+467+472
show sticks, 6MET and chain G and resi 323+326+418+465+467+472
```

All listed hotspot residues should be at the protein-protein interface in
both structures with conserved side-chain orientation.

---

## 2. Published Mutagenesis Cross-Reference

The following published mutagenesis studies provide direct experimental
evidence for the functional importance of the listed residues. These papers
are freely accessible and contain quantitative binding/entry data:

### gp120 Residues

| Residue | What to look up | Reference |
|---------|----------------|-----------|
| K282 | Bridging sheet; conserved structure critical for CCR5 binding | Rizzuto et al. 1998, Science 280:1949 |
| G323, R326 | V3 loop mutations abolish CCR5 binding and entry | Huang et al. 2005, Science 310:1025 |
| I418 | β20-β21 base; I423A (HXB2 numbering) destabilizes State 1 | Herschhorn et al. 2016, mBio 7:e01598-16 |
| R465, G467, R472 | gp120 C-terminal interface with gp41 | Pancera et al. 2010, PNAS 107:1166 |

### CD4 Residues

| Residue | What to look up | Reference |
|---------|----------------|-----------|
| E13, T15 | Domain 1 mutations reduce gp120 binding | Moebius et al. 1992, J Exp Med 176:507 |
| H27, N66 | CD4 D1 structural integrity | Kwong et al. 1998, Nature 393:648 |

### CCR5 Residues

| Residue | What to look up | Reference |
|---------|----------------|-----------|
| S7, A29 | N-terminal truncation/substitution eliminates entry | Dragic et al. 1998, J Virol 72:279 |
| F109 | Direct maraviroc contact in PDB 4MBS crystal structure | Tan et al. 2013, Science 341:1387 |
| L174, Y184 | ECL2 primary gp120 contact surface; mutagenesis affects HIV-1 entry | Dragic et al. 2000, PNAS 97:5639 |
| Q194 | ECL2/TM5 junction; CCR5 mutagenesis identifies role in signaling | Zhou et al. 2000, Eur J Immunol 30:164 |
| Y251 | ECL3; affects HIV-1 and chemokine binding | Maeda et al. 2006, J Biol Chem 281:12688 |
| D276 | TM7; site-directed mutagenesis identifies role in signaling | Zhou et al. 2000, Eur J Immunol 30:164 |

---

## 3. Maraviroc Binding Site Comparison

Download the maraviroc-bound CCR5 crystal structure:

```
curl -O https://files.rcsb.org/download/4MBS.pdb
```

Align the CCR5 chain from 4MBS onto chain B of 6MEO and inspect
the spatial relationship between the maraviroc binding pocket and
the CCR5 hotspot residues:

```pymol
fetch 6MEO
fetch 4MBS
align 4MBS and chain A, 6MEO and chain B
show sticks, 4MBS and resn MRV
show sticks, 6MEO and chain B and resi 109+174+184+194+251+255+258+276
```

Key observations to verify:
- **F109** (TM3) directly contacts maraviroc in the transmembrane pocket
- The ECL2 residues (L174, Y184, Q194) sit on the **extracellular surface**
  above the maraviroc pocket -- a different but structurally connected binding site
- Maraviroc works from inside the transmembrane domain; the listed ECL2 hotspots
  describe the external gp120-contact surface

---

## 4. Pharmacophore Feature Verification

Load the pharmacophore features from `results/01-pocket-discovery/hiv1_pharmacophore_data.json` and
verify that each feature center falls within the binding pocket void and
at appropriate distance from the pocket wall:

```python
import json

with open("results/01-pocket-discovery/hiv1_pharmacophore_data.json") as f:
    data = json.load(f)

for target in data["pharmacophore_targets"]:
    print(f"\n{target['name']}:")
    for feat in target["pharmacophore_features"]:
        xyz = feat["center_xyz_angstrom"]
        wall = feat["wall_distance_angstrom"]
        print(f"  Feature {feat['feature_id']}: "
              f"({xyz[0]:.1f}, {xyz[1]:.1f}, {xyz[2]:.1f}) "
              f"wall_dist={wall:.1f}Å type={feat['feature_type']}")
```

Features should be:
- Inside the pocket void (not clashing with protein atoms)
- 2-5 Å from the nearest pocket wall (appropriate for drug-like interactions)
- Classified as hydrophobic, H-bond donor/acceptor, or aromatic based on
  surrounding residue chemistry

---

## 5. Alanine Scanning Replication

For the most rigorous verification, the hotspot residues can be tested
experimentally via alanine scanning mutagenesis. Priority targets for
experimental validation:

1. **R465, G467, R472** (gp120 C-terminal relay) -- these are the least-studied
   residues and would provide the strongest new evidence
2. **L255, N258** (CCR5 ECL3) -- structural data supports their role but
   direct mutagenesis data is limited

Standard assays: pseudovirus entry assay with Env mutants, or cell-cell
fusion assay measuring syncytia formation with mutant CCR5.

---

## 6. B-factor Verification

Verify the reported B-factor values directly from the PDB file:

```python
with open("structures/hiv_6meo.pdb") as f:
    for line in f:
        if not line.startswith("ATOM"): continue
        if line[21] != "G": continue
        resnum = int(line[22:26].strip())
        if 465 <= resnum <= 472 and line[12:16].strip() == "CA":
            bfactor = float(line[60:66].strip())
            resname = line[17:20].strip()
            print(f"{resname}{resnum} CA B-factor: {bfactor:.2f}")
```

Expected output (from PDB 6MEO):
- ARG465: 86.97, PRO466: 87.53, GLY467: 82.85, GLY468: 77.85
- GLY469: 72.76, ASP470: 68.52, MET471: 69.09, ARG472: 63.63

Chain G mean should be approximately 79.9. The relay pocket mean (72.8)
should be below the chain average, confirming the 38th percentile claim.

---

## 7. Druggability Criteria Verification

The druggability assessment uses standard Halgren (2009) criteria.
Verify each metric from the pharmacophore data:

```python
import json

with open("results/01-pocket-discovery/hiv1_pharmacophore_data.json") as f:
    data = json.load(f)

sc = data["structural_characterization"]

# B-factor
bf = sc["bfactor_analysis"]["relay_pocket"]
print(f"Pocket B-factor mean: {bf['pocket_mean_all_atoms']} (chain mean: {bf['gp120_chain_mean']})")

# Druggability
dr = sc["druggability_scoring"]
print(f"Volume: {dr['ellipsoid_volume_angstrom3']} A^3")
print(f"Hydrophobic fraction: {dr['hydrophobic_fraction']}")
print(f"Verdict: {dr['verdict']}")

# Conservation
cv = sc["sequence_conservation"]
print(f"Mean entropy: {cv['mean_entropy']}")
print(f"HXB2 identity: {cv['hxb2_identity']}")
for r in cv["per_residue_entropy"]:
    print(f"  {r['aa']}{r['position']}: entropy={r['entropy']} ({r['class']})")

# Geometry
pg = sc["pocket_geometry"]
print(f"Burial ratio: {pg['burial_ratio']}")
print(f"Depth/span: {pg['depth_span_ratio']}")
```

All values should match those reported in the preprint text.

---

## 8. Sequence Conservation Verification

The RPGGGDMR motif can be verified against the HXB2 reference sequence
(GenBank K03455) and the Los Alamos HIV Sequence Database:

1. **HXB2 identity:** Download GenBank K03455 and extract gp120 positions
   461-472. The sequence should be NNTSRPGGGDMR.

2. **Clade conservation:** Visit the Los Alamos HIV Sequence Database
   (https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html),
   download Group M consensus alignments, and compute Shannon entropy
   for positions 465-472. Published analyses consistently show C5 as
   the most conserved gp120 domain (mean entropy ~0.25), with the relay
   segment (0.15) even more conserved than the C5 average.

3. **Quick check:** Search for "RPGGGDMR" in any large HIV-1 Env
   alignment -- the motif should appear with >96% identity across clades.

---

## 9. fpocket Independent Pocket Detection Verification

The fpocket analysis can be reproduced by running fpocket 4.2 on PDB 6MEO:

```bash
# Install fpocket (Linux)
git clone https://github.com/Discngine/fpocket.git
cd fpocket && make && cd ..

# Run on 6MEO
fpocket -f structures/hiv_6meo.pdb
```

Then verify relay residue detection:

```python
import os, json

relay = set(range(465, 473))
pdir = "6MEO_out/pockets"

for fn in sorted(os.listdir(pdir)):
    if not fn.endswith("_atm.pdb"):
        continue
    found = set()
    with open(os.path.join(pdir, fn)) as f:
        for line in f:
            if line[:4] == "ATOM":
                chain = line[21]
                resnum = int(line[22:26])
                if chain == "G" and resnum in relay:
                    found.add((line[17:20].strip(), resnum))
    if found:
        pocket_num = fn.replace("pocket", "").replace("_atm.pdb", "")
        print(f"Pocket {pocket_num}: {sorted(found, key=lambda x: x[1])}")
```

Expected: residues 465-468 in fpocket pocket 5, residues 470-472 in pocket 19,
both spanning chains A (CD4) and G (gp120). Gly469 is expected to be absent
(glycine contributes no pocket-lining side chain). Compare to
`results/01-pocket-discovery/fpocket_analysis.json` for full details including druggability scores and
the maraviroc positive control (fpocket pocket 13, druggability score ~0.97).

---

## 10. FEP Binding Free Energy Verification

The alchemical free energy perturbation (FEP) results can be verified from
the per-replicate analysis JSONs. The thermodynamic cycle is:

```
ΔG_bind = ΔG_solvent − ΔG_complex
```

where ΔG_complex is the decoupling free energy of BMS-626529 in the protein
complex (higher = stronger binding) and ΔG_solvent is the decoupling free
energy in water. Standard convention: negative ΔG_bind = favorable binding.

### Verify binding free energies

```python
import json, statistics

systems = {
    "ternary": [],  # gp120 + CD4 + drug
    "binary": [],   # gp120 + drug (no CD4)
    "solvent": [],  # drug + water
}

for sys_name in systems:
    for rep in [1, 2, 3]:
        with open(f"results/07-free-energy/{sys_name}_rep{rep}_analysis.json") as f:
            data = json.load(f)
        systems[sys_name].append(data["MBAR_decorrelated"]["dG_kcal"])

for sys_name, vals in systems.items():
    mean = statistics.mean(vals)
    sd = statistics.stdev(vals)
    print(f"{sys_name}: MBAR = {mean:.4f} ± {sd:.4f} kcal/mol (reps: {vals})")

# Binding free energies (standard convention)
solv_mean = statistics.mean(systems["solvent"])
bind_binary = solv_mean - statistics.mean(systems["binary"])
bind_ternary = solv_mean - statistics.mean(systems["ternary"])
cd4_effect = bind_ternary - bind_binary

print(f"\nΔG_bind(gp120 only) = {bind_binary:.4f} kcal/mol")
print(f"ΔG_bind(gp120+CD4)  = {bind_ternary:.4f} kcal/mol")
print(f"CD4 effect (ΔΔG)    = {cd4_effect:+.4f} kcal/mol")
```

Expected results:
- ΔG_bind(gp120 only): approximately **−3.86** kcal/mol (favorable binding via MBAR; TI: −4.03)
- ΔG_bind(gp120+CD4): approximately **−0.04** kcal/mol (near-zero binding via MBAR; TI: −0.09)
- CD4 effect: approximately **+3.82** kcal/mol (negative allosteric coupling -- drug selectively stabilizes pre-CD4 conformation; TI: +3.95)

### Verify estimator convergence

All three estimators (TI, BAR, MBAR) should agree within ~0.3 kcal/mol for
each leg. The solvent leg should show <1% coefficient of variation across
replicates (internal consistency control).

### Verify from the master summary

```python
with open("results/07-free-energy/fep_multirep_summary.json") as f:
    summary = json.load(f)

print(json.dumps(summary["binding_free_energy_standard"], indent=2))
print(json.dumps(summary["cd4_effect_on_drug_binding"], indent=2))
```

### Protocol verification

The FEP protocol uses staged decoupling (Coulomb first, then VdW with
soft-core potentials). Verify from the MDP files in `topology/`:

- `couple-lambda0 = vdw-q` (drug fully coupled at λ=0)
- `couple-lambda1 = none` (drug fully decoupled at λ=1)
- 32 λ windows, 5 ns production each
- Soft-core: α=0.5, p=1, σ=0.3 nm
- SD integrator, 300 K, PME, AMBER99SB-ILDN + GAFF2

---

---

## Revision History

- **v1.0 (Sep 1, 2026):** Initial supplementary data package submitted with PCOMPBIOL-D-26-02168.
- **v1.1 (Sep 28, 2026):** `results/08-binding-energetics/mmgbsa_ligand_results.json` updated: added `pdb_residue` and `pdb_domain` fields to per-residue decomposition entries mapping GROMACS sequential numbering to PDB 6MEO chain G numbering. Added `numbering_note` field documenting the piecewise GROMACS→PDB offset (+30 C1, +85 after V1/V2, +93 after V4). Corrected note from "relay pocket residues" to "gp120 residues lining the relay pocket cavity (C1, C2, C3 domains)." SHA256SUMS updated. No changes to numerical values, methods, or conclusions.

---

**Author:** Jonathan Michael Retting, Independent Researcher, Los Angeles, CA
**License:** CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
