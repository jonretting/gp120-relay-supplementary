import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Fix MM-PBSA atom mismatch by matching residues (same order, same count).

Both structures have 105 residues. Compare per-residue atom counts to find
the extra atoms that tleap added at terminal positions.
"""
import parmed
import os

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
os.chdir(WORKDIR)

parm = parmed.load_file('COM.prmtop', '_GMXMMPBSA_COM.inpcrd')
ref = parmed.load_file('_GMXMMPBSA_COM.pdb')

print(f"AMBER: {len(parm.atoms)} atoms, {len(parm.residues)} residues")
print(f"GROMACS: {len(ref.atoms)} atoms, {len(ref.residues)} residues")

# Compare per-residue atom counts
extra_atoms = []
amber_idx = 0

print(f"\n{'Res#':>4} {'AMBER':>20} {'#A':>4} {'GROMACS':>20} {'#G':>4} {'Diff':>5}")
for ri in range(len(parm.residues)):
    ar = parm.residues[ri]
    gr = ref.residues[ri]
    
    a_atoms = list(ar.atoms)
    g_atoms = list(gr.atoms)
    
    diff = len(a_atoms) - len(g_atoms)
    
    if diff != 0:
        print(f"{ri:4d} {ar.name:>4}{ar.number:>5} {len(a_atoms):4d}   {gr.name:>4}{gr.number:>5} {len(g_atoms):4d} {diff:>+5d}")
        
        # The extra atoms are the ones in AMBER but not in GROMACS
        # They are terminal patch atoms (H1, H2, H3 at N-term, OXT at C-term)
        # Find them by name: atoms in AMBER residue not present in GROMACS residue
        gro_names = set(a.name for a in g_atoms)
        for a in a_atoms:
            if a.name not in gro_names:
                extra_atoms.append(a.idx + 1)  # 1-indexed
                print(f"      -> extra: {a.name} (idx {a.idx+1})")
            else:
                # Remove from set to handle duplicate names
                gro_names.discard(a.name)

    amber_idx += len(a_atoms)

print(f"\nTotal extra atoms: {len(extra_atoms)}")
print(f"Expected: 37, Found: {len(extra_atoms)}")

if len(extra_atoms) > 0:
    # Strip from complex
    mask = '@' + ','.join(str(i) for i in extra_atoms)
    com = parmed.load_file('COM.prmtop', '_GMXMMPBSA_COM.inpcrd')
    com.strip(mask)
    com.save('COM_fixed.prmtop', overwrite=True)
    print(f"\nCOM_fixed.prmtop: {len(com.atoms)} atoms")
    
    # Find receptor/ligand split
    # From gmx_MMPBSA: receptor = chain G (gp120), ligand = chain A (CD4)
    # Receptor is first 82 residues, ligand is last 23
    REC_NRES = 82
    
    # Receptor extra atoms (only those in first 82 residues)
    rec_parm = parmed.load_file('REC.prmtop', '_GMXMMPBSA_REC.inpcrd')
    rec_ref = parmed.load_file('_GMXMMPBSA_REC.pdb')
    print(f"\nREC: {len(rec_parm.atoms)} -> {len(rec_ref.atoms)}")
    
    rec_extra = []
    for ri in range(len(rec_parm.residues)):
        ar = rec_parm.residues[ri]
        gr = rec_ref.residues[ri]
        gro_names = set(a.name for a in gr.atoms)
        for a in ar.atoms:
            if a.name not in gro_names:
                rec_extra.append(a.idx + 1)
            else:
                gro_names.discard(a.name)
    
    if rec_extra:
        rec_parm.strip('@' + ','.join(str(i) for i in rec_extra))
    rec_parm.save('REC_fixed.prmtop', overwrite=True)
    print(f"REC_fixed.prmtop: {len(rec_parm.atoms)} atoms")
    
    # Ligand
    lig_parm = parmed.load_file('LIG.prmtop', '_GMXMMPBSA_LIG.inpcrd')
    lig_ref = parmed.load_file('_GMXMMPBSA_LIG.pdb')
    print(f"\nLIG: {len(lig_parm.atoms)} -> {len(lig_ref.atoms)}")
    
    lig_extra = []
    for ri in range(len(lig_parm.residues)):
        ar = lig_parm.residues[ri]
        gr = lig_ref.residues[ri]
        gro_names = set(a.name for a in gr.atoms)
        for a in ar.atoms:
            if a.name not in gro_names:
                lig_extra.append(a.idx + 1)
            else:
                gro_names.discard(a.name)
    
    if lig_extra:
        lig_parm.strip('@' + ','.join(str(i) for i in lig_extra))
    lig_parm.save('LIG_fixed.prmtop', overwrite=True)
    print(f"LIG_fixed.prmtop: {len(lig_parm.atoms)} atoms")
    
    # Verify
    c = parmed.load_file('COM_fixed.prmtop')
    r = parmed.load_file('REC_fixed.prmtop')
    l = parmed.load_file('LIG_fixed.prmtop')
    print(f"\n=== FINAL VERIFICATION ===")
    print(f"COM: {len(c.atoms)} (target: 1715)")
    print(f"REC: {len(r.atoms)} (target: 1322)")
    print(f"LIG: {len(l.atoms)} (target: 393)")
    print(f"REC+LIG: {len(r.atoms)+len(l.atoms)} (target: 1715)")
    
    if len(c.atoms) == 1715:
        print("\nSUCCESS - Ready for MMPBSA.py!")
