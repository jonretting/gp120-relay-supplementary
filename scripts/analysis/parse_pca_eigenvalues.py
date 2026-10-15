#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse PCA eigenvalues from gmx anaeig output into JSON."""
import json, sys, math

pdb_id = sys.argv[1] if len(sys.argv) > 1 else "UNKNOWN"
schlitter = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

vals = []
with open("eigenvalues_pocket.xvg") as f:
    for line in f:
        if line.startswith(('#', '@')):
            continue
        parts = line.split()
        if len(parts) >= 2:
            vals.append(float(parts[1]))

total = sum(vals) if vals else 1.0
pcs = []
cum = 0
for i, v in enumerate(vals[:10]):
    pct = 100.0 * v / total
    cum += pct
    pcs.append({
        "pc": i + 1,
        "eigenvalue_nm2": round(v, 4),
        "pct_variance": round(pct, 1),
        "cumulative_pct": round(cum, 1)
    })

result = {
    "pdb_id": pdb_id,
    "method": "PCA via gmx covar + gmx anaeig on relay pocket, 500ns, dt=50ps",
    "total_variance_nm2": round(total, 4),
    "top_3_pct": round(sum(p["pct_variance"] for p in pcs[:3]), 1),
    "top_5_pct": round(sum(p["pct_variance"] for p in pcs[:5]), 1),
    "schlitter_entropy_J_mol_K": schlitter,
    "principal_components": pcs
}

outfile = f"{pdb_id}_500ns_pca.json"
with open(outfile, "w") as f:
    json.dump(result, f, indent=2)

print(json.dumps(result, indent=2))
