# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
import statistics
import sys, os

os.chdir(sys.argv[1] if len(sys.argv) > 1 else '.')

def p(fn):
    return [float(l.split()[1]) for l in open(fn) if l[0] not in '#@']

print("RMSF Relay Pocket (per residue):")
names = {465:'ARG', 466:'PRO', 467:'GLY', 468:'GLY', 469:'GLY', 470:'ASP', 471:'MET', 472:'ARG'}
for l in open("rmsf_relay.xvg"):
    if l[0] in '#@': continue
    r = int(float(l.split()[0]))
    v = float(l.split()[1])
    nm = names.get(r, 'UNK')
    print(f"  {nm}{r}: {v:.4f} nm = {v*10:.2f} A")

print("\nGyration Relay Pocket:")
v = p("gyrate_relay.xvg")
print(f"  Mean Rg: {statistics.mean(v):.4f} nm = {statistics.mean(v)*10:.2f} A")
print(f"  SD: {statistics.stdev(v):.4f} nm")

print("\nSASA Relay Pocket:")
v = p("sasa_relay.xvg")
print(f"  Mean: {statistics.mean(v):.4f} nm2 = {statistics.mean(v)*100:.2f} A2")
print(f"  SD: {statistics.stdev(v):.4f} nm2")

print("\nH-bonds within Relay Pocket:")
v = p("hbond_relay.xvg")
print(f"  Mean: {statistics.mean(v):.1f}, SD: {statistics.stdev(v):.1f}, Range: {min(v):.0f}-{max(v):.0f}")

print("\nH-bonds Surrounding Protein -> Relay:")
v = p("hbond_surround_relay.xvg")
print(f"  Mean: {statistics.mean(v):.1f}, SD: {statistics.stdev(v):.1f}, Range: {min(v):.0f}-{max(v):.0f}")

print("\nGyration Full Protein:")
v = p("gyrate_protein.xvg")
print(f"  Mean Rg: {statistics.mean(v):.4f} nm = {statistics.mean(v)*10:.2f} A")
print(f"  SD: {statistics.stdev(v):.4f} nm")

print("\nRMSF Top 10 Most Mobile Residues:")
all_rmsf = []
for l in open("rmsf_calpha.xvg"):
    if l[0] in '#@': continue
    parts = l.split()
    all_rmsf.append((int(float(parts[0])), float(parts[1])))
all_rmsf.sort(key=lambda x: -x[1])
for r, v in all_rmsf[:10]:
    print(f"  Res {r}: {v:.4f} nm = {v*10:.2f} A")
print("RMSF Relay in context (percentile):")
relay_vals = [v for r, v in all_rmsf if 465 <= r <= 472]
all_vals = [v for _, v in all_rmsf]
for r in range(465, 473):
    rv = [v for ri, v in all_rmsf if ri == r][0]
    pct = sum(1 for v in all_vals if v <= rv) / len(all_vals) * 100
    nm = names.get(r, 'UNK')
    print(f"  {nm}{r}: {rv:.4f} nm = {rv*10:.2f} A (percentile: {pct:.0f}%)")
