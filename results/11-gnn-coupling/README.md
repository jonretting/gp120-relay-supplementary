# 11-gnn-coupling — AE-PocketMiner GNN Allosteric Coupling Prediction

**Paper section:** Methods §2.11, Results §3.21, Table (GNN coupling), Discussion §4.1/§4.2

**Tool:** AE-PocketMiner (Zhang et al. 2026, bioRxiv 10.64898/2026.05.21.726899)
- Attention-enabled GVP-GNN (Geometric Vector Perceptron Graph Neural Network)
- Trained on >3,000 Folding@home MD simulations of 16 non-HIV proteins + SARS-CoV-2 proteome
- Predecessor: PocketMiner (Meller et al. 2023, Nat Commun 14:1177, 263 citations, ROC-AUC 0.87)
- Repository: https://github.com/bowman-lab/ae-pocketminer (MIT license)

**HIV gp120 was absent from the training set — this is a genuine out-of-distribution test.**

## Files

### Analysis Results

| File | Description |
|------|-------------|
| `ev14_aepocketminer_analysis.json` | Combined analysis: per-residue pocket probabilities, attention-weight coupling scores, top-40 coupled residues mapped to functional regions, for all 5 structures |

### Raw Predictions (per structure)

| File | Shape | Description |
|------|-------|-------------|
| `{PDBID}_chainG-preds.npy` | (N,) | Per-residue cryptic pocket probability (0-1) |
| `{PDBID}_chainG-attention_weights.npy` | (1, 2, N, N) | Two-head attention weight matrix (batch=1, heads=2, residues, residues) |

Structures: 6MEO (398 res), 5FYJ (468 res), 5U7O (443 res), 6MU6 (441 res), 4TVP (452 res)

## Key Results

| Structure | Pocket prob. (mean) | #1 coupled residue | Region | Score | gp41 (top 40) | CD4 (top 40) |
|-----------|--------------------:|-------------------:|--------|------:|--------------:|-------------:|
| 6MEO      | 0.509               | 435                | Fostemsavir | 0.033 | 4 | 6 |
| **5FYJ (G471A)** | 0.464       | **487**            | **gp41-contact** | **0.409** | 7 | 4 |
| 5U7O      | 0.317               | 427                | CD4-binding | 0.090 | 7 | 5 |
| 6MU6      | 0.211               | 488                | gp41-contact | 0.113 | 9 | 6 |
| 4TVP      | 0.272               | 481                | gp41-contact | 0.266 | 7 | 5 |

**Critical finding:** G471A variant (5FYJ) attention score 0.409 to gp41 residue 487 = 20x over wildtype baseline. Independent ML confirmation of the steric-gridlock gate mechanism.

## Reproduction

### Step 1: Download RCSB PDB files

```bash
# From the supplementary package root:
for PDB in 6MEO 5FYJ 5U7O 6MU6 4TVP; do
    wget -O structures/${PDB}.pdb https://files.rcsb.org/download/${PDB}.pdb
done
```

### Step 2: Extract chain G

```bash
python scripts/analysis/extract_chain_g_aepocketminer.py
# Creates structures/aepocketminer/{PDBID}_chainG.pdb
```

### Step 3: Run AE-PocketMiner

```bash
# Clone and set up AE-PocketMiner
git clone https://github.com/bowman-lab/ae-pocketminer.git
cd ae-pocketminer

# Create conda environment (see ae-pocketminer/environment.yml)
conda env create -f environment.yml
conda activate aepocketminer

# Create config.yaml pointing to the chain G PDBs
# (see ae-pocketminer README for config format)
python src/xtal_predict.py config.yaml
# Outputs: {PDBID}_chainG-preds.npy and {PDBID}_chainG-attention_weights.npy
```

Docker alternative (CPU, ~9 seconds for all 5 structures):
```bash
docker run --rm -v $(pwd)/structures/aepocketminer:/app/inputs \
    -v $(pwd)/results/11-gnn-coupling:/app/results \
    ae-pocketminer:latest
```

### Step 4: Analyze results

```bash
python scripts/analysis/analyze_aepocketminer_results.py
# Reads preds + attention .npy files, writes ev14_aepocketminer_analysis.json
```

### Software

- AE-PocketMiner: commit `main` (May 2026)
- TensorFlow 2.10 (CPU inference sufficient)
- NumPy (for .npy file reading and analysis)
- No GROMACS or MD software required

### Docker Build (optional)

A Dockerfile for reproducible inference:

```dockerfile
FROM nvidia/cuda:11.2.2-cudnn8-runtime-ubuntu20.04
# Install Miniforge (conda-forge only)
RUN apt-get update && apt-get install -y wget && \
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/conda && rm Miniforge3-Linux-x86_64.sh
ENV PATH="/opt/conda/bin:${PATH}"
# Create environment from ae-pocketminer environment.yml
COPY environment.yml /app/environment.yml
RUN conda env create -f /app/environment.yml
ENV LD_LIBRARY_PATH="/opt/conda/envs/aepocketminer/lib:${LD_LIBRARY_PATH}"
WORKDIR /app
COPY . /app
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "aepocketminer"]
CMD ["python", "src/xtal_predict.py", "/app/config.yaml"]
```

## Attention Matrix Interpretation

The attention weight matrix (N × N) captures learned allosteric coupling between every residue pair. For each relay pocket residue j ∈ {465, ..., 472}, the mean attention weight to each non-pocket residue i is:

    ā_i = (1/8) Σ_j a_ij

where a_ij is the attention weight from relay residue j to residue i, averaged over the 2 attention heads. The top-40 residues ranked by ā_i are mapped to known functional regions (CD4 binding site, gp41-contact surface, fostemsavir binding site).
