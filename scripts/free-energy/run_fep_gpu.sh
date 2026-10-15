#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
# FEP Production Run Script -- GPU-accelerated
# Usage: ./run_fep_gpu.sh <leg> <replica> [start_lambda] [end_lambda]
#
# Same protocol as run_fep.sh but uses GPU offloading.
# Requires GROMACS built with CUDA support.

set -euo pipefail

GMX="${GMX:-gmx}"
NTOMP="${NTOMP:-10}"
GPU_ID="${GPU_ID:-0}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

LEG="${1:?Usage: $0 <leg> <replica> [start_lambda] [end_lambda]}"
REP="${2:?Usage: $0 <leg> <replica> [start_lambda] [end_lambda]}"
START_LAMBDA="${3:-0}"
END_LAMBDA="${4:-31}"

RUNS_DIR="$BASE_DIR/runs"
SYSTEMS_DIR="$BASE_DIR/systems"
MDP_DIR="$BASE_DIR/mdp"
REP_DIR="$RUNS_DIR/$LEG/rep$REP"
LOG_DIR="$BASE_DIR/logs"
mkdir -p "$LOG_DIR"

LOGFILE="$LOG_DIR/${LEG}_rep${REP}_gpu_$(date +%Y%m%d_%H%M%S).log"

GPU_MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu -gpu_id $GPU_ID"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

get_mdp() {
    local phase=$1
    if [ "$LEG" = "solvent" ]; then
        echo "$MDP_DIR/solvent_${phase}.mdp"
    else
        echo "$MDP_DIR/${phase}.mdp"
    fi
}

run_window() {
    local LAM=$1
    local WIN_DIR="$REP_DIR/lambda_$(printf '%02d' $LAM)"
    mkdir -p "$WIN_DIR"

    log "=== Lambda $LAM / $LEG / rep$REP (GPU) ==="

    local COORD TOPOL
    if [ "$LEG" = "solvent" ]; then
        COORD="$SYSTEMS_DIR/solvent/solvent.gro"
        TOPOL="$SYSTEMS_DIR/solvent/topol_solvent.top"
    elif [ "$LEG" = "binary" ]; then
        COORD="$SYSTEMS_DIR/binary/binary_neutral.gro"
        TOPOL="$SYSTEMS_DIR/binary/topol_complex.top"
    else
        COORD="$SYSTEMS_DIR/ternary/md_prod.gro"
        [ -f "$COORD" ] || COORD="$SYSTEMS_DIR/ternary/complex_ions.gro"
        TOPOL="$SYSTEMS_DIR/ternary/topol_complex.top"
    fi

    cd "$WIN_DIR"
    local SEED=$((10000 + RANDOM))

    if [ ! -f em.gro ]; then
        log "  EM..."
        sed "s/LAMBDA_STATE/$LAM/" "$(get_mdp em)" > em.mdp
        $GMX grompp -f em.mdp -c "$COORD" -p "$TOPOL" -o em.tpr -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm em -ntmpi 1 -ntomp "$NTOMP" 2>&1 | tail -3
    fi

    if [ ! -f nvt.gro ]; then
        log "  NVT (GPU)..."
        sed "s/LAMBDA_STATE/$LAM/; s/RANDOM_SEED/$SEED/" "$(get_mdp nvt)" > nvt.mdp
        $GMX grompp -f nvt.mdp -c em.gro -r em.gro -p "$TOPOL" -o nvt.tpr -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm nvt -ntmpi 1 -ntomp "$NTOMP" $GPU_MDRUN_FLAGS 2>&1 | tail -3
    fi

    if [ ! -f npt.gro ]; then
        log "  NPT (GPU)..."
        sed "s/LAMBDA_STATE/$LAM/" "$(get_mdp npt)" > npt.mdp
        $GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p "$TOPOL" -o npt.tpr -t nvt.cpt -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm npt -ntmpi 1 -ntomp "$NTOMP" $GPU_MDRUN_FLAGS 2>&1 | tail -3
    fi

    if [ ! -f prod.gro ]; then
        log "  Production 5 ns (GPU)..."
        sed "s/LAMBDA_STATE/$LAM/" "$(get_mdp prod)" > prod.mdp
        $GMX grompp -f prod.mdp -c npt.gro -p "$TOPOL" -o prod.tpr -t npt.cpt -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm prod -ntmpi 1 -ntomp "$NTOMP" $GPU_MDRUN_FLAGS 2>&1 | tail -3
        log "  Production complete."
    fi

    if [ -f prod.log ]; then
        local PERF
        PERF=$(grep "Performance:" prod.log | tail -1 || echo "N/A")
        log "  Performance: $PERF"
    fi
}

log "Starting FEP (GPU): $LEG rep$REP, lambdas $START_LAMBDA-$END_LAMBDA"
log "GMX: $GMX | GPU: $GPU_ID | Threads: $NTOMP"

for LAM in $(seq "$START_LAMBDA" "$END_LAMBDA"); do
    run_window "$LAM"
done

log "FEP leg $LEG rep$REP complete (GPU)"
