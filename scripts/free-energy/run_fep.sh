#!/bin/bash
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
# FEP Production Run Script -- CPU-only
# Usage: ./run_fep.sh <leg> <replica> [start_lambda] [end_lambda]
#   e.g.: ./run_fep.sh ternary 1 0 31
#
# Runs the full protocol for each lambda window:
#   1. Energy minimization
#   2. NVT equilibration (100 ps)
#   3. NPT equilibration (500 ps)
#   4. Production (5 ns)

set -euo pipefail

GMX="${GMX:-gmx}"
NTOMP="${NTOMP:-20}"
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

LOGFILE="$LOG_DIR/${LEG}_rep${REP}_$(date +%Y%m%d_%H%M%S).log"

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

    log "=== Lambda $LAM / $LEG / rep$REP ==="

    local COORD TOPOL TOPOL_DIR
    if [ "$LEG" = "solvent" ]; then
        COORD="$SYSTEMS_DIR/solvent/solvent.gro"
        TOPOL="$SYSTEMS_DIR/solvent/topol_solvent.top"
        TOPOL_DIR="$SYSTEMS_DIR/solvent"
    elif [ "$LEG" = "binary" ]; then
        COORD="$SYSTEMS_DIR/binary/binary_neutral.gro"
        TOPOL="$SYSTEMS_DIR/binary/topol_complex.top"
        TOPOL_DIR="$SYSTEMS_DIR/binary"
    else
        COORD="$SYSTEMS_DIR/ternary/md_prod.gro"
        [ -f "$COORD" ] || COORD="$SYSTEMS_DIR/ternary/complex_ions.gro"
        TOPOL="$SYSTEMS_DIR/ternary/topol_complex.top"
        TOPOL_DIR="$SYSTEMS_DIR/ternary"
    fi

    cd "$WIN_DIR"

    local SEED=$((10000 + RANDOM))

    # Step 1: Energy minimization
    if [ ! -f em.gro ]; then
        log "  EM..."
        local EM_MDP=$(get_mdp em)
        sed "s/LAMBDA_STATE/$LAM/" "$EM_MDP" > em.mdp
        $GMX grompp -f em.mdp -c "$COORD" -p "$TOPOL" -o em.tpr -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm em -ntmpi 1 -ntomp "$NTOMP" 2>&1 | tail -3
    else
        log "  EM already done, skipping"
    fi

    # Step 2: NVT equilibration
    if [ ! -f nvt.gro ]; then
        log "  NVT..."
        local NVT_MDP=$(get_mdp nvt)
        sed "s/LAMBDA_STATE/$LAM/; s/RANDOM_SEED/$SEED/" "$NVT_MDP" > nvt.mdp
        $GMX grompp -f nvt.mdp -c em.gro -r em.gro -p "$TOPOL" -o nvt.tpr -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm nvt -ntmpi 1 -ntomp "$NTOMP" 2>&1 | tail -3
    else
        log "  NVT already done, skipping"
    fi

    # Step 3: NPT equilibration
    if [ ! -f npt.gro ]; then
        log "  NPT..."
        local NPT_MDP=$(get_mdp npt)
        sed "s/LAMBDA_STATE/$LAM/" "$NPT_MDP" > npt.mdp
        $GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p "$TOPOL" -o npt.tpr -t nvt.cpt -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm npt -ntmpi 1 -ntomp "$NTOMP" 2>&1 | tail -3
    else
        log "  NPT already done, skipping"
    fi

    # Step 4: Production
    if [ ! -f prod.gro ]; then
        log "  Production (5 ns)..."
        local PROD_MDP=$(get_mdp prod)
        sed "s/LAMBDA_STATE/$LAM/" "$PROD_MDP" > prod.mdp
        $GMX grompp -f prod.mdp -c npt.gro -p "$TOPOL" -o prod.tpr -t npt.cpt -maxwarn 2 2>&1 | tail -3
        $GMX mdrun -deffnm prod -ntmpi 1 -ntomp "$NTOMP" 2>&1 | tail -3
        log "  Production complete."
    else
        log "  Production already done, skipping"
    fi

    if [ -f prod.log ]; then
        local PERF
        PERF=$(grep "Performance:" prod.log | tail -1 || echo "N/A")
        log "  Performance: $PERF"
    fi
}

log "Starting FEP: $LEG rep$REP, lambdas $START_LAMBDA-$END_LAMBDA"
log "GMX: $GMX"
log "Threads: $NTOMP"
log "Base: $BASE_DIR"

for LAM in $(seq "$START_LAMBDA" "$END_LAMBDA"); do
    run_window "$LAM"
done

log "FEP leg $LEG rep$REP complete (lambdas $START_LAMBDA-$END_LAMBDA)"
