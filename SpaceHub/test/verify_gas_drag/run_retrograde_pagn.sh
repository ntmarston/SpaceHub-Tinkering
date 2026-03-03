#!/usr/bin/env bash
# run_retrograde_pagn.sh
#
# Re-runs ONLY the retrograde inclinations (i >= 90°) that were previously
# simulated incorrectly due to the prograde-remap bug. Correct outputs
# (i < 90° for pagn/stae321; i = 20°, 45° for Zeng) are left untouched.
#
# Covers four executables in a single global job pool:
#   stae321-11-pagn  → fig11-pagn/AeroDrag-i{inc}.dat
#   stae321-12-pagn  → fig12-pagn/DynFriction-i{inc}.dat
#   ZengLowEccentricity  → ZengLowEcc/incl-{inc}.dat
#   ZengHighEccentricity → ZengHighEcc/incl-{inc}.dat
#
# Error handling: if a simulation exits non-zero OR its output contains
# "Reach max iter", the .dat file is renamed to *_FAILED.dat and a
# .log file is written alongside it.
#
# Overnight usage (from verify_gas_drag/):
#   nohup ./run_retrograde_pagn.sh > run_retrograde_pagn.log 2>&1 &
#
# Test usage (kills each sim after 30 s to verify file placement):
#   ./run_retrograde_pagn.sh --test

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

N_WORKERS=6
TEST_MODE=false

for arg in "$@"; do
    [[ "$arg" == "--test" ]] && TEST_MODE=true
done

$TEST_MODE && echo "[TEST MODE] Each simulation will be killed after 30 s."
echo "Working directory: $SCRIPT_DIR"
echo "Started: $(date)"

# ── Worker ────────────────────────────────────────────────────────────────────
# run_one <exe> <inc> <outfile> <label>
# Runs ./exe inc, captures output to <outfile base>.log.
# On failure (non-zero exit OR "Reach max iter" in log), renames/creates
# <outfile base>_FAILED.dat and prints a FAIL line; otherwise prints OK.
run_one() {
    local exe="$1"
    local inc="$2"
    local outfile="$3"
    local label="$4"
    local logfile="${outfile%.dat}.log"
    local rc=0

    if $TEST_MODE; then
        timeout 30 ./"${exe}" "${inc}" > "$logfile" 2>&1 || rc=$?
    else
        ./"${exe}" "${inc}" > "$logfile" 2>&1 || rc=$?
    fi

    if [[ $rc -ne 0 ]] || grep -qi "reach max iter" "$logfile"; then
        echo "[FAIL] ${label}  i=${inc}°  (exit=${rc})  →  ${logfile}"
        if [[ -f "$outfile" ]]; then
            mv "$outfile" "${outfile%.dat}_FAILED.dat"
        else
            printf "FAILED: simulation produced no output (exit=%d).\nSee: %s\n" \
                "$rc" "$logfile" > "${outfile%.dat}_FAILED.dat"
        fi
    else
        echo "[OK]   ${label}  i=${inc}°"
    fi
}

# ── Compile ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Compiling ==="
g++ -std=c++17 -O3 -pthread stae321-11-pagn.cpp    -o stae321-11-pagn    \
    && echo "  stae321-11-pagn    OK" || { echo "  stae321-11-pagn    FAILED"; exit 1; }
g++ -std=c++17 -O3 -pthread stae321-12-pagn.cpp    -o stae321-12-pagn    \
    && echo "  stae321-12-pagn    OK" || { echo "  stae321-12-pagn    FAILED"; exit 1; }
g++ -std=c++17 -O3 -pthread ZengLowEccentricity.cpp  -o ZengLowEccentricity  \
    && echo "  ZengLowEccentricity  OK" || { echo "  ZengLowEccentricity  FAILED"; exit 1; }
g++ -std=c++17 -O3 -pthread ZengHighEccentricity.cpp -o ZengHighEccentricity \
    && echo "  ZengHighEccentricity OK" || { echo "  ZengHighEccentricity FAILED"; exit 1; }
echo "=== Compilation complete ==="

# ── Ensure output directories exist ──────────────────────────────────────────
mkdir -p fig11-pagn fig12-pagn ZengLowEcc ZengHighEcc

# ── Task list ─────────────────────────────────────────────────────────────────
# Format: "exe|inc|outfile|label"
# Only the bugged inclinations (i >= 90°) are listed.
TASKS=(
    # stae321-11-pagn (aerodynamic drag, stars)
    "stae321-11-pagn|90|fig11-pagn/AeroDrag-i90.dat|Fig11-AeroDrag-pagn"
    "stae321-11-pagn|107|fig11-pagn/AeroDrag-i107.dat|Fig11-AeroDrag-pagn"
    "stae321-11-pagn|124|fig11-pagn/AeroDrag-i124.dat|Fig11-AeroDrag-pagn"
    "stae321-11-pagn|141|fig11-pagn/AeroDrag-i141.dat|Fig11-AeroDrag-pagn"
    "stae321-11-pagn|158|fig11-pagn/AeroDrag-i158.dat|Fig11-AeroDrag-pagn"
    "stae321-11-pagn|175|fig11-pagn/AeroDrag-i175.dat|Fig11-AeroDrag-pagn"

    # stae321-12-pagn (dynamical friction, BHs)
    "stae321-12-pagn|90|fig12-pagn/DynFriction-i90.dat|Fig12-DynFric-pagn"
    "stae321-12-pagn|107|fig12-pagn/DynFriction-i107.dat|Fig12-DynFric-pagn"
    "stae321-12-pagn|124|fig12-pagn/DynFriction-i124.dat|Fig12-DynFric-pagn"
    "stae321-12-pagn|141|fig12-pagn/DynFriction-i141.dat|Fig12-DynFric-pagn"
    "stae321-12-pagn|158|fig12-pagn/DynFriction-i158.dat|Fig12-DynFric-pagn"
    "stae321-12-pagn|175|fig12-pagn/DynFriction-i175.dat|Fig12-DynFric-pagn"

    # ZengLowEccentricity (e=0.3, dynamical friction)
    "ZengLowEccentricity|90|ZengLowEcc/incl-90.dat|ZengLow-DynFric"
    "ZengLowEccentricity|105|ZengLowEcc/incl-105.dat|ZengLow-DynFric"
    "ZengLowEccentricity|120|ZengLowEcc/incl-120.dat|ZengLow-DynFric"
    "ZengLowEccentricity|135|ZengLowEcc/incl-135.dat|ZengLow-DynFric"
    "ZengLowEccentricity|170|ZengLowEcc/incl-170.dat|ZengLow-DynFric"

    # ZengHighEccentricity (e=0.7, dynamical friction)
    "ZengHighEccentricity|90|ZengHighEcc/incl-90.dat|ZengHigh-DynFric"
    "ZengHighEccentricity|105|ZengHighEcc/incl-105.dat|ZengHigh-DynFric"
    "ZengHighEccentricity|120|ZengHighEcc/incl-120.dat|ZengHigh-DynFric"
    "ZengHighEccentricity|135|ZengHighEcc/incl-135.dat|ZengHigh-DynFric"
    "ZengHighEccentricity|170|ZengHighEcc/incl-170.dat|ZengHigh-DynFric"
)

echo ""
echo "=== Submitting ${#TASKS[@]} simulations (${N_WORKERS} workers) ==="

# ── Single global parallel pool ───────────────────────────────────────────────
for task in "${TASKS[@]}"; do
    # Throttle: wait for a free slot
    while (( $(jobs -rp | wc -l) >= N_WORKERS )); do
        sleep 0.5
    done
    IFS='|' read -r exe inc outfile label <<< "$task"
    echo "  starting ${label}  i=${inc}°  at $(date +%T)"
    run_one "$exe" "$inc" "$outfile" "$label" &
done

wait || true

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo "=== All simulations finished at $(date) ==="
echo ""

for dir in fig11-pagn fig12-pagn ZengLowEcc ZengHighEcc; do
    n_ok=$(ls "${dir}"/*.dat 2>/dev/null | grep -v FAILED | wc -l)
    n_fail=$(ls "${dir}"/*_FAILED.dat 2>/dev/null | wc -l)
    echo "  ${dir}/  →  ${n_ok} OK,  ${n_fail} FAILED"
    if (( n_fail > 0 )); then
        ls "${dir}"/*_FAILED.dat 2>/dev/null | sed 's/^/    /'
    fi
done
