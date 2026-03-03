#!/usr/bin/env bash
# run_overnight_pagn.sh – compile and run pagn-based stae321 + Zeng & Pan simulations.
#
# Replaces the old stae321 (disktab) simulations with pagn versions.
# Retains the Zeng & Pan simulations which use a different inclination array.
#
# Overnight usage (from verify_gas_drag/):
#   nohup ./run_overnight_pagn.sh > run_overnight_pagn.log 2>&1 &
#
# Test usage (terminates each sim after 30 s to check file placement):
#   ./run_overnight_pagn.sh --test

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

# ── Parallel runner ───────────────────────────────────────────────────────
# run_all <executable> <label>
# Runs the executable once per inclination, up to N_WORKERS simultaneously.
run_all() {
    local exe="$1"
    local label="$2"
    local t0=$SECONDS

    echo ""
    echo "=== ${label}: ${#INCLINATIONS[@]} sims × ${N_WORKERS} workers ==="

    for inc in "${INCLINATIONS[@]}"; do
        # Throttle: wait for a slot
        while (( $(jobs -rp | wc -l) >= N_WORKERS )); do
            sleep 0.5
        done

        echo "[${label}] i=${inc}° starting at $(date +%T)"

        if $TEST_MODE; then
            ( timeout 30 ./"${exe}" "${inc}" || true ) &
        else
            ./"${exe}" "${inc}" &
        fi
    done

    wait || true   # wait for all jobs launched in this block
    local elapsed=$(( SECONDS - t0 ))
    printf "=== ${label} done at $(date) — batch took %dh %02dm %02ds ===\n" \
        $(( elapsed/3600 )) $(( (elapsed%3600)/60 )) $(( elapsed%60 ))
}

# ── pagn stae321 simulations ─────────────────────────────────────────────
INCLINATIONS=(5 22 39 56 73 90 107 124 141 158 175)

echo "Clearing previous pagn outputs and executables..."
rm -rf fig11-pagn fig12-pagn stae321-11-pagn stae321-12-pagn
mkdir -p fig11-pagn fig12-pagn

echo "Compiling stae321-11-pagn.cpp..."
g++ -std=c++17 -O3 -pthread stae321-11-pagn.cpp -o stae321-11-pagn

echo "Compiling stae321-12-pagn.cpp..."
g++ -std=c++17 -O3 -pthread stae321-12-pagn.cpp -o stae321-12-pagn

echo "Compilation complete."

run_all "stae321-11-pagn" "Fig11-AeroDrag-pagn"
run_all "stae321-12-pagn" "Fig12-DynFric-pagn"

echo ""
echo "fig11-pagn/ output (${#INCLINATIONS[@]} expected):"
ls -lh fig11-pagn/*.dat 2>/dev/null || echo "  WARNING: no .dat files found in fig11-pagn/"
echo ""
echo "fig12-pagn/ output (${#INCLINATIONS[@]} expected):"
ls -lh fig12-pagn/*.dat 2>/dev/null || echo "  WARNING: no .dat files found in fig12-pagn/"

# ── Zeng & Pan batches ────────────────────────────────────────────────────
INCLINATIONS=(20 45 90 105 120 135 170)

echo ""
echo "=== Zeng & Pan simulations ==="
echo "Clearing previous Zeng outputs and executables..."
rm -rf ZengLowEcc ZengHighEcc ZengLowEccentricity ZengHighEccentricity
mkdir -p ZengLowEcc ZengHighEcc

echo "Compiling ZengLowEccentricity.cpp..."
g++ -std=c++17 -O3 -pthread ZengLowEccentricity.cpp -o ZengLowEccentricity

echo "Compiling ZengHighEccentricity.cpp..."
g++ -std=c++17 -O3 -pthread ZengHighEccentricity.cpp -o ZengHighEccentricity

echo "Zeng compilation complete."

run_all "ZengLowEccentricity"  "ZengLow-DynFric"
run_all "ZengHighEccentricity" "ZengHigh-DynFric"

echo ""
echo "All simulations finished at $(date)."
echo ""
echo "ZengLowEcc/ output (${#INCLINATIONS[@]} expected):"
ls -lh ZengLowEcc/*.dat 2>/dev/null || echo "  WARNING: no .dat files found in ZengLowEcc/"
echo ""
echo "ZengHighEcc/ output (${#INCLINATIONS[@]} expected):"
ls -lh ZengHighEcc/*.dat 2>/dev/null || echo "  WARNING: no .dat files found in ZengHighEcc/"
