#!/usr/bin/env bash
# run_overnight.sh – compile and run all stae321 simulations in parallel.
#
# Overnight usage (from verify_gas_drag/):
#   nohup ./run_overnight.sh > run_overnight.log 2>&1 &
#
# Test usage (terminates each sim after 30 s to check file placement):
#   ./run_overnight.sh --test

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INCLINATIONS=(5 22 39 56 73 90 107 124 141 158 175)
N_WORKERS=6
TEST_MODE=false

for arg in "$@"; do
    [[ "$arg" == "--test" ]] && TEST_MODE=true
done

$TEST_MODE && echo "[TEST MODE] Each simulation will be killed after 30 s."
echo "Working directory: $SCRIPT_DIR"
echo "Started: $(date)"

# ── Clean previous run ────────────────────────────────────────────────────────
echo "Clearing previous outputs and executables..."
rm -rf fig11 fig12 stae321-11 stae321-12
mkdir -p fig11 fig12

# ── Compile ───────────────────────────────────────────────────────────────────
echo "Compiling stae321-11.cpp..."
g++ -std=c++17 -O3 -pthread stae321-11.cpp -o stae321-11

echo "Compiling stae321-12.cpp..."
g++ -std=c++17 -O3 -pthread stae321-12.cpp -o stae321-12

echo "Compilation complete."

# ── Parallel runner ───────────────────────────────────────────────────────────
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

run_all "stae321-11" "Fig11-AeroDrag"
run_all "stae321-12" "Fig12-DynFric"

echo ""
echo "fig11/ output (${#INCLINATIONS[@]} expected):"
ls -lh fig11/*.dat 2>/dev/null || echo "  WARNING: no .dat files found in fig11/"
echo ""
echo "fig12/ output (${#INCLINATIONS[@]} expected):"
ls -lh fig12/*.dat 2>/dev/null || echo "  WARNING: no .dat files found in fig12/"

# ── Zeng & Pan batches ────────────────────────────────────────────────────────
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
