#!/usr/bin/env bash
# run_clocktests.sh – compile and run all 6 clock-test integrator methods.
#
# Overnight usage (from clocking_tests/):
#   nohup ./run_clocktests.sh > run_clocktests.log 2>&1 &
#
# Test usage (kills each sim after 30 s to verify it starts):
#   ./run_clocktests.sh --test

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

METHODS=(AR_Radau BS Radau Sym6 Sym10 Sym6_Plus)
INCLINATION=5.0
N_WORKERS=3
TEST_MODE=false
TIMEOUT=57600  # 16 hours

for arg in "$@"; do
    [[ "$arg" == "--test" ]] && TEST_MODE=true
done

if $TEST_MODE; then
    TIMEOUT=30
    echo "[TEST MODE] Each simulation will be killed after 30 s."
fi
echo "Working directory: $SCRIPT_DIR"
echo "Started: $(date)"

# ── Worker ────────────────────────────────────────────────────────────────────
run_one() {
    local method="$1"
    local exe="ct_${method}"
    local rc=0

    echo "  [START] ${method}  at $(date +%T)"
    timeout "$TIMEOUT" ./"${exe}" "${INCLINATION}" || rc=$?

    # Exit 124 = timeout; treat as OK in test mode (sim was still running)
    if $TEST_MODE && [[ $rc -eq 124 ]]; then
        echo "  [OK]    ${method}  (timed out as expected in test mode)"
        return
    fi

    if [[ $rc -eq 124 ]]; then
        echo "  [TIMEOUT] ${method}  (>${TIMEOUT}s)"
    elif [[ $rc -ne 0 ]]; then
        echo "  [FAIL]  ${method}  (exit=${rc})"
    else
        echo "  [OK]    ${method}"
    fi
}

# ── Compile ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Compiling ==="
for method in "${METHODS[@]}"; do
    g++ -std=c++17 -O3 -pthread "inclined_${method}.cpp" -o "ct_${method}" \
        && echo "  ct_${method}    OK" \
        || { echo "  ct_${method}    FAILED"; exit 1; }
done
echo "=== Compilation complete ==="

# ── Submit tasks (3 parallel workers) ─────────────────────────────────────────
echo ""
echo "=== Submitting ${#METHODS[@]} simulations (${N_WORKERS} workers, timeout=${TIMEOUT}s) ==="

for method in "${METHODS[@]}"; do
    # Throttle: wait for a free worker slot
    while (( $(jobs -rp | wc -l) >= N_WORKERS )); do
        sleep 0.5
    done
    run_one "$method" &
done

wait || true

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo "=== All simulations finished at $(date) ==="
echo ""

# C++ ostringstream strips trailing .0: "5.0" -> "5"
INCL_STR=$(printf "%g" "$INCLINATION")

for method in "${METHODS[@]}"; do
    dat_file="${method}-clocktest-${INCL_STR}.dat"
    log_file="${method}-clocktest-${INCL_STR}.log"
    if [[ -f "$dat_file" ]]; then
        elapsed=$(grep "^STOP" "$log_file" 2>/dev/null \
                  | grep -o 'elapsed=[0-9.]*' | cut -d= -f2 || echo "?")
        echo "  ${method}: OK  (${elapsed}s)"
    else
        echo "  ${method}: FAILED or incomplete  (no .dat file)"
    fi
done
