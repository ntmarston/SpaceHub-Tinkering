#!/usr/bin/env bash
# run_runtime_demo.sh – compile and run runtime_demo_e0 and runtime_demo_e09
#                       for a sweep of inclinations.
#
# The C++ executables write their own .log files; this script captures
# each process's stdout/stderr to .errorlog files instead.
#
# Overnight usage:
#   cd SpaceHub/test/verify_gas_drag && nohup ./run_runtime_demo.sh > out/RuntimeDemo/run.log 2>&1 &
#
# Test usage (kills each sim after 30 s to verify file placement):
#   ./run_runtime_demo.sh --test

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INCLINATIONS=(0 30 45 70 90 135 180)
N_WORKERS=6
TEST_MODE=false
TIMEOUT=14400 # 8 hours

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
# run_one <exe> <inc> <outfile> <label>
run_one() {
    local exe="$1"
    local inc="$2"
    local outfile="$3"
    local label="$4"
    local errorlog="${outfile%.dat}.errorlog"
    local rc=0

    timeout "$TIMEOUT" "${exe}" "${inc}" > "$errorlog" 2>&1 || rc=$?

    # Exit 124 = timeout
    if [[ $rc -eq 124 ]]; then
        if $TEST_MODE; then
            echo "[OK]   ${label}  i=${inc}°  (timed out as expected in test mode)"
        else
            echo "[TIMEOUT] ${label}  i=${inc}°  (exit=124)"
            if [[ -f "$outfile" ]]; then
                mv "$outfile" "$(dirname "$outfile")/TO_$(basename "$outfile")"
            fi
        fi
        return
    fi

    if [[ $rc -ne 0 ]] || grep -qi "reach max iter" "$errorlog"; then
        echo "[FAIL] ${label}  i=${inc}°  (exit=${rc})  →  ${errorlog}"
        if [[ -f "$outfile" ]]; then
            mv "$outfile" "${outfile%.dat}_FAILED.dat"
        else
            printf "FAILED: simulation produced no output (exit=%d).\nSee: %s\n" \
                "$rc" "$errorlog" > "${outfile%.dat}_FAILED.dat"
        fi
    else
        echo "[OK]   ${label}  i=${inc}°"
    fi
}

# ── Compile ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Compiling ==="
mkdir -p runtime_demo_out/bin
g++ -std=c++17 -O3 -pthread runtime_demo_e0.cpp  -o runtime_demo_out/bin/runtime_demo_e0 \
    && echo "  runtime_demo_e0   OK" || { echo "  runtime_demo_e0   FAILED"; exit 1; }
g++ -std=c++17 -O3 -pthread runtime_demo_e09.cpp -o runtime_demo_out/bin/runtime_demo_e09 \
    && echo "  runtime_demo_e09  OK" || { echo "  runtime_demo_e09  FAILED"; exit 1; }
echo "=== Compilation complete ==="

# ── Ensure output directory exists ───────────────────────────────────────────
mkdir -p out/RuntimeDemo

# ── Task list ─────────────────────────────────────────────────────────────────
# Format: "exe|inc|outfile|label"
TASKS=()
for inc in "${INCLINATIONS[@]}"; do
    TASKS+=("runtime_demo_out/bin/runtime_demo_e0|${inc}|out/RuntimeDemo/incl-${inc}_ecc-0.0.dat|e0")
done
for inc in "${INCLINATIONS[@]}"; do
    TASKS+=("runtime_demo_out/bin/runtime_demo_e09|${inc}|out/RuntimeDemo/incl-${inc}_ecc-0.9.dat|e09")
done

echo ""
echo "=== Submitting ${#TASKS[@]} simulations (${N_WORKERS} workers, timeout=${TIMEOUT}s) ==="

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

n_ok=$(ls out/RuntimeDemo/*.dat 2>/dev/null | grep -cv FAILED || true)
n_fail=$(ls out/RuntimeDemo/*_FAILED.dat 2>/dev/null | wc -l || true)
echo "  out/RuntimeDemo/  →  ${n_ok} OK,  ${n_fail} FAILED"
if (( n_fail > 0 )); then
    ls out/RuntimeDemo/*_FAILED.dat 2>/dev/null | sed 's/^/    /'
fi
