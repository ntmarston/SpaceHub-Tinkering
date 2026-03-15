#!/usr/bin/env bash
# run_zengpan_reprod.sh – compile and run Zeng & Pan reproduction simulations.
#
# Reproduces figures 9 (e=0.3) and 10 (e=0.7) from arXiv:2601.11925v2
# using SG_01Edd.csv (?)
#
# Overnight usage:
#   cd SpaceHub/test/verify_gas_drag && nohup ./run_zengpan_reprod.sh > out/run_zengpan.log 2>&1 &
#
# Test usage (kills each sim after 30 s to verify file placement):
#   ./run_zengpan_reprod.sh --test

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

#INCLINATIONS=(90 105 120 135)
INCLINATIONS=(120 135)
N_WORKERS=2
TEST_MODE=false
TIMEOUT=115200 # 32 hours

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
    local logfile="${outfile%.dat}.log"
    local rc=0

    timeout "$TIMEOUT" "${exe}" "${inc}" > "$logfile" 2>&1 || rc=$?

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
# g++ -std=c++17 -O3 -pthread ZengLowEcc.cpp -o simulations/bin/ZengLowEcc \
#     && echo "  ZengLowEcc    OK" || { echo "  ZengLowEcc    FAILED"; exit 1; }
g++ -std=c++17 -O3 -pthread ZengHighEcc.cpp -o simulations/bin/ZengHighEcc \
    && echo "  ZengHighEcc   OK" || { echo "  ZengHighEcc   FAILED"; exit 1; }
echo "=== Compilation complete ==="

# ── Ensure output directories exist ──────────────────────────────────────────
# mkdir -p out/ZengLowEcc out/ZengHighEcc
mkdir -p out/ZengHighEcc

# ── Task list ─────────────────────────────────────────────────────────────────
# Format: "exe|inc|outfile|label"
TASKS=()
# for inc in "${INCLINATIONS[@]}"; do
#     TASKS+=("simulations/bin/ZengLowEcc|${inc}|out/ZengLowEcc/incl-${inc}.dat|ZengLow")
# done
for inc in "${INCLINATIONS[@]}"; do
    TASKS+=("simulations/bin/ZengHighEcc|${inc}|out/ZengHighEcc/incl-${inc}.dat|ZengHigh")
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

for dir in out/ZengHighEcc; do  # out/ZengLowEcc (to uncomment put before out/ZengHighEcc on this line)
    n_ok=$(ls "${dir}"/*.dat 2>/dev/null | grep -cv FAILED || true)
    n_fail=$(ls "${dir}"/*_FAILED.dat 2>/dev/null | wc -l || true)
    echo "  ${dir}/  →  ${n_ok} OK,  ${n_fail} FAILED"
    if (( n_fail > 0 )); then
        ls "${dir}"/*_FAILED.dat 2>/dev/null | sed 's/^/    /'
    fi
done
