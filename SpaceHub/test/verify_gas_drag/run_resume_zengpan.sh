#!/usr/bin/env bash
# run_resume_zengpan.sh – resume timed-out ZengLowEcc-pagn simulations.
#
# Resumes inclinations 45, 120, 135 from their last snapshot for another
# 12 hours. Writes to incl-{inc}-resumed.dat (original files preserved).
#
# Overnight usage (from verify_gas_drag/):
#   nohup ./run_resume_zengpan.sh > run_resume_zengpan.log 2>&1 &
#
# Test usage (kills each sim after 30 s to verify file placement):
#   ./run_resume_zengpan.sh --test

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INCLINATIONS=(45 120 135)
N_WORKERS=3
TEST_MODE=false
TIMEOUT=43200  # 12 hours

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

    timeout "$TIMEOUT" ./"${exe}" "${inc}" > "$logfile" 2>&1 || rc=$?

    # Exit 124 = timeout; treat as OK in test mode (sim was still running)
    if $TEST_MODE && [[ $rc -eq 124 ]]; then
        echo "[OK]   ${label}  i=${inc}°  (timed out as expected in test mode)"
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
g++ -std=c++17 -O3 -pthread resume_ZengLowEcc.cpp -o sim-resume-ZengLowEcc \
    && echo "  sim-resume-ZengLowEcc  OK" || { echo "  sim-resume-ZengLowEcc  FAILED"; exit 1; }
echo "=== Compilation complete ==="

# ── Task list ─────────────────────────────────────────────────────────────────
# Format: "exe|inc|outfile|label"
TASKS=()
for inc in "${INCLINATIONS[@]}"; do
    TASKS+=("sim-resume-ZengLowEcc|${inc}|ZengLowEcc-pagn/incl-${inc}-resumed.dat|ResumeLow")
done

echo ""
echo "=== Submitting ${#TASKS[@]} resumed simulations (${N_WORKERS} workers, timeout=${TIMEOUT}s) ==="

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
echo "=== All resumed simulations finished at $(date) ==="
echo ""

n_ok=$(ls ZengLowEcc-pagn/*-resumed.dat 2>/dev/null | grep -cv FAILED || true)
n_fail=$(ls ZengLowEcc-pagn/*-resumed*_FAILED.dat 2>/dev/null | wc -l || true)
echo "  ZengLowEcc-pagn/ (resumed)  →  ${n_ok} OK,  ${n_fail} FAILED"
if (( n_fail > 0 )); then
    ls ZengLowEcc-pagn/*-resumed*_FAILED.dat 2>/dev/null | sed 's/^/    /'
fi
