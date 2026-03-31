#!/usr/bin/env bash
# run_migration_tests.sh – compile and run DiskMigration smoketests.
#
# Runs Jimenez, CN06, and Zhu tests in parallel with a 12-hour timeout.
#
# Usage:
#   cd SpaceHub/test/migration_test && nohup ./run_migration_tests.sh > out/run_migration.log 2>&1 &

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TIMEOUT=43200  # 12 hours
TESTS=(JimenezTest CN06Test ZhuTest)

echo "Working directory: $SCRIPT_DIR"
echo "Started: $(date)"

# ── Compile ──────────────────────────────────────────────────────────────────
echo ""
echo "=== Compiling ==="
for test in "${TESTS[@]}"; do
    g++ -std=c++17 -O3 -pthread "${test}.cpp" -o "${test}" \
        && echo "  ${test}  OK" || { echo "  ${test}  FAILED"; exit 1; }
done
echo "=== Compilation complete ==="

mkdir -p out

# ── Worker ───────────────────────────────────────────────────────────────────
run_one() {
    local name="$1"
    local logfile="out/${name}.log"
    local rc=0
    local t_start t_end elapsed

    t_start=$(date +%s)
    timeout "$TIMEOUT" "./${name}" > "$logfile" 2>&1 || rc=$?
    t_end=$(date +%s)
    elapsed=$(( t_end - t_start ))

    if [[ $rc -eq 124 ]]; then
        echo "[TIMEOUT] ${name}  (${elapsed}s, killed after ${TIMEOUT}s)"
    elif [[ $rc -ne 0 ]]; then
        echo "[FAIL]    ${name}  (exit=${rc}, ${elapsed}s)  →  ${logfile}"
    else
        echo "[OK]      ${name}  (${elapsed}s)"
    fi
}

# ── Run all in parallel ──────────────────────────────────────────────────────
echo ""
echo "=== Running ${#TESTS[@]} tests (timeout=${TIMEOUT}s) ==="
for test in "${TESTS[@]}"; do
    echo "  starting ${test} at $(date +%T)"
    run_one "$test" &
done

wait || true

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "=== All tests finished at $(date) ==="
n_ok=$(ls out/*.dat 2>/dev/null | wc -l || true)
echo "  Output files in out/: ${n_ok}"
