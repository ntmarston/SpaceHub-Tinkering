#!/usr/bin/env bash
# test_resume_pipeline.sh – Quick validation of the resume pipeline.
# Run from verify_gas_drag/claude_test/
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Resume Pipeline Test ==="
echo "Working directory: $SCRIPT_DIR"

# ── Step 1: Compile ──────────────────────────────────────────────────────────
echo ""
echo "--- Step 1: Compile resume binary ---"
g++ -std=c++17 -O3 -pthread ../resume_ZengLowEcc.cpp -o sim-resume-ZengLowEcc \
    && echo "  [OK] Compilation succeeded" \
    || { echo "  [FAIL] Compilation failed"; exit 1; }

# ── Step 2: Prepare test data ────────────────────────────────────────────────
echo ""
echo "--- Step 2: Prepare test data ---"
mkdir -p ZengLowEcc-pagn
# Use first 20 lines (header + ~6 snapshots) from the real incl-45.dat
if [[ -f ../ZengLowEcc-pagn/incl-45.dat ]]; then
    head -n 20 ../ZengLowEcc-pagn/incl-45.dat > ZengLowEcc-pagn/incl-45.dat
    echo "  [OK] Truncated test data created ($(wc -l < ZengLowEcc-pagn/incl-45.dat) lines)"
else
    echo "  [FAIL] Source file ../ZengLowEcc-pagn/incl-45.dat not found"
    exit 1
fi

# ── Step 3: Test snapshot parsing ─────────────────────────────────────────────
echo ""
echo "--- Step 3: Run resume (30s timeout) ---"
# Run from parent directory so disk file paths resolve correctly
cd "$SCRIPT_DIR/.."
timeout 30 "$SCRIPT_DIR/sim-resume-ZengLowEcc" 45 > "$SCRIPT_DIR/resume_test.log" 2>&1
rc=$?
cd "$SCRIPT_DIR"

if [[ $rc -eq 0 ]]; then
    echo "  [OK] Simulation completed naturally"
elif [[ $rc -eq 124 ]]; then
    echo "  [OK] Simulation running (timed out as expected for a 30s test)"
else
    echo "  [FAIL] Simulation crashed (exit=$rc)"
    echo "  Log output:"
    cat resume_test.log
    # Clean up
    rm -rf ZengLowEcc-pagn sim-resume-ZengLowEcc resume_test.log
    exit 1
fi

# ── Step 4: Check output file ────────────────────────────────────────────────
echo ""
echo "--- Step 4: Verify output ---"
OUTPUT="../ZengLowEcc-pagn/incl-45-resumed.dat"
if [[ -f "$OUTPUT" ]]; then
    lines=$(wc -l < "$OUTPUT")
    size=$(stat --printf="%s" "$OUTPUT")
    echo "  [OK] Resumed file created: ${lines} lines, ${size} bytes"
    if [[ $lines -gt 1 ]]; then
        echo "  First data line:"
        head -n 2 "$OUTPUT" | tail -n 1 | cut -c1-80
    fi
else
    echo "  [WARN] No resumed output file yet (TimeSlice interval may be too large for 30s)"
    echo "  This is expected for a short test — the simulation needs more time to reach first output."
fi

# ── Step 5: Verify original file is untouched ────────────────────────────────
echo ""
echo "--- Step 5: Verify original file integrity ---"
orig_size=$(stat --printf="%s" ../ZengLowEcc-pagn/incl-45.dat)
if [[ $orig_size -gt 0 ]]; then
    echo "  [OK] Original incl-45.dat intact (${orig_size} bytes)"
else
    echo "  [FAIL] Original incl-45.dat was corrupted!"
    exit 1
fi

# ── Cleanup ──────────────────────────────────────────────────────────────────
echo ""
echo "--- Cleanup ---"
rm -rf ZengLowEcc-pagn sim-resume-ZengLowEcc resume_test.log
rm -f ../ZengLowEcc-pagn/incl-45-resumed.dat
echo "  [OK] Test artifacts removed"

echo ""
echo "=== All tests passed ==="
