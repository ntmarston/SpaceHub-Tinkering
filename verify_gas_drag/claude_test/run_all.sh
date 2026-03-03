#!/bin/bash
# Run all diagnostic and isolation tests
# Each test gets a 2-minute wall-clock timeout

set -e
cd "$(dirname "$0")"
mkdir -p output

TIMEOUT=120
CXX="g++ -std=c++17 -O3 -pthread"

echo "==============================="
echo "Phase A: Diagnostic test"
echo "==============================="
echo "Compiling test_diagnostic.cpp..."
$CXX test_diagnostic.cpp -o diag
echo "Running with ${TIMEOUT}s timeout..."
timeout $TIMEOUT ./diag 2>output/diagnostic.log && echo "RESULT: COMPLETED" || {
    code=$?
    if [ $code -eq 124 ]; then
        echo "RESULT: TIMEOUT (hung for ${TIMEOUT}s)"
    else
        echo "RESULT: FAILED (exit code $code)"
    fi
}
echo ""

echo "==============================="
echo "Phase B0: Gravity only (baseline)"
echo "==============================="
echo "Compiling test_b0_gravity_only.cpp..."
$CXX test_b0_gravity_only.cpp -o test_b0
echo "Running..."
timeout $TIMEOUT ./test_b0 && echo "RESULT: COMPLETED" || {
    code=$?
    if [ $code -eq 124 ]; then echo "RESULT: TIMEOUT"; else echo "RESULT: FAILED ($code)"; fi
}
echo ""

echo "==============================="
echo "Phase B1: grad_P = 0 (no sub-Keplerian)"
echo "==============================="
echo "Compiling test_b1_no_gradP.cpp..."
$CXX test_b1_no_gradP.cpp -o test_b1
echo "Running..."
timeout $TIMEOUT ./test_b1 && echo "RESULT: COMPLETED" || {
    code=$?
    if [ $code -eq 124 ]; then echo "RESULT: TIMEOUT"; else echo "RESULT: FAILED ($code)"; fi
}
echo ""

echo "==============================="
echo "Phase B2: No Newton's 3rd law reaction"
echo "==============================="
echo "Compiling test_b2_no_reaction.cpp..."
$CXX test_b2_no_reaction.cpp -o test_b2
echo "Running..."
timeout $TIMEOUT ./test_b2 && echo "RESULT: COMPLETED" || {
    code=$?
    if [ $code -eq 124 ]; then echo "RESULT: TIMEOUT"; else echo "RESULT: FAILED ($code)"; fi
}
echo ""

echo "==============================="
echo "Phase B3: Individual interp() calls"
echo "==============================="
echo "Compiling test_b3_individual_interp.cpp..."
$CXX test_b3_individual_interp.cpp -o test_b3
echo "Running..."
timeout $TIMEOUT ./test_b3 && echo "RESULT: COMPLETED" || {
    code=$?
    if [ $code -eq 124 ]; then echo "RESULT: TIMEOUT"; else echo "RESULT: FAILED ($code)"; fi
}
echo ""

echo "==============================="
echo "Phase B4: Q-dependent vertical density"
echo "==============================="
echo "Compiling test_b4_q_density.cpp..."
$CXX test_b4_q_density.cpp -o test_b4
echo "Running..."
timeout $TIMEOUT ./test_b4 && echo "RESULT: COMPLETED" || {
    code=$?
    if [ $code -eq 124 ]; then echo "RESULT: TIMEOUT"; else echo "RESULT: FAILED ($code)"; fi
}
echo ""

echo "==============================="
echo "All tests complete. Check output/ for data files."
echo "Diagnostic log: output/diagnostic.log"
echo "==============================="
