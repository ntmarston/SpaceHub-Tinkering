#!/bin/bash
# Regression test script: compile, run 3 tests, compare against baselines
set -e

echo "=== Compiling ==="
g++ -std=c++17 -O3 -pthread SpaceHub/test/nick_test/disk_model/KeplerSimple-DiskModel.cpp -o test_dm1
g++ -std=c++17 -O3 -pthread SpaceHub/test/nick_test/disk_model/KeplerSimple-DiskModelSelfReg.cpp -o test_dm2
g++ -std=c++17 -O3 -pthread SpaceHub/test/nick_test/disk_model/KeplerSimple-DiskModeRetrograde.cpp -o test_dm3

echo "=== Test 1: Standard zone, low ecc (10k yr) ==="
timeout 60 ./test_dm1
cp SpaceHub/test/nick_test/disk_model/AnalyticalTests/DiskModel-Retrograde.dat /tmp/current_test1.dat

echo "=== Test 2: Self-reg zone (1M yr) ==="
timeout 120 ./test_dm2

echo "=== Test 3: Standard zone, high ecc (1M yr) ==="
timeout 120 ./test_dm3
cp SpaceHub/test/nick_test/disk_model/AnalyticalTests/DiskModel-Retrograde.dat /tmp/current_test3.dat

echo ""
echo "=== Comparing against baselines ==="
PASS=true
for i in 1 2 3; do
    if [ "$i" = "2" ]; then
        CURRENT="SpaceHub/test/nick_test/disk_model/AnalyticalTests/DiskModel-SelfReg.dat"
    else
        CURRENT="/tmp/current_test${i}.dat"
    fi
    BASELINE="/tmp/baseline_test${i}.dat"
    if diff -q "$CURRENT" "$BASELINE" > /dev/null 2>&1; then
        echo "Test $i: IDENTICAL"
    else
        echo "Test $i: DIFFERS"
        diff "$CURRENT" "$BASELINE" | head -5
        PASS=false
    fi
done

if $PASS; then
    echo "ALL TESTS PASS (bit-identical)"
else
    echo "SOME TESTS DIFFER - check output above"
fi
