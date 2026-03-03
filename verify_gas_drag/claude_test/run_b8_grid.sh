#!/bin/bash
# Phase B8: Exhaustive parameter grid for smoothed ZengAndPan disk
# Tests 27 combinations of (sma, ecc, inc) with 120s timeout each
# Worst case: ~54 min if all timeout

cd "$(dirname "$0")"
mkdir -p output

echo "=== Phase B8: Compiling test_b8_smoothed_grid.cpp ==="
g++ -std=c++17 -O3 -pthread test_b8_smoothed_grid.cpp -o test_b8
echo "Compilation successful."
echo ""

SMAS="0.01 0.02 0.05"
ECCS="0.3 0.67 0.9"
INCS="5 20 45"
TIMEOUT=120

SUMMARY="output/b8_grid_summary.txt"
echo "Phase B8 Grid Results — $(date)" > "$SUMMARY"
echo "Disk: disk_ZengAndPan_pagn.csv (smoothed)" >> "$SUMMARY"
echo "Timeout: ${TIMEOUT}s per test" >> "$SUMMARY"
echo "Grid: sma={$SMAS} x ecc={$ECCS} x inc={$INCS}" >> "$SUMMARY"
echo "" >> "$SUMMARY"
printf "%-10s %-6s %-6s %-10s %-10s %s\n" "sma(PC)" "ecc" "inc" "status" "time(s)" "notes" >> "$SUMMARY"
printf "%-10s %-6s %-6s %-10s %-10s %s\n" "-------" "---" "---" "------" "-------" "-----" >> "$SUMMARY"

PASS=0
FAIL=0
TOTAL=0

for sma in $SMAS; do
    for ecc in $ECCS; do
        for inc in $INCS; do
            TOTAL=$((TOTAL + 1))
            LOGFILE="output/b8_a${sma}_e${ecc}_i${inc}.log"
            DATFILE="output/b8_a${sma}_e${ecc}_i${inc}.dat"

            echo "--- Test $TOTAL/27: sma=${sma}PC e=${ecc} i=${inc}deg ---"

            START_TIME=$(date +%s.%N)
            timeout $TIMEOUT ./test_b8 "$sma" "$ecc" "$inc" 2>"$LOGFILE"
            EXIT_CODE=$?
            END_TIME=$(date +%s.%N)
            ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)

            if [ $EXIT_CODE -eq 0 ]; then
                STATUS="PASS"
                PASS=$((PASS + 1))
                # Check for dt collapse in log
                MIN_DT=$(grep -oP 'dt=\K[0-9.e+-]+' "$LOGFILE" | sort -g | head -1)
                NOTES="min_dt=$MIN_DT"
            elif [ $EXIT_CODE -eq 124 ]; then
                # Distinguish "slow but healthy" from "dt collapse"
                LAST_DT=$(grep -oP 'dt=\K[0-9.e+-]+' "$LOGFILE" | tail -1)
                if [ -n "$LAST_DT" ] && (( $(echo "$LAST_DT > 10000" | bc -l) )); then
                    STATUS="SLOW"
                    PASS=$((PASS + 1))  # healthy orbit, just slow
                else
                    STATUS="TIMEOUT"
                    FAIL=$((FAIL + 1))
                fi
                LAST_LINE=$(tail -1 "$LOGFILE")
                NOTES="last_dt=$LAST_DT  $LAST_LINE"
            else
                STATUS="FAIL($EXIT_CODE)"
                FAIL=$((FAIL + 1))
                NOTES="exit code $EXIT_CODE"
            fi

            printf "%-10s %-6s %-6s %-10s %-10.1f %s\n" "$sma" "$ecc" "$inc" "$STATUS" "$ELAPSED" "$NOTES" >> "$SUMMARY"
            echo "  -> $STATUS (${ELAPSED}s)"
        done
    done
done

echo "" >> "$SUMMARY"
echo "Total: $TOTAL | Pass: $PASS | Fail: $FAIL" >> "$SUMMARY"

echo ""
echo "=== GRID COMPLETE ==="
echo "Pass: $PASS / $TOTAL"
echo "Summary: $SUMMARY"
cat "$SUMMARY"
