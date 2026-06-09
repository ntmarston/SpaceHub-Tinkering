#!/usr/bin/env bash
set -e

OUTDIR="SpaceHub/test/migration_test/AnalyticalValidation/CN08/out"
CN08="SpaceHub/test/migration_test/AnalyticalValidation/CN08"

echo "Clearing old outputs..."
rm -f "$OUTDIR"/*

for fig in 1 2 3; do
    echo "--- Compiling CN08_Fig${fig} ---"
    g++ -std=c++17 -O3 -pthread "$CN08/CN08_Fig${fig}.cpp" -o "$CN08/CN08_Fig${fig}"
    echo "--- Running CN08_Fig${fig} ---"
    "$CN08/CN08_Fig${fig}"
done

echo "Done."