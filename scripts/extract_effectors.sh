#!/bin/bash

mkdir -p results/effectors

for file in results/prokka/*/*.faa; do
    base=$(basename $file .faa)
    
    # Example: extract secreted proteins (simple motif-based placeholder)
    grep -A1 "signal" $file > results/effectors/${base}_effectors.faa
done
