#!/bin/bash

# Directory paths
input_dir="data/heuristic-public"
output_file="solution.sol"

# Clean up previous output file if it exists
rm -f "$output_file"

# Iterate over the input files
for i in $(seq 1 100); do
    echo $i
    input_file="$input_dir/$i.gr"

    ./runlim --time-limit=60 ./build/heiCross < $input_file > temp.txt

    # Read the output from heiCross
    solution=$(cat temp.txt)

    # Append the result along with the elapsed time to the output file
    echo $i.gr,$solution >> "$output_file"

    # Clean up temporary files
    rm -f temp.txt
done
