#!/bin/bash

# Directory paths
input_dir="heuristic-public"
output_file="output.txt"

# Clean up previous output file if it exists
rm -f "$output_file"

# Iterate over the input files
for i in $(seq 1 100); do
    input_file="$input_dir/$i.gr"

    # Measure the time taken by heiCross, kill if it takes longer than 5 minutes (300 seconds)
    { /usr/bin/time -f "%e" timeout 300 ./build/heiCross < "$input_file" > temp_output.txt ; } 2> time_output.txt

    # Check if the program was killed
    if [ $? -eq 124 ]; then
        echo "KILLED" >> "$output_file"
    else
        # Read the output from heiCross
        solution=$(cat temp_output.txt)

        # Read the time taken
        time_taken=$(cat time_output.txt)

        # Append the result to the output file
        echo "$solution,$time_taken" >> "$output_file"
    fi

    # Clean up temporary files
    rm -f temp_output.txt time_output.txt
done
