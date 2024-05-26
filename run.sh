#!/bin/bash

# Directory paths
input_dir="heuristic-public"
output_file="output.txt"

# Clean up previous output file if it exists
rm -f "$output_file"

# Iterate over the input files
for i in $(seq 1 100); do
    input_file="$input_dir/$i.gr"

    # Start the timer
    start_time=$(($(date +%s%N)/1000000)) # Convert nanoseconds to milliseconds

    # Run the algorithm and capture its PID
    ./build/heiCross < "$input_file" > temp_output.txt &
    pid=$!

    # Start the timer for 5 minutes (300 seconds)
    (sleep 300 && kill -s SIGTERM $pid) &

    # Wait for the process to finish
    wait $pid

    # Calculate the elapsed time
    end_time=$(($(date +%s%N)/1000000)) # Convert nanoseconds to milliseconds
    elapsed_time=$((end_time - start_time))

    # Read the output from heiCross
    solution=$(cat temp_output.txt)

    # Append the result along with the elapsed time to the output file
    echo "$solution,$elapsed_time" >> "$output_file"

    # Clean up temporary files
    rm -f temp_output.txt
done
