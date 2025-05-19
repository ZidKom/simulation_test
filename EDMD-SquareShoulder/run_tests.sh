#!/bin/bash

# Compile and run the simulation with various test cases

# Set up working directory
cd /home/zid/simulation_test/EDMD-SquareShoulder

# Run all tests using pos.xyz as input
INPUT_FILE="pos.xyz"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file $INPUT_FILE not found!"
    exit 1
fi

# Build the project
make clean && make -j4
if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

# Test 1: Basic simulation with default parameters
echo "Running basic simulation test..."
./edmd_sim -n 100 -t 10.0 -L 10.0 --snapshot_file snapshot_basic

# Test 2: LAMMPS input format test (if input file exists)
if [ -f "test_data.lammps" ]; then
    echo "Running LAMMPS input format test..."
    ./edmd_sim -i test_data.lammps --lammps -t 5.0 --snapshot_file snapshot_lammps
fi

# Test 3: Benchmark mode
echo "Running benchmark test..."
./edmd_sim -n 200 -t 20.0 -L 15.0 -b

# Test 4: Validation mode (if reference file exists)
if [ -f "reference.snap" ]; then
    echo "Running validation test..."
    ./edmd_sim -n 100 -t 10.0 -L 10.0 -v -r reference.snap --snapshot_file snapshot_validation
fi

# Test 5: Square shoulder potential test
echo "Running square shoulder potential test..."
./edmd_sim -n 50 -t 5.0 -L 8.0 --shoulder_U 1.0 --shoulder_sig 0.3 --snapshot_file snapshot_shoulder

# Run the main simulation test
./edmd_sim --input "$INPUT_FILE" --sim_end_time 1.0 --snapshot_file test_output
if [ $? -eq 0 ]; then
    echo "Simulation completed successfully!"
else
    echo "Simulation failed!"
fi

echo "All tests completed!"
