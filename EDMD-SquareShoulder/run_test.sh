#!/bin/bash

# Build script for EDMD-SquareShoulder with all features

# Set colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Building EDMD-SquareShoulder simulation...${NC}"

# Clean and rebuild
make clean
make -j4

if [ $? -ne 0 ]; then
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi

echo -e "${GREEN}Build successful!${NC}"

# Use pos.xyz as the input file
INPUT_FILE="pos.xyz"

if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Input file $INPUT_FILE not found!${NC}"
    exit 1
fi

echo -e "${BLUE}Running simulation with $INPUT_FILE...${NC}"
./edmd_sim --input "$INPUT_FILE" --sim_end_time 1.0 --snapshot_file test_output

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Simulation completed successfully!${NC}"
else
    echo -e "${RED}Simulation failed!${NC}"
fi
