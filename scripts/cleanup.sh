#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUTPUT_DIR="$SCRIPT_DIR/../output"

if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"/*
    echo "All files in the output directory have been removed."
else
    echo "Output directory does not exist at: $OUTPUT_DIR"
fi