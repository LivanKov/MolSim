#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUTPUT_DIR="$SCRIPT_DIR/../output"
DOCS_DIR = "$SCRIPT_DIR/../docs"

if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"/*
    echo "All files in the output directory have been removed."
else
    echo "Output directory does not exist at: $OUTPUT_DIR"
fi

if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"/*
    echo "All files in the docs directory have been removed."
else
    echo "Docs directory does not exist at: $OUTPUT_DIR"
fi