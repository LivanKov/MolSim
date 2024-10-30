#!/bin/bash
# Get the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define the path to the project root and src directory
PROJECT_ROOT="$SCRIPT_DIR/.."
SRC_DIR="$PROJECT_ROOT/src"
CLANG_FORMAT_FILE="$PROJECT_ROOT/.clang-format"

# Ensure .clang-format file exists in the root
if [ ! -f "$CLANG_FORMAT_FILE" ]; then
    echo "Error: .clang-format file not found in the project root."
    exit 1
fi

# Find and format all .cpp and .h files in the src directory
find "$SRC_DIR" -name '*.cpp' -o -name '*.h' | xargs clang-format -i --style=file

echo "All source files in $SRC_DIR have been formatted according to $CLANG_FORMAT_FILE."
