#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="$SCRIPT_DIR/.."
SRC_DIR="$PROJECT_ROOT/src"
TEST_DIR="$PROJECT_ROOT/tests"
CLANG_FORMAT_FILE="$PROJECT_ROOT/.clang-format"

if [ ! -f "$CLANG_FORMAT_FILE" ]; then
    echo "Error: .clang-format file not found in the project root."
    exit 1
fi

find "$SRC_DIR" -name '*.cpp' -o -name '*.h' | xargs clang-format -i --style=file

echo "All source files in $SRC_DIR have been formatted according to $CLANG_FORMAT_FILE."

find "$TEST_DIR" -name '*.cc' | xargs clang-format -i --style=file

echo "All source files in $TEST_DIR have been formatted according to $CLANG_FORMAT_FILE."
