#!/bin/bash
PRESET="${1:-amd}"

echo "Compilando com preset: $PRESET"
cmake --preset "$PRESET"
cmake --build --preset "$PRESET" --target toyprob cutestma -j

echo "Pronto. Os executaveis estao em build/$PRESET/"