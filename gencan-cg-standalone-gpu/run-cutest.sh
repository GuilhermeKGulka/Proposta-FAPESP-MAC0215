#!/bin/bash
PRESET="${1:-amd}"
LIST="${2:-problems.txt}"
BUILD_DIR="build/${PRESET}"
OUTPUT_TABLE="table-gencan-gpu.txt"

if [ -z "$CUTEST" ]; then
    echo "Erro: rode 'source ../env_cutest.sh' primeiro."
    exit 1
fi

echo "Compilando a parte fixa para $PRESET..."
cmake --preset "$PRESET" > /dev/null
cmake --build --preset "$PRESET" --target gencan_cpu gencan_gpu kernels > /dev/null

mkdir -p "$BUILD_DIR/logs"
if [ ! -f "$OUTPUT_TABLE" ]; then
    echo "# Problem      SC IERR      CPU(s)         n                        f    gpsupn     iter     fcnt     gcnt     hcnt" > "$OUTPUT_TABLE"
fi

for p in $(grep -v '^#' "$LIST"); do
    echo "===== Compilando $p ====="
    log="$BUILD_DIR/logs/$p.log"
    
    if ! cmake --preset "$PRESET" -DPROBLEM="$p" > "$log" 2>&1; then
        echo "[falhou] $p"; continue
    fi
    if ! cmake --build --preset "$PRESET" --target cutestma >> "$log" 2>&1; then
        echo "[falhou] $p"; continue
    fi

    echo "Rodando $p..."
    ( cd "$BUILD_DIR" && ./cutestma ) 2>&1 | tee -a "$log"
    
    if [ -f "$BUILD_DIR/tabline.txt" ]; then
        cat "$BUILD_DIR/tabline.txt" >> "$OUTPUT_TABLE"
    fi
done
echo "Resultados salvos em $OUTPUT_TABLE"