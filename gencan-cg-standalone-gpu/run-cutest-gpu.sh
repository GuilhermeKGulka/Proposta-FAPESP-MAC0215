#!/bin/bash
# Uso: PRESET=amd ./run-cutest-gpu.sh [problema]
PRESET="${PRESET:-amd}"
BUILD_DIR="build/${PRESET}"
LIST="${1:-problems.txt}"

if [ -z "$CUTEST" ]; then
    echo "Erro: rode 'source ../env_cutest.sh' primeiro."
    exit 1
fi

# Compila a parte fixa
cmake --preset "$PRESET" > /dev/null
cmake --build --preset "$PRESET" --target gencan_cpu gencan_gpu kernels > /dev/null

mkdir -p "$BUILD_DIR/logs"

for p in $(grep -v '^#' "$LIST"); do
    echo "===== Compilando $p ====="
    log="$BUILD_DIR/logs/$p.log"
    
    if ! cmake --preset "$PRESET" -DPROBLEM="$p" > "$log" 2>&1; then
        echo "[falhou] $p"; continue
    fi
    if ! cmake --build --preset "$PRESET" --target cutestma_gpu >> "$log" 2>&1; then
        echo "[falhou] $p"; continue
    fi

    echo "Rodando $p..."
    ( cd "$BUILD_DIR" && ./cutestma_gpu ) 2>&1 | tee -a "$log"
    
    # Salva os resultados no tabline
    if [ -f "$BUILD_DIR/tabline.txt" ]; then
        cat "$BUILD_DIR/tabline.txt" >> "table_gencan_gpu.txt"
    fi
done