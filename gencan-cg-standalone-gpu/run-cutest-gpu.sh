#!/bin/bash
# Uso: source ../env_cutest.sh ; PRESET=amd ./run_cutest_gpu.sh [lista.txt]
PRESET="${PRESET:-amd}"
BUILD="build/${PRESET}"
LIST="${1:-problems.txt}"
mkdir -p "$BUILD/logs"
[ -z "$CUTEST" ] && { echo "rode: source ../env_cutest.sh"; exit 1; }

HDR="# Problem      SC IERR      CPU(s)         n                        f    gpsupn     iter     fcnt     gcnt     hcnt"
[ -f table_cpu.txt ] || echo "$HDR" > table_cpu.txt
[ -f table_gpu.txt ] || echo "$HDR" > table_gpu.txt
touch done_gencan.txt

cmake --preset "$PRESET" > "$BUILD/logs/_core.log" 2>&1 && \
cmake --build --preset "$PRESET" --target kernels gencan_gpu >> "$BUILD/logs/_core.log" 2>&1 || \
  { echo "falha na parte fixa"; tail -20 "$BUILD/logs/_core.log"; exit 1; }

for p in $(grep -v '^#' "$LIST"); do
  grep -qx "$p" done_gencan.txt && continue
  echo "===== $p ====="
  log="$BUILD/logs/$p.log"
  if ! cmake --preset "$PRESET" -DPROBLEM="$p" > "$log" 2>&1 || \
     ! cmake --build --preset "$PRESET" --target cutestma_cpu cutestma_gpu >> "$log" 2>&1; then
    echo "[falhou compilar] $p"; continue
  fi
  for v in cpu gpu; do
    rm -f "$BUILD/tabline.txt" "$BUILD/solver-interrupted-tabline.txt"
    ( cd "$BUILD" && ulimit -St 600 && ./cutestma_$v ) >> "$log" 2>&1
    if   [ -f "$BUILD/tabline.txt" ]; then cat "$BUILD/tabline.txt" >> "table_$v.txt"
    else printf " %-10s  FALHOU/TIMEOUT\n" "$p" >> "table_$v.txt"; fi
  done
  echo "$p" >> done_gencan.txt
done