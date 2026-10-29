#!/bin/bash
# Roda o benchmark Newton-CG (CPU x GPU) em uma lista de problemas CUTEst.
#
# Uso:   PRESET=amd|nvidia ./run_benchmark.sh [problems.txt]
#        NOPRE=1 PRESET=nvidia ./run_benchmark.sh      (sem SpMV_preprocess)
# Requer: source <raiz do repositorio>/env_cutest.sh  (CUTEST, SIFDECODE, MASTSIF, MYARCH)
#
# Pode interromper (Ctrl-C) e rodar de novo: cada problema executado e anotado em
# done_problems.txt (nome<TAB>codigo de saida) e e pulado nas proximas execucoes.
# Para repetir um problema, apague a linha dele nesse arquivo.

PRESET="${PRESET:-amd}"
BUILD_DIR="build/${PRESET}"
LIST="${1:-problems.txt}"
OUTPUT_CSV="$PWD/results_newton_cg.csv"
DONE_FILE="done_problems.txt"

for v in CUTEST SIFDECODE MASTSIF MYARCH; do
    if [ -z "${!v}" ]; then
        echo "Erro: variavel $v nao definida. Rode: source <raiz do repositorio>/env_cutest.sh"
        exit 1
    fi
done
if [ ! -f "$LIST" ]; then
    echo "Erro: lista de problemas '$LIST' nao encontrada."
    exit 1
fi

if [ ! -f "$OUTPUT_CSV" ] || [ ! -s "$OUTPUT_CSV" ]; then
    echo "Problem,N,NNZ,Density,NI_CPU,NI_GPU,fcnt_CPU,fcnt_GPU,gcnt_CPU,gcnt_GPU,hcnt_CPU,hcnt_GPU,istop_CPU,istop_GPU,cg_fails_CPU,cg_fails_GPU,first_cg_fail_CPU,first_cg_fail_GPU,last_cg_iter_CPU,last_cg_iter_GPU,t_total_CPU_ms,t_total_GPU_ms,t_avg_CPU_ms,t_avg_GPU_ms,t_cg_CPU_ms,t_cg_GPU_ms,Speedup,Speedup_CG,Sol_Diff,t_GPU_prep_ms,Gnorm_CPU,Gnorm_GPU" > "$OUTPUT_CSV"
fi
mkdir -p "$BUILD_DIR/logs"
touch "$DONE_FILE"

# Parte fixa (modulos, kernels, API generica): compilada uma vez
echo "Compilando a parte fixa (preset $PRESET)..."
CORE_LOG="$BUILD_DIR/logs/_core.log"
if ! { cmake --preset "$PRESET" && cmake --build --preset "$PRESET" --target newtoncore kernels; } > "$CORE_LOG" 2>&1; then
    echo "Erro ao configurar/compilar a parte fixa. Ultimas linhas de $CORE_LOG:"
    tail -n 20 "$CORE_LOG"
    exit 1
fi

args=()
[ -n "$NOPRE" ] && args=(nopre)
trap 'echo; echo "Interrompido. Rode de novo para continuar de onde parou."; exit 130' INT

ran=0; skipped=0; failed=()
for p in $(grep -v '^#' "$LIST"); do
    if cut -f1 "$DONE_FILE" | grep -qx "$p"; then
        skipped=$((skipped + 1))
        continue
    fi
    echo "===== $p ====="
    log="$BUILD_DIR/logs/$p.log"

    if ! cmake --preset "$PRESET" -DPROBLEM="$p" > "$log" 2>&1; then
        echo "[falhou: decodificar/configurar] $p (log: $log)"; tail -n 3 "$log"
        failed+=("$p"); continue
    fi
    if ! cmake --build --preset "$PRESET" --target benchmark_newton_cg >> "$log" 2>&1; then
        echo "[falhou: compilar] $p (log: $log)"; tail -n 3 "$log"
        failed+=("$p"); continue
    fi

    # Limite de tempo removido
    ( cd "$BUILD_DIR" && RESULTS_CSV="$OUTPUT_CSV" ./benchmark_newton_cg "${args[@]}" ) 2>&1 | tee -a "$log"
    rc=${PIPESTATUS[0]}
    [ "$rc" -ne 0 ] && echo "[aviso] $p terminou com codigo $rc"
    printf '%s\t%s\n' "$p" "$rc" >> "$DONE_FILE"
    ran=$((ran + 1))
done

echo "=========================================================="
echo " Executados: $ran | ja feitos (pulados): $skipped | falhas de decodificar/compilar: ${#failed[@]}"
[ "${#failed[@]}" -gt 0 ] && echo " Falharam: ${failed[*]}"
echo " Resultados em $OUTPUT_CSV"
echo "=========================================================="