#!/bin/bash
# Roda o benchmark Newton-CG em uma lista de problemas CUTEst.
problems=$(grep -v '^#' "${1:-problems.txt}")

OUTPUT_CSV="results_newton_cg.csv"

if [ ! -f "$OUTPUT_CSV" ] || [ ! -s "$OUTPUT_CSV" ]; then
    echo "Problem,N,NNZ,Density,NI_CPU,NI_GPU,fcnt_CPU,fcnt_GPU,gcnt_CPU,gcnt_GPU,hcnt_CPU,hcnt_GPU,istop_CPU,istop_GPU,cg_fails_CPU,cg_fails_GPU,first_cg_fail_CPU,first_cg_fail_GPU,last_cg_iter_CPU,last_cg_iter_GPU,t_total_CPU_ms,t_total_GPU_ms,t_avg_CPU_ms,t_avg_GPU_ms,t_cg_CPU_ms,t_cg_GPU_ms,Speedup,Speedup_CG,Sol_Diff" > "$OUTPUT_CSV"
fi

for p in $problems; do
    echo "===== $p ====="
    make --no-print-directory PROBLEM=$p run
    make --no-print-directory PROBLEM=$p clean
done
echo "Resultados em results_newton_cg.csv"