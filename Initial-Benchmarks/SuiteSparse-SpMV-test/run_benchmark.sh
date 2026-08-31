#!/bin/bash

# Arquivo de resultados em CSV
OUTPUT_CSV="results.csv"

# Diretorio contendo as subpastas dos grupos/matrizes
MM_DIR="MM"

# 1. Verifica se o executavel do benchmark existe
if [ ! -f "./benchmark" ]; then
    echo "Erro: Executavel './benchmark' nao encontrado."
    echo "Compile primeiro com:"
    echo "  hipfc -O3 mmio.f90 benchmark_suitesparse.f90 -lhipsparse -o benchmark"
    exit 1
fi

# 2. Verifica se a pasta MM existe
if [ ! -d "$MM_DIR" ]; then
    echo "Erro: Diretorio '$MM_DIR' nao encontrado no caminho atual."
    exit 1
fi

# 3. Gerencia o cabecalho e inicializa o contador com base nas entradas existentes
if [ ! -f "$OUTPUT_CSV" ] || [ ! -s "$OUTPUT_CSV" ]; then
    echo "Problema,Dim,nnz,Density,niter,t_GPU_copy_ms,t_GPU_SpMV_total_ms,t_GPU_SpMV_avg_ms,t_CPU_SpMV_total_ms,t_CPU_SpMV_avg_ms,Speedup,Speedup_with_copy" > "$OUTPUT_CSV"
    count=0
else
    # Conta as linhas existentes subtraindo o cabecalho
    count=$(tail -n +2 "$OUTPUT_CSV" | wc -l)
fi

echo "=========================================================="
echo " Iniciando varredura recursiva e benchmark na pasta '$MM_DIR'..."
echo " Matrizes ja processadas encontradas no CSV: $count"
echo "=========================================================="

# Busca recursiva por arquivos .mtx dentro de MM/
find -L "$MM_DIR" -type f -name "*.mtx" | sort | while read -r matrix_file; do
    # Nome do arquivo sem a extensao .mtx (ex: dielFilterV2real)
    filename=$(basename "$matrix_file" .mtx)
    
    # Nome do diretorio pai onde o arquivo esta contido
    parent_dir=$(basename "$(dirname "$matrix_file")")
    
    # REGRA DE FILTRAGEM DE ARQUIVOS AUXILIARES:
    if [[ "$filename" =~ _(b|x|z|c)$ ]]; then
        continue
    fi

    # Valida se e a matriz principal pelo padrao esperado
    is_target_matrix=0
    if [ "$filename" = "$parent_dir" ]; then
        is_target_matrix=1
    elif [[ ! "$filename" =~ _[bxz]$ ]]; then
        is_target_matrix=1
    fi

    if [ $is_target_matrix -eq 1 ]; then
        # VERIFICACAO DE DUPLICIDADE NO CSV:
        # Checa se o nome da matriz ja existe isolado na primeira coluna do CSV
        if tail -n +2 "$OUTPUT_CSV" | awk -F',' '{print $1}' | grep -qx "$filename"; then
            echo "[Pulando] Matriz '$filename' ja registrada no CSV."
            continue
        fi

        count=$((count + 1))
        echo "----------------------------------------------------------"
        echo "[$count] Processando: $matrix_file"
        echo "----------------------------------------------------------"
        ./benchmark_cg "$matrix_file"
    fi
done

echo "=========================================================="
echo " Benchmark CG concluído! Resultados salvos em $OUTPUT_CSV"
echo "=========================================================="