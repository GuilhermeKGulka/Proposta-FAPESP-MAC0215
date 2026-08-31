import csv
import ssgetpy
import sys

if len(sys.argv) < 2:
    print("Erro: Nenhum arquivo informado.")
    print("Uso correto: python add_kind_to_csv.py <arquivo_resultados> (SEM .csv)")
    sys.exit(1)

input_csv = sys.argv[1]+".csv"
output_csv = sys.argv[1]+"_with_kinds.csv"

# 1. Carrega o índice completo da SuiteSparse via ssgetpy (limit=0 corrige o aviso anterior)
print("Buscando metadados na SuiteSparse via ssgetpy (isso pode levar alguns segundos)...")
all_matrices = {m.name: m.kind for m in ssgetpy.search(limit=0)}

print(f"Lendo '{input_csv}' e inserindo a coluna 'Kind'...")

try:
    with open(input_csv, mode='r', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        header = next(reader)
        
        # Insere a coluna 'Kind' logo após a coluna 'Problem' (índice 1)
        header.insert(1, 'Kind')
        
        rows = []
        count_sucesso = 0
        count_falha = 0
        
        for row in reader:
            if not row:
                continue  # Pula linhas vazias, se houver
                
            problem_name = row[0].strip()
            kind = "Desconhecido"
            
            # Busca a matriz especificamente pelo nome na API do ssgetpy
            resultado = ssgetpy.search(name=problem_name)
            
            if resultado:
                # Pega o tipo da primeira correspondência encontrada
                kind = resultado[0].kind
                count_sucesso += 1
                print(f"[OK] {problem_name} -> {kind}")
            else:
                print(f"[AVISO] Matriz '{problem_name}' não encontrada no ssgetpy.")
                count_falha += 1
                
            # Insere o valor 'kind' na posição correspondente
            row.insert(1, kind)
            rows.append(row)

    # 2. Escreve os resultados atualizados no novo arquivo CSV
    with open(output_csv, mode='w', encoding='utf-8', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        writer.writerows(rows)

    print("-" * 50)
    print(f"Processo concluído! Arquivo gerado: '{output_csv}'")
    print(f"Matrizes classificadas com sucesso: {count_sucesso}")
    if count_falha > 0:
        print(f"Matrizes não encontradas: {count_falha}")

except FileNotFoundError:
    print(f"Erro: O arquivo '{input_csv}' não foi encontrado no diretório atual.")
except Exception as e:
    print(f"Ocorreu um erro inesperado: {e}")