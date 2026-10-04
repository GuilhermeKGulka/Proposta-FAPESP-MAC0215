import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Nome do arquivo de log de texto de saída
arquivo_txt = 'relatorio_analise_cg.txt'

# Classe auxiliar para escrever simultaneamente no terminal e no arquivo .txt
class Tee(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, 'w', encoding='utf-8')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

# Redireciona a saída padrão
sys.stdout = Tee(arquivo_txt)

try:
    # 1. Carregamento dos dados
    arquivo_csv = sys.argv[1]
    df = pd.read_csv(arquivo_csv)

    # 2. Identificação dinâmica das colunas essenciais
    istop_cpu_col = [c for c in df.columns if 'istop_cpu' in c.lower()][0]
    istop_gpu_col = [c for c in df.columns if 'istop_gpu' in c.lower()][0]
    speedup_col = [c for c in df.columns if c.lower() == 'speedup'][0]
    speedup_copy_col = [c for c in df.columns if 'speedup' in c.lower() and 'copy' in c.lower()][0]

    # 3. Classificação de Status de Convergência
    cpu_conv = (df[istop_cpu_col] == 0)
    gpu_conv = (df[istop_gpu_col] == 0)

    df.loc[cpu_conv & gpu_conv, 'Status'] = 'Ambas Convergiram'
    df.loc[cpu_conv & ~gpu_conv, 'Status'] = 'Apenas CPU Convergiu'
    df.loc[~cpu_conv & gpu_conv, 'Status'] = 'Apenas GPU Convergiu'
    df.loc[~cpu_conv & ~gpu_conv, 'Status'] = 'Nenhuma Convergiu'

    # Diferença relativa percentual de iterações: ((GPU - CPU) / CPU) * 100
    df['iter_diff_pct'] = ((df['CG_Iter_GPU'] - df['CG_Iter_CPU']) / df['CG_Iter_CPU']) * 100

    total_matrizes = len(df)
    contagem_status = df['Status'].value_counts()
    percentuais = (contagem_status / total_matrizes) * 100

    # Subconjunto estrito onde ambas convergem
    df_conv = df[df['Status'] == 'Ambas Convergiram'].copy()
    total_conv = len(df_conv)

    # ==========================================
    # RELATÓRIOS (Salvos no .txt e exibidos no terminal)
    # ==========================================
    print("=" * 70)
    print(f"ANÁLISE DE CONVERGÊNCIA GLOBAL (Arquivo: {arquivo_csv} | Total: {total_matrizes})")
    print("=" * 70)

    categorias = ['Ambas Convergiram', 'Apenas CPU Convergiu', 'Apenas GPU Convergiu', 'Nenhuma Convergiu']
    for cat in categorias:
        cont = contagem_status.get(cat, 0)
        perc = percentuais.get(cat, 0.0)
        print(f"{cat:<25}: {cont:>4} ({perc:>5.1f}%)")

    print("\n" + "=" * 70)
    print(f"ESTATÍSTICAS DE SPEEDUP (Filtrado: Apenas {total_conv} matrizes onde AMBAS convergiram)")
    print("=" * 70)

    if not df_conv.empty:
        stats_puro = df_conv[speedup_col].describe()
        stats_copia = df_conv[speedup_copy_col].describe()
        
        print(f"Speedup Puro    - Média: {stats_puro['mean']:>6.2f}x | Mediana: {stats_puro['50%']:>6.2f}x | Máx: {stats_puro['max']:>6.2f}x")
        print(f"Speedup c/ Cópia - Média: {stats_copia['mean']:>6.2f}x | Mediana: {stats_copia['50%']:>6.2f}x | Máx: {stats_copia['max']:>6.2f}x")
    else:
        print("Nenhuma matriz atendeu ao critério de convergência em ambas.")

    print("\n" + "=" * 70)
    print(f"DIFERENÇA RELATIVA DE ITERAÇÕES (%) [GPU vs CPU] (Filtrado: {total_conv} matrizes)")
    print("=" * 70)

    if not df_conv.empty:
        diff_stats = df_conv['iter_diff_pct'].describe()
        print(f"Média do Acréscimo   : {diff_stats['mean']:>8.2f}%")
        print(f"Mediana do Acréscimo : {diff_stats['50%']:>8.2f}%")
        print(f"Mínimo (GPU fez menos): {diff_stats['min']:>8.2f}%")
        print(f"Máximo Desvio        : {diff_stats['max']:>8.2f}%")
        print(f"Desvio Padrão        : {diff_stats['std']:>8.2f}%")
    print("=" * 70)

    # ==========================================
    # VISUALIZAÇÃO GRÁFICA UNIFICADA (3 PAINÉIS)
    # ==========================================
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(20, 5))

    cores_status = ["#2ecc71", "#3498db", "#e67e22", "#e74c3c"]

    # Painel 1: Status de Convergência Geral
    ax1 = sns.countplot(data=df, x='Status', order=categorias, palette=cores_status, ax=axes[0])
    axes[0].set_title('Status Geral de Convergência', fontsize=12, pad=10)
    axes[0].set_ylabel('Número de Matrizes')
    axes[0].set_xlabel('')
    axes[0].set_xticklabels(['Ambas', 'Só CPU', 'Só GPU', 'Nenhuma'], rotation=15)

    for p in ax1.patches:
        altura = int(p.get_height())
        if altura > 0:
            ax1.annotate(f'{altura}', 
                        (p.get_x() + p.get_width() / 2., p.get_height()), 
                        ha='center', va='center', 
                        xytext=(0, 6), textcoords='offset points', 
                        fontweight='bold', fontsize=9)

    # Painel 2: Speedup Puro vs com Cópia (Exclusivo Ambas Convergiram)
    if not df_conv.empty:
        sns.scatterplot(
            data=df_conv, x=speedup_col, y=speedup_copy_col,
            ax=axes[1], color='#2980b9', alpha=0.7, s=35
        )
        max_sp = max(df_conv[speedup_col].max(), df_conv[speedup_copy_col].max())
        axes[1].plot([0, max_sp], [0, max_sp], 'k--', linewidth=1.2, label='Ideal (Sem custo de cópia)')
        axes[1].set_title('Speedup (Apenas Ambas Convergiram)', fontsize=12, pad=10)
        axes[1].set_xlabel('Speedup Puro (Kernel)')
        axes[1].set_ylabel('Speedup com Cópia de Memória')
        axes[1].legend(fontsize=8, loc='upper left')

    # Painel 3: Diferença Relativa Percentual vs NNZ (Exclusivo Ambas Convergiram)
    if not df_conv.empty:
        sns.scatterplot(
            data=df_conv, x='NNZ', y='iter_diff_pct',
            ax=axes[2], color='#8e44ad', alpha=0.7, s=35
        )
        axes[2].axhline(0, color='black', linestyle='--', linewidth=1.2, label='0% (Mesmo nº de iterações)')
        axes[2].set_xscale('log')
        axes[2].set_yscale('symlog')
        axes[2].set_title('Dif. Relativa de Iterações % (Ambas Convergiram)', fontsize=12, pad=10)
        axes[2].set_xlabel('NNZ [escala log]')
        axes[2].set_ylabel('Diferença Relativa (%) [escala symlog]')
        axes[2].legend(fontsize=8)

    plt.tight_layout()
    nome_figura = 'analise_completa_cg.png'
    plt.savefig(nome_figura, dpi=300)
    print(f"\n[!] Painel gráfico atualizado exportado com sucesso como '{nome_figura}'.")

finally:
    # Restaura a saída padrão e fecha o arquivo de log
    sys.stdout.log.close()
    sys.stdout = sys.stdout.terminal
    print(f"[!] Relatório em texto gravado com sucesso em '{arquivo_txt}'.")
