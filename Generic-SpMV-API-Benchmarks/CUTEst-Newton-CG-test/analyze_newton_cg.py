#!/usr/bin/env python3
"""
Analise grafica do benchmark Newton-CG focada na validacao da proposta FAPESP.
Uso: python analyze_newton_cg.py [prefixo] [--csv caminho_alternativo]

Exemplo: python analyze_newton_cg.py amd
Lera 'amd_results_newton_cg.csv' e salvará tudo em 'plots_amd/'.
"""
import warnings
import glob
import os
import sys
import argparse
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

warnings.filterwarnings('ignore')
plt.rcParams.update({'figure.dpi': 150, 'axes.grid': True,
                     'grid.alpha': 0.3, 'font.size': 11})

# ---------------------------------------------------------------------------
# Argumentos
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Gera graficos do benchmark Newton-CG.")
parser.add_argument("prefix", nargs="?", default="amd",
                    help="Prefixo do hardware (ex: amd, nvidia). Default: amd")
parser.add_argument("--csv", default=None,
                    help="Caminho do CSV. Se omitido, usa '<prefix>_results_newton_cg.csv'")
args = parser.parse_args()

# Configuração de caminhos
outdir = f"plots_{args.prefix}"
csv_path = args.csv if args.csv else f"{args.prefix}_results_newton_cg.csv"

if not os.path.exists(csv_path):
    print(f"Erro: Arquivo '{csv_path}' nao encontrado.")
    print(f"Certifique-se de que o arquivo esteja na pasta ou passe --csv <caminho>")
    sys.exit(1)

os.makedirs(outdir, exist_ok=True)

# ---------------------------------------------------------------------------
# Classe Tee para salvar o terminal num arquivo de texto
# ---------------------------------------------------------------------------
class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

log_file = open(os.path.join(outdir, f"{args.prefix}_resumo.txt"), "w")
original_stdout = sys.stdout
sys.stdout = Tee(sys.stdout, log_file)

# ---------------------------------------------------------------------------
# Carga e Tratamento (Robusto contra CSV mal formatado)
# ---------------------------------------------------------------------------
with open(csv_path, encoding="utf-8-sig") as f:
    lines = [ln.rstrip().rstrip(",").rstrip() for ln in f if ln.strip()]

raw = pd.read_csv(io.StringIO("\n".join(lines)), header=None, dtype=str, skipinitialspace=True)
hdr = [str(c).strip() for c in raw.iloc[0]]
data = raw.iloc[1:].reset_index(drop=True)

# Completa o header se houver menos nomes que colunas
while len(hdr) < data.shape[1]:
    hdr.append("Max_err_rel" if len(hdr) == data.shape[1] - 1 else f"Unnamed_{len(hdr)}")

data.columns = hdr[:data.shape[1]]

# Converte para numérico
for c in data.columns:
    if c != 'Problem':
        data[c] = pd.to_numeric(data[c], errors='coerce')

# Filtra apenas onde ambos convergiram
conv = data[(data['istop_CPU'] == 0) & (data['istop_GPU'] == 0)].copy()

def save_plot(name):
    filename = f"{args.prefix}_{name}.png"
    path = os.path.join(outdir, filename)
    plt.savefig(path, bbox_inches='tight')
    plt.close()

# ---------------------------------------------------------------------------
# Resumo e Validacao (Texto)
# ---------------------------------------------------------------------------
print("=" * 70)
print(f" RESULTADOS DA VALIDACAO — {args.prefix.upper()}")
print("=" * 70)
print(f"Total de problemas testados: {len(data)}")
print(f"Convergiram nos dois:        {len(conv)}")
print(f"Divergiram/Falharam:         {len(data) - len(conv)}")

# Criterio 2: Trajetoria Identica
if len(conv) > 0:
    ni_diff = (conv['NI_CPU'] - conv['NI_GPU']).abs()
    cg_diff = (conv['last_cg_iter_CPU'] - conv['last_cg_iter_GPU']).abs()
    print("\n--- Criterio 2: Caminho de convergencia identico ---")
    print(f"Problemas com iteracoes Newton diferentes: {(ni_diff > 0).sum()} de {len(conv)}")
    print(f"Problemas com iteracoes CG (ultimo passo) diferentes: {(cg_diff > 0).sum()} de {len(conv)}")

    # Criterio 1: Precisao Numerica
    valid_sol = conv['Sol_Diff'].dropna()
    print("\n--- Criterio 1: Precisao da solucao ---")
    print(f"Max |x_cpu - x_gpu|: {valid_sol.max():.2e}")
    print(f"Mediana |x_cpu - x_gpu|: {valid_sol.median():.2e}")

    # Criterio 3: Speedup para N >= 10^5
    large = conv[conv['N'] >= 100000]
    print("\n--- Criterio 3: Speedup para N >= 10^5 ---")
    if len(large) > 0:
        print(f"Problemas de grande porte: {len(large)}")
        print(f"Speedup medio (Total): {large['Speedup'].mean():.2f}x")
        print(f"Speedup medio (CG):    {large['Speedup_CG'].mean():.2f}x")
    else:
        print("Nenhum problema com N >= 10^5 testado.")


# ---------------------------------------------------------------------------
# Plot 1: Validacao do Newton (Iteracoes Medias CG)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# NI_CPU vs NI_GPU
axes[0].scatter(conv['NI_CPU'], conv['NI_GPU'], s=40, alpha=0.7, color='C0', edgecolor='k')
axes[0].plot([0, 120], [0, 120], 'k--', label='y = x')
axes[0].set(xlabel='Iteracoes Newton (CPU)', ylabel='Iteracoes Newton (GPU)',
            title='Iteracoes Newton: CPU vs GPU')
axes[0].legend()
axes[0].set_aspect('equal')

# MÉDIA das iteracoes CG por passo de Newton
if 'avg_cg_iter_CPU' in conv.columns and 'avg_cg_iter_GPU' in conv.columns:
    mask = (conv['avg_cg_iter_CPU'] > 0) & (conv['avg_cg_iter_GPU'] > 0)
    axes[1].scatter(conv.loc[mask, 'avg_cg_iter_CPU'], conv.loc[mask, 'avg_cg_iter_GPU'],
                    s=40, alpha=0.7, color='C1', edgecolor='k')
    axes[1].plot([1, 1e5], [1, 1e5], 'k--', label='y = x')
    axes[1].set(xscale='log', yscale='log', xlabel='Media de iteracoes CG por Newton (CPU)',
                ylabel='Media de iteracoes CG por Newton (GPU)',
                title='Media de Iteracoes CG: CPU vs GPU')
    axes[1].legend()
else:
    axes[1].text(0.5, 0.5, 'Colunas avg_cg_iter ausentes no CSV', ha='center', va='center')

plt.suptitle('Validacao: A GPU nao alterou a trajetoria de convergencia', fontsize=14, y=1.02)
plt.tight_layout()
save_plot('plot1_validacao_convergencia')


# ---------------------------------------------------------------------------
# Plot 2: Escalabilidade (Speedup vs N)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(conv['N'], conv['Speedup'], s=40, alpha=0.7, label='Speedup Total', color='C0', edgecolor='k')
ax.scatter(conv['N'], conv['Speedup_CG'], s=40, alpha=0.5, label='Speedup CG', color='C1', marker='^', edgecolor='k')

ax.axhline(1.0, color='k', linestyle='--', linewidth=1)
ax.axvline(1e5, color='r', linestyle='--', linewidth=1, label='Meta Proposta (N=10^5)')

ax.set(xscale='log', yscale='log',
       xlabel='Numero de variaveis (N) [log]',
       ylabel='Speedup CPU/GPU [log]',
       title='Escalabilidade: Ganho da GPU em funcao do porte do problema')
ax.legend(loc='upper left')
plt.tight_layout()
save_plot('plot2_escalabilidade')


# ---------------------------------------------------------------------------
# Plot 3: Lei de Amdahl (O Gargalo Serial)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
frac_cg = conv['t_cg_CPU_ms'] / conv['t_total_CPU_ms']
sc = ax.scatter(conv['Speedup_CG'], conv['Speedup'],
                c=frac_cg, cmap='viridis',
                s=50, alpha=0.8, edgecolor='k', norm=LogNorm())

lims = [1e-3, 1e3]
ax.plot(lims, lims, 'k--', linewidth=1, label='Limite de Amdahl (y=x)')
ax.set(xscale='log', yscale='log', xlim=lims, ylim=[1e-3, 1e3],
       xlabel='Speedup do CG apenas',
       ylabel='Speedup Total da aplicacao',
       title='Lei de Amdahl: O custo serial (CUTEst) limita o ganho total')
ax.axhline(1, color='gray', linestyle=':', linewidth=0.8)
ax.axvline(1, color='gray', linestyle=':', linewidth=0.8)

cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Fracao do tempo da CPU gasta no CG')
ax.legend(loc='upper left')
plt.tight_layout()
save_plot('plot3_amdahl')


# ---------------------------------------------------------------------------
# Plot 4: Precisao Numerica (Norma do Gradiente)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

gn_c, gn_g = conv['Gnorm_CPU'], conv['Gnorm_GPU']
v = (gn_c > 0) & (gn_g > 0) & gn_c.notna() & gn_g.notna()
axes[0].scatter(gn_c[v], gn_g[v], s=40, alpha=0.7, color='C2', edgecolor='k')
axes[0].plot([1e-15, 1e5], [1e-15, 1e5], 'k--', label='y = x')
axes[0].axhline(1e-6, color='r', linestyle='--', linewidth=0.8, label='Tolerancia (1e-6)')
axes[0].axvline(1e-6, color='r', linestyle='--', linewidth=0.8)
axes[0].set(xscale='log', yscale='log', xlim=[1e-15, 1e5], ylim=[1e-15, 1e5],
            xlabel='||grad|| final — CPU', ylabel='||grad|| final — GPU',
            title='Norma do gradiente final')
axes[0].legend()

valid = conv['Sol_Diff'].dropna()
valid = valid[(valid > 0) & (valid < 1e30)]
axes[1].hist(np.log10(valid), bins=30, color='C3', alpha=0.7, edgecolor='k')
axes[1].axvline(-6, color='r', linestyle='--', linewidth=1, label='Tolerancia (1e-6)')
axes[1].set(xlabel='log10 |x_cpu - x_gpu|', ylabel='Numero de problemas',
            title='Diferenca na solucao final')
axes[1].legend()

plt.suptitle('Validacao: Alta precisao numerica preservada', fontsize=14, y=1.02)
plt.tight_layout()
save_plot('plot4_precisao')


# ---------------------------------------------------------------------------
# Plot 5: Tempo de Transferencia (Overhead) vs Tamanho
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))
if 't_GPU_prep_ms' in conv.columns:
    overhead_ratio = (conv['t_GPU_prep_ms'] / conv['t_cg_GPU_ms']) * 100
    v = overhead_ratio > 0
    ax.scatter(conv.loc[v, 'NNZ'], overhead_ratio[v], s=40, alpha=0.7, color='C4', edgecolor='k')
    ax.set(xscale='log', yscale='log',
           xlabel='NNZ da Hessiana [log]',
           ylabel='Overhead do SpMV (% do tempo do CG)',
           title='Overhead de inicializacao do SpMV (hipSPARSE)')
plt.tight_layout()
save_plot('plot5_overhead')


# ---------------------------------------------------------------------------
# Top e Bottom Tabelas no terminal
# ---------------------------------------------------------------------------
cols = ['Problem', 'N', 'NNZ', 'NI_CPU', 'NI_GPU', 'Speedup', 'Speedup_CG']
print("\n" + "=" * 70 + "\nTOP 10 SPEEDUP TOTAL\n" + "=" * 70)
print(conv.nlargest(10, 'Speedup')[cols].to_string(index=False))

print("\n" + "=" * 70 + "\nTOP 10 SPEEDUP CG (Apenas o Solver)\n" + "=" * 70)
print(conv.nlargest(10, 'Speedup_CG')[cols].to_string(index=False))

print(f"\nPlots salvos em '{outdir}':")
for f in sorted(glob.glob(os.path.join(outdir, '*.png'))):
    print(f"  {f}")

# ---------------------------------------------------------------------------
# Restaura a saída padrão e fecha o arquivo de log
# ---------------------------------------------------------------------------
sys.stdout = original_stdout
log_file.close()
print(f"Log do terminal salvo em: {os.path.join(outdir, f'{args.prefix}_resumo.txt')}")