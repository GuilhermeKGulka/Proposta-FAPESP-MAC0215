import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

if len(sys.argv) < 2:
    print("Uso: python plot_newton.py results_newton.csv")
    sys.exit(1)

df = pd.read_csv(sys.argv[1])

# Filtragem e normalização de nomes de colunas
col_n = [c for c in df.columns if c.lower() == 'n'][0]
col_nnz = [c for c in df.columns if c.lower() == 'nnz'][0]
col_time = [c for c in df.columns if 'time' in c.lower()][0]
col_cg_iters = [c for c in df.columns if 'cg_iters' in c.lower()][0]

df['density'] = df[col_nnz] / (df[col_n] * df[col_n])

# Mapeamento do tamanho (Área dos pontos baseada nas iterações totais do CG)
r_min, r_max = 3.0, 15.0
log_cg = np.log10(df[col_cg_iters].replace(0, 1)) # evita log(0)
cg_norm = (log_cg - log_cg.min()) / (log_cg.max() - log_cg.min() + 1e-9)

df['radius'] = r_min + cg_norm * (r_max - r_min)
df['area_pts'] = (df['radius']) ** 2

df = df.sort_values(by=col_cg_iters, ascending=False).reset_index(drop=True)

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 11})
fig, ax = plt.subplots(figsize=(10, 6))

scatter = ax.scatter(
    df[col_n],
    df['density'],
    c=df[col_time],
    cmap='viridis',
    s=df['area_pts'],
    alpha=0.8,
    edgecolor='black',
    linewidth=0.5
)

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Tempo Total (ms)')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Dimensão (N) [escala log]')
ax.set_ylabel(r'Densidade da Hessiana ($NNZ / N^2$) [escala log]')
ax.set_title('Performance do Newton-CG em GPU (Problemas CUTEst)')

# Legenda baseada nas iterações de CG
cg_ticks = [10, 100, 1000, 10000]
size_handles = []
for val in cg_ticks:
    log_val = np.log10(val)
    norm_val = (log_val - log_cg.min()) / (log_cg.max() - log_cg.min() + 1e-9)
    r_val = r_min + norm_val * (r_max - r_min)
    size_handles.append(
        Line2D([0], [0], marker='o', color='w', label=str(val),
               markerfacecolor='gray', markersize=r_val, markeredgecolor='black', alpha=0.8)
    )

custom_legend = [Line2D([], [], color='none', label=r'$\bf{Total\ CG\ Iters}$')] + size_handles
ax.legend(handles=custom_legend, loc='upper right', bbox_to_anchor=(1.25, 1))

plt.tight_layout()
plt.savefig('plot_newton_scaled.png', dpi=300)
plt.show()