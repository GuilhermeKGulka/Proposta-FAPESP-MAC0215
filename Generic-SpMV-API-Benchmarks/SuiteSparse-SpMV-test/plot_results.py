import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

# 1. Carregamento do dataset
df = pd.read_csv('results_cg_1_with_kinds.csv')

# Identificação dinâmica das colunas
nnz_col = [c for c in df.columns if c.lower() == 'nnz'][0]
density_col = [c for c in df.columns if c.lower() == 'density'][0]
speedup_col = [c for c in df.columns if 'speedup' in c.lower() and 'copy' in c.lower()][0]
kind_col = [c for c in df.columns if c.lower() == 'kind'][0]

# 3. Mapeamento dos tamanhos focado no intervalo [1, 100]
# Raio ajustado para equilibrar visualmente o gráfico e a legenda
r_min, r_max = 2.0, 10.0 

log_speedup = np.log10(df[speedup_col])
speedup_norm = log_speedup + 1.0

df['radius'] = r_min + speedup_norm * (r_max - r_min)
df['area_pts'] = (df['radius']) ** 2  # Área em pontos do scatter

# 4. Ordenação para que as bolinhas menores fiquem sobrepostas no topo
df = df.sort_values(by=speedup_col, ascending=False).reset_index(drop=True)

# 5. Renderização do gráfico
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 11})
fig, ax = plt.subplots(figsize=(11, 7))

kinds = sorted(df[kind_col].dropna().unique())
palette = sns.color_palette("tab20", n_colors=len(kinds))
color_map = dict(zip(kinds, palette))

for _, row in df.iterrows():
    ax.scatter(
        row[nnz_col],
        row[density_col],
        s=row['area_pts'],
        color=color_map[row[kind_col]],
        alpha=0.7,
        edgecolor='black',
        linewidth=0.5
    )

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Elementos Não-Nulos (NNZ) [escala log]')
ax.set_ylabel(r'Densidade ($NNZ / N^2$) [escala log]')
ax.set_title('Performance vs Densidade e Estrutura da Matriz')

# 6. Legendas customizadas
kind_handles = [
    Line2D([0], [0], marker='o', color='w', label=str(k),
           markerfacecolor=color_map[k], markersize=7.5, markeredgecolor='black')
    for k in kinds
]

# Ticks destacados para a legenda do speedup
speedup_ticks = [0.1, 1, 10, 100]
size_handles = []
for val in speedup_ticks:
    log_val = np.log10(val)
    norm_val = log_val + 1.0
    r_val = r_min + norm_val * (r_max - r_min)
    label_txt = val
    size_handles.append(
        Line2D([0], [0], marker='o', color='w', label=label_txt,
               markerfacecolor='gray', markersize=r_val, markeredgecolor='black', alpha=1)
    )

custom_legend = []
custom_legend += (
    [Line2D([], [], color='none', label=r'$\bf{Tipo\ de\ Matriz}$')] + kind_handles +
    [Line2D([], [], color='none', label='')] +
    [Line2D([], [], color='none', label=r'$\bf{Speedup\ (1\ a\ 100)}$')] + size_handles
)

ax.legend(
    handles=custom_legend,
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    borderaxespad=0,
    fontsize=8.5,
    labelspacing=0.65,
    handletextpad=0.8,
    borderpad=0.5
)

plt.tight_layout()
plt.savefig('plot_results_scaled.png', dpi=300)
plt.show()