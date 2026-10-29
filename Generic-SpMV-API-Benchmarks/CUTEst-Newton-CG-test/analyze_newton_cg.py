#!/usr/bin/env python3
"""
Analise grafica do benchmark Newton-CG CPU vs GPU.
Le results_newton_cg.csv e gera PNGs + resumo de texto.
"""
import warnings
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

warnings.filterwarnings('ignore')
plt.rcParams.update({'figure.dpi': 150, 'axes.grid': True,
                     'grid.alpha': 0.3, 'font.size': 10})


# ---------------------------------------------------------------------------
# Carga
# ---------------------------------------------------------------------------
df = pd.read_csv('results_newton_cg.csv')
for c in df.columns:
    if c != 'Problem':
        df[c] = pd.to_numeric(df[c], errors='coerce')

for c in ['t_total_CPU_ms', 't_total_GPU_ms', 't_cg_CPU_ms', 't_cg_GPU_ms', 'NNZ']:
    df.loc[df[c] <= 0, c] = np.nan

df['cg_frac_CPU'] = (df['t_cg_CPU_ms'] / df['t_total_CPU_ms']).clip(0, 1)
df['cg_frac_GPU'] = (df['t_cg_GPU_ms'] / df['t_total_GPU_ms']).clip(0, 1)
df['host_CPU'] = df['t_total_CPU_ms'] - df['t_cg_CPU_ms']
df['host_GPU'] = df['t_total_GPU_ms'] - df['t_cg_GPU_ms'] - df['t_GPU_prep_ms'].fillna(0)

conv = df[(df['istop_CPU'] == 0) & (df['istop_GPU'] == 0)].copy()


# ---------------------------------------------------------------------------
# Resumo
# ---------------------------------------------------------------------------
print("=" * 70)
print(f"Total de problemas no CSV: {len(df)}")
print(f"Convergiram nos dois:      {len(conv)}")
print(f"Divergiram:                {len(df) - len(conv)}")
print("=" * 70)
for label, col in [("Speedup total", "Speedup"), ("Speedup CG", "Speedup_CG")]:
    s = conv[col].dropna()
    print(f"{label}: mediana={s.median():.3f}  media={s.mean():.3f}  "
          f"min={s.min():.3f}  max={s.max():.2f}")
print()
for lim in [0.5, 1.0, 2.0, 5.0, 10.0]:
    n = (conv['Speedup'] > lim).sum()
    print(f"  Speedup > {lim:>4}: {n:>3} / {len(conv)}")
print()
cross = conv[conv['Speedup'] > 1]
if len(cross):
    print(f"Crossover (Speedup>1): menor NNZ={cross['NNZ'].min():.0f}  "
          f"mediana NNZ={cross['NNZ'].median():.0f}")


# ---------------------------------------------------------------------------
# Plot 1 — Speedup × NNZ
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 6))
ax.scatter(conv['NNZ'], conv['Speedup'], s=25, alpha=0.7, label='Speedup total',
           color='C0', edgecolor='k', linewidth=0.3)
ax.scatter(conv['NNZ'], conv['Speedup_CG'], s=25, alpha=0.5, label='Speedup CG',
           color='C1', marker='^', edgecolor='k', linewidth=0.3)
ax.axhline(1.0, color='k', linestyle='--', linewidth=1, label='paridade')
ax.set(xscale='log', yscale='log',
       xlabel='NNZ da Hessiana (log)',
       ylabel='Speedup CPU/GPU (log)',
       title='Speedup em funcao do tamanho: ganho cresce com a Hessiana')
ax.legend(loc='upper left', framealpha=0.95)
plt.tight_layout(); plt.savefig('plot1_speedup_vs_nnz.png'); plt.close()


# ---------------------------------------------------------------------------
# Plot 2 — Amdahl
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(conv['Speedup_CG'], conv['Speedup'],
                c=conv['cg_frac_CPU'], cmap='viridis',
                s=30, alpha=0.75, edgecolor='k', linewidth=0.3, norm=LogNorm())
lims = [1e-5, 1e3]
ax.plot(lims, lims, 'k--', linewidth=0.8, label='y = x (limite Amdahl)')
ax.set(xscale='log', yscale='log', xlim=lims, ylim=[1e-2, 1e2],
       xlabel='Speedup CG apenas',
       ylabel='Speedup total da aplicacao',
       title='Amdahl: CG mais rapido nao implica app mais rapida')
ax.axhline(1, color='gray', linestyle=':', linewidth=0.6)
ax.axvline(1, color='gray', linestyle=':', linewidth=0.6)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Fracao do tempo da CPU em CG')
ax.legend(loc='upper left')
plt.tight_layout(); plt.savefig('plot2_amdahl.png'); plt.close()


# ---------------------------------------------------------------------------
# Plot 3 — Crossover do CG
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))
ax.scatter(conv['NNZ'], conv['cg_frac_CPU'] * 100, s=25, alpha=0.7,
           color='C0', edgecolor='k', linewidth=0.3)
ax.axhline(50, color='r', linestyle='--', linewidth=0.8, label='50% do tempo')
ax.set(xscale='log', yscale='log',
       xlabel='NNZ (log)', ylabel='% tempo da CPU em CG (log)',
       title='A partir de qual NNZ o CG domina o tempo da CPU')
ax.legend()
plt.tight_layout(); plt.savefig('plot3_cg_fraction.png'); plt.close()


# ---------------------------------------------------------------------------
# Plot 4 — Trabalho serializado igual nos dois backends
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7))
ax.scatter(conv['host_CPU'], conv['host_GPU'], s=25, alpha=0.7,
           color='C2', edgecolor='k', linewidth=0.3)
lims = [min(conv[['host_CPU', 'host_GPU']].min()),
        conv[['host_CPU', 'host_GPU']].max().max()]
ax.plot(lims, lims, 'k--', linewidth=0.8, label='y = x')
ax.set(xscale='log', yscale='log', xlim=lims, ylim=lims,
       xlabel='Tempo nao-CG na execucao CPU (ms)',
       ylabel='Tempo nao-CG na execucao GPU (ms)',
       title='Trabalho serializado (CUTEst) e igual nos dois backends')
ax.legend()
plt.tight_layout(); plt.savefig('plot4_host_time.png'); plt.close()


# ---------------------------------------------------------------------------
# Plot 5 — Distribuicao (hist + CDF)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
bins = np.logspace(-3, 2, 35)
axes[0].hist(conv['Speedup'], bins=bins, alpha=0.65, label='Speedup total',
             color='C0', edgecolor='k')
axes[0].hist(conv['Speedup_CG'], bins=bins, alpha=0.55, label='Speedup CG',
             color='C1', edgecolor='k')
axes[0].axvline(1, color='k', linestyle='--', linewidth=1)
axes[0].set(xscale='log', xlabel='Speedup (log)',
            ylabel='Numero de problemas', title='Histograma de speedup')
axes[0].legend()

s1 = np.sort(conv['Speedup'].dropna())
s2 = np.sort(conv['Speedup_CG'].dropna())
axes[1].plot(s1, np.arange(1, len(s1) + 1) / len(s1), label='Speedup total',
             color='C0', linewidth=1.5)
axes[1].plot(s2, np.arange(1, len(s2) + 1) / len(s2), label='Speedup CG',
             color='C1', linewidth=1.5)
axes[1].axvline(1, color='k', linestyle='--', linewidth=1)
axes[1].set(xscale='log', xlabel='Speedup (log)',
            ylabel='Fracao acumulada', title='CDF de speedup')
axes[1].legend(loc='lower right')
plt.tight_layout(); plt.savefig('plot5_distribution.png'); plt.close()


# ---------------------------------------------------------------------------
# Plot 6 — Corretude
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
valid = conv['Sol_Diff'].dropna()
valid = valid[(valid > 0) & (valid < 1e30)]
axes[0].hist(np.log10(valid), bins=30, color='C0', alpha=0.7, edgecolor='k')
axes[0].axvline(-8, color='r', linestyle='--', linewidth=1, label='tol 1e-8')
axes[0].set(xlabel='log10 |x_cpu - x_gpu|', ylabel='Numero de problemas',
            title='Diferenca na solucao final (convergiram)')
axes[0].legend()

gn_c, gn_g = conv['Gnorm_CPU'], conv['Gnorm_GPU']
v = (gn_c > 0) & (gn_g > 0) & gn_c.notna() & gn_g.notna()
axes[1].scatter(gn_c[v], gn_g[v], s=25, alpha=0.7,
                color='C2', edgecolor='k', linewidth=0.3)
lims = [1e-15, 1e5]
axes[1].plot(lims, lims, 'k--', linewidth=0.8, label='y = x')
axes[1].axhline(1e-8, color='r', linestyle='--', linewidth=0.6)
axes[1].axvline(1e-8, color='r', linestyle='--', linewidth=0.6)
axes[1].set(xscale='log', yscale='log', xlim=lims, ylim=lims,
            xlabel='||grad|| final — CPU', ylabel='||grad|| final — GPU',
            title='Norma do gradiente final')
axes[1].legend()
plt.tight_layout(); plt.savefig('plot6_correctness.png'); plt.close()


# ---------------------------------------------------------------------------
# Plot 7 — Iteracoes CG
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
axes[0].scatter(conv['last_cg_iter_CPU'], conv['last_cg_iter_GPU'],
                s=25, alpha=0.7, color='C0', edgecolor='k', linewidth=0.3)
lims = [1, 1e5]
axes[0].plot(lims, lims, 'k--', linewidth=0.8)
axes[0].set(xscale='log', yscale='log', xlim=lims, ylim=lims,
            xlabel='Ultimo CG iter — CPU', ylabel='Ultimo CG iter — GPU',
            title='Convergencia equivalente entre backends')

sc = axes[1].scatter(conv['last_cg_iter_CPU'], conv['Speedup_CG'],
                     c=conv['NNZ'], cmap='viridis', s=25, alpha=0.75,
                     edgecolor='k', linewidth=0.3, norm=LogNorm())
axes[1].axhline(1, color='k', linestyle='--', linewidth=0.8)
axes[1].set(xscale='log', yscale='log',
            xlabel='Iteracoes CG do ultimo passo Newton',
            ylabel='Speedup CG',
            title='Mais iteracoes = melhor amortizacao do overhead')
cbar = plt.colorbar(sc, ax=axes[1]); cbar.set_label('NNZ')
plt.tight_layout(); plt.savefig('plot7_cg_behavior.png'); plt.close()


# ---------------------------------------------------------------------------
# Top e bottom
# ---------------------------------------------------------------------------
cols = ['Problem', 'N', 'NNZ', 'Speedup', 'Speedup_CG',
        'cg_frac_CPU', 't_total_CPU_ms', 't_total_GPU_ms']

print("\n" + "=" * 70 + "\nTOP 10 SPEEDUP\n" + "=" * 70)
top = conv.nlargest(10, 'Speedup')[cols].copy()
top['cg_frac_CPU'] = top['cg_frac_CPU'].map(lambda x: f"{x*100:5.1f}%")
print(top.to_string(index=False))

print("\n" + "=" * 70 + "\nBOTTOM 10 SPEEDUP (GPU nao compensa)\n" + "=" * 70)
bot = conv.nsmallest(10, 'Speedup')[cols].copy()
bot['cg_frac_CPU'] = bot['cg_frac_CPU'].map(lambda x: f"{x*100:5.1f}%")
print(bot.to_string(index=False))


# ---------------------------------------------------------------------------
# Verificacao da lei de Amdahl
# ---------------------------------------------------------------------------
print("\n" + "=" * 70 + "\nVALIDACAO DA LEI DE AMDAHL\n" + "=" * 70)
pred = 1 + conv['cg_frac_CPU'] * (conv['Speedup_CG'] - 1)
mask = pred.notna() & conv['Speedup'].notna() & (pred > 0) & (conv['Speedup'] > 0)
err = np.abs(np.log(conv['Speedup'][mask]) - np.log(pred[mask]))
print("Modelo: Speedup_total = 1 + f_CPU * (Speedup_CG - 1)")
print(f"  Erro log medio:   {err.mean():.3f}")
print(f"  Erro log maximo:  {err.max():.3f}")
print(f"  Acuracia (erro<20%): {(err < np.log(1.2)).mean()*100:.1f}%")

print("\nPlots gerados:")
for f in sorted(glob.glob('plot*.png')):
    print(f"  {f}")