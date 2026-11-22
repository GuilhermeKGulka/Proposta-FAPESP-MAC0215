# CHANGELOG GPU — diferenças entre o Gencan de referência e o Gencan CPU/GPU

> Objetivo: reduzir a dívida cognitiva. Quem conhece o Gencan original (`gencan-cg-standalone/`) deve conseguir entender, por este arquivo, tudo o que mudou em `gencan-cg-standalone-gpu/` e por quê.
> Começou na fase 6 (22/11/2026, branch `fase6`). Plano: `docs/FASE6.md`. Decisões gerais: `docs/PLANO.md` §5.

## Como usar

- **Referência:** `gencan-cg-standalone/` fica intocado. Toda entrada compara com ele.
- **Quando atualizar:** todo commit que mexe no código GPU, ou em arquivos compartilhados com a CPU, acrescenta uma entrada aqui.
- **Ordem:** as entradas ficam em ordem cronológica, a mais nova no fim. O §"Estado atual" resume o resultado líquido e é reescrito quando muda.

Formato de cada entrada:

```
### AAAA-MM-DD — título curto (commit abc1234)
- Arquivo/rotina: ...
- O que mudou (em relação a gencan-cg-standalone/): ...
- Por quê: ...
- Validação: teste(s) rodado(s) e resultado.
```

## Linha de base (antes da fase 6)

O código da fase 5 validada, no commit `ddf53fd` (`main`), é o ponto de partida. Nele, a GPU vive em arquivos paralelos:
- `gencan_gpu.f90`: cópia de `gencan.f90`, trocando `preparehv` + `cgm` por `hess_csr_update` + `cgm_gpu` e passando um `ctx`;
- `cutestma_gpu.f90`: cópia de `cutestma.f90`;
- `hess_csr.f90`: monta o CSR da face na GPU;
- `cgm_gpu.f90`, `cgm_kernels.cpp`, `hipsparse_generic.f90`.

`lssit.f90`, `cgm.f90`, `gencan.f90` e `cutestma.f90` dentro de `gencan-cg-standalone-gpu/sources/` são idênticos aos da referência.

Resultados de referência (os da fase 6 serão comparados com eles):
- **CUTEst:** `pequena.txt` (7 problemas) com `precond` igual nos dois lados, em `gencan-cg-standalone-gpu/resultados/table_{cpu,gpu}_fase5_{nopre,pre}.txt`.
  - ROSENBR, BDQRTIC, CRAGGLVY, EDENSCH, TORSION1 e OBSTCLAE têm `iter`/`fcnt`/`gcnt`/`hcnt` idênticos entre CPU e GPU.
  - BQPGAUSS diverge por arredondamento (97/90 sem precond, 89/118 com), com `f` igual a 1e-12. Ver PLANO §6 e §8.
- **Testes sintéticos:**
  - `test_hess_csr`: 5/5, erro ≤ 1.6e-16;
  - `test_cgm_gpu`: 9/9, erro relativo em `d` ≤ 1.5e-15.

## Estado atual

Igual à linha de base: nenhuma mudança de código na fase 6 até agora.

## Entradas

### 2026-11-22 — início da fase 6
- Arquivo/rotina: `docs/CHANGELOG_GPU.md` (novo), branch `fase6`.
- O que mudou: nenhum código. Este arquivo foi criado e a linha de base registrada.
- Por quê: passo 1 do `docs/FASE6.md`.
- Validação: não se aplica.
