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

- **Fontes Fortran:** iguais à linha de base. Nenhum `.f90` mudou até agora na fase 6.
- **Build:**
  - a GPU é opcional (`GENCAN_GPU` AUTO/ON/OFF), assim como o CUTEst;
  - preset `cpu` para compilar sem HIP;
  - `toyprob_cpu` voltou ao build.

## Entradas

### 2026-11-22 — início da fase 6
- Arquivo/rotina: `docs/CHANGELOG_GPU.md` (novo), branch `fase6`.
- O que mudou: nenhum código. Este arquivo foi criado e a linha de base registrada.
- Por quê: passo 1 do `docs/FASE6.md`.
- Validação: não se aplica.

### 2026-11-22 — CMake: GPU e CUTEst opcionais (fase 6, passo 2)
- **Arquivo/rotina:** `gencan-cg-standalone-gpu/CMakeLists.txt`, `CMakePresets.json`, comentário de `run-tests.sh`.
- **O que mudou (em relação a `gencan-cg-standalone/`):**
  - A referência não tem CMake: compila com `gfortran -O3` pelos scripts `compile.sh`/`run-cutest.sh`.
  - Aqui a opção `GENCAN_GPU` (`AUTO`/`ON`/`OFF`) decide se a parte GPU entra no build. `AUTO` procura o backend de `HIP_PLATFORM` e, se faltar algo, compila só a CPU com aviso.
  - A busca do backend é a mesma dos `Generic-SpMV-API-Benchmarks` (ROCm e CUDA), agora não fatal.
  - CUTEst opcional: sem `env_cutest.sh`, não há `cutestma_*`.
  - Novo alvo `toyprob_cpu`.
  - `lssit/cgm/gencan` são compilados com `-cpp` (`Fortran_PREPROCESS`) e grava seus módulos em `mod/cpu`. A parte GPU grava os seus em `mod/gpu`.
  - Presets: `amd`/`nvidia` com `GENCAN_GPU=ON` e o novo `cpu` com `OFF`.
- **Por quê:** pedido do orientador de que o mesmo pacote compile sem HIP. É também a base para compilar os mesmos fontes com e sem `-DGENCAN_GPU` (D2(b) do `docs/FASE6.md`). O CUTEst opcional permite compilar e rodar o `toyprob` sem instalar nada além do gfortran.
- **Validação (22/11/2026, RX 9070 XT):**
  - Preset `cpu` sem o ambiente do CUTEst: compila só `gencan_cpu` + `toyprob_cpu`, sem nada de HIP no `ldd`.
  - Bit a bit: a saída de `toyprob_cpu` e de `cutestma_cpu` (build `cpu`) é idêntica à da referência compilada pelos scripts dela, nos 7 problemas da `pequena.txt`. Só a linha de tempo foi ignorada.
  - `HIP_PLATFORM=nvidia` sem CUDA: em AUTO compila só a CPU e avisa "falta CUDA Toolkit"; em ON aborta.
  - amd AUTO com o pacote `hip` desabilitado (simulando ROCm incompleto): compila só a CPU, com aviso.
  - Preset `amd`: todos os alvos compilam; `run-tests.sh` dá 5/5 e 9/9; ROSENBR tem o mesmo `f` em CPU e GPU.
  - **[NÃO VERIFICADO]** build NVIDIA com CUDA (precisa da máquina com RTX 5070 Ti).
