# Gencan CPU/GPU

Gencan (minimização com restrições de caixa, `l ≤ x ≤ u`) com o gradiente conjugado do passo de Newton rodando na GPU, via HIP. Funciona em AMD (ROCm) e NVIDIA (CUDA).

Há um só código. Sem a macro `GENCAN_GPU`, os fontes se comportam exatamente como o Gencan de referência (`../gencan-cg-standalone/`). Com ela, a GPU entra em três pontos, por analogia com o original. Por enquanto, só o caso de Hessiana explícita (`truehp = .true.`) usa a GPU.

## Como a GPU entra

O fluxo do Gencan não muda: `evalh` → `comphpre` → `preparehv` → `cgm`.

| Onde | O que faz no modo GPU |
|---|---|
| `lssit.f90` (`#ifdef GENCAN_GPU`) | O `lssit_type` ganha `usegpu` (padrão `.true.`) e o estado da GPU. O `preparehv` manda a Hessiana da face para a GPU, no lugar de compactá-la na CPU. |
| `cgm.f90` (`#ifdef GENCAN_GPU`) | Se a Hessiana da face está na GPU, chama `cgm_gpu` (mesma assinatura) e retorna. |
| `lssit_gpu.f90` | Hessiana da face na GPU: CSR montado a partir do COO do `evalh` só quando a estrutura ou a face mudam; a cada passo, envio dos valores + scatter; SpMV. Também guarda os vetores de trabalho do CG. Toda a memória da GPU é alocada uma vez. |
| `cgm_gpu.f90` | CG com a mesma ordem de testes do `cgm`, vetores no device e uma sincronização por iteração. |

O `gencan.f90`, o `cutestma.f90` e o `toyprob.f90` são iguais aos da referência.

Detalhes:
- **diferenças em relação à referência, mudança por mudança:** `../docs/CHANGELOG_GPU.md`;
- **conversão para CSR (corretude e eficiência):** `../docs/CONVERSAO_CSR.md`;
- **plano e decisões:** `../docs/FASE6.md` e `../docs/PLANO.md`.

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `sources/gencan.f90`, `lssit.f90`, `cgm.f90` | Gencan. A GPU entra só sob `#ifdef GENCAN_GPU`. |
| `sources/lssit_gpu.f90`, `cgm_gpu.f90` | Parte GPU (Hessiana da face e CG). |
| `sources/cgm_kernels.cpp` | Kernels HIP do CG (reduções, atualizações, scatter). |
| `sources/hipsparse_generic.f90` | Bindings `bind(C)` para a API genérica do hipSPARSE (SpMV). |
| `sources/cutestma.f90`, `toyprob.f90` | Drivers: CUTEst e o problema de exemplo. |
| `sources/hess_csr.f90` | Versão da fase 5 da Hessiana na GPU, mantida só para comparação de tempos no `test_hess_csr`. |
| `tests/test_hess_csr.f90` | `preparehv` com GPU contra `preparehv` + `computehv` da CPU, mais a medição de tempos. |
| `tests/test_cgm_gpu.f90` | O mesmo `cgm` com e sem GPU, em 9 cenários. |
| `resultados/` | Tabelas CUTEst de referência de cada fase. |

## Compilar

Requisitos:
- gfortran e CMake ≥ 3.21;
- para a GPU: ROCm (com hipfort, hipSPARSE e hipBLAS) ou CUDA + hipfort/hipSPARSE compilados para CUDA. O hipfort precisa ter sido compilado com o mesmo gfortran;
- para os problemas CUTEst: o ambiente do CUTEst carregado.

```bash
source ../env_cutest.sh          # opcional; sem ele os alvos cutestma_* nao sao criados

cmake --preset cpu    && cmake --build --preset cpu      # so CPU (nao procura HIP)
cmake --preset amd    && cmake --build --preset amd      # CPU + GPU AMD (exige ROCm)
cmake --preset nvidia && cmake --build --preset nvidia   # CPU + GPU NVIDIA (exige CUDA)
```

Opções do CMake:
- **`GENCAN_GPU=AUTO|ON|OFF`:**
  - sem preset, vale `AUTO`: a GPU entra se o backend for encontrado; se não, compila só a CPU e avisa o que faltou;
  - os presets `amd`/`nvidia` usam `ON`, que aborta se faltar algo;
  - o preset `cpu` usa `OFF`.
- **`HIP_PLATFORM=amd|nvidia`:** o backend.
- **`KERNEL_FLAGS`:** flags extras do `hipcc`, por exemplo `-DKERNEL_FLAGS="--offload-arch=gfx1201"`.
- **Se algo não for encontrado:** `-DHIPFORT_ROOT=...` (ou `-DHIPFORT_MOD_DIR`/`-DHIPFORT_LIB`) e `-DCUDAToolkit_ROOT=...`.
- **`PROBLEM=<NOME>`:** o problema CUTEst. Ele é decodificado quando o CMake configura, então trocar de problema exige reconfigurar.

Executáveis (a versão `_gpu` só existe com a GPU ligada):
- `toyprob_cpu` e `toyprob_gpu`: o problema de exemplo.
- `cutestma_cpu` e `cutestma_gpu`: um problema CUTEst. As duas versões saem do mesmo `cutestma.f90`.
- `test_hess_csr` e `test_cgm_gpu`: testes sintéticos.

## Rodar

```bash
./build/amd/toyprob_gpu

# um problema CUTEst, CPU x GPU
cmake --preset amd -DPROBLEM=ROSENBR && cmake --build --preset amd
(cd build/amd && ./cutestma_cpu | tail -15 && ./cutestma_gpu | tail -15)

# lista de problemas (retomavel; pula o que esta em done_gencan.txt)
PRESET=amd ./run-cutest-gpu.sh pequena.txt
paste -d'\n' table_cpu.txt table_gpu.txt | less -S

# testes sinteticos
./run-tests.sh

# build cpu contra a referencia (exigido depois de mexer em lssit/cgm/gencan.f90)
./compare-cpu-ref.sh
```

Listas de problemas:
- `pequena.txt`: 7 problemas, para validação rápida.
- `grandes_bound.txt`: 18 problemas com caixa e `n ≥ 10000`.
- `problems.txt`: lista completa.

Para gerar outras listas, use `../CUTEst-aux/list_cutest_problems.py`, que aceita `--min-n` pelo `n` padrão do SIF.

O `run-cutest-gpu.sh` limita cada execução a 600 s de CPU (`ulimit -St`). Para medir desempenho em problemas grandes, aumente ou retire esse limite.

## Estado

- **AMD (RX 9070 XT):**
  - testes sintéticos OK;
  - build `cpu` idêntico à referência;
  - CUTEst `pequena.txt` com resultado idêntico ao da fase 5, que já tinha sido validada contra a CPU.
- **NVIDIA:** o build com CUDA ainda não foi testado com o código atual.
- **Fora do escopo por enquanto:** `truehp = .false.` (quocientes incrementais) e restrições gerais (Algencan). Nesses casos o Gencan usa o caminho da CPU.
- **Obsoletos:** `compile.sh` e `run-cutest.sh` deste diretório usam alvos antigos (`toyprob`, `cutestma`). Use o CMake.
