# Benchmark Newton-CG truncado (CUTEst) — HIP/Fortran

Compara o Newton truncado com CG na GPU (hipSPARSE generico + hipBLAS + kernels escalares HIP)
contra o CG na CPU, em problemas irrestritos do CUTEst. Funciona nos backends AMD (ROCm) e
NVIDIA (CUDA) via hipfort. O SpMV do CG usa a API generica do hipSPARSE
(`bufferSize` + `preprocess` + `SpMV`): descritores, buffer e preprocess sao criados **uma vez**
(padrao e nnz da Hessiana sao fixos; so os valores mudam, no mesmo ponteiro de device).

## Requisitos

- CMake >= 3.21, gfortran e `hipcc` (para `scalar_kernels.cpp`)
- hipfort (`libhipfort-amdgcn.a` / `libhipfort-nvptx.a` e os `.mod`), compilado pelo **mesmo** gfortran
- AMD: ROCm (hip, hipsparse, hipblas). NVIDIA: CUDA Toolkit + hipSPARSE/hipBLAS compilados para CUDA
- CUTEst, SIFDecode, ARCHDefs (submodulos do repositorio) e os `.SIF` em `externals/mastsif`,
  com o ambiente carregado por `source <raiz do repositorio>/env_cutest.sh`

## Uso

```bash
source ../../env_cutest.sh               # CUTEST, SIFDECODE, MASTSIF, MYARCH
cp ../../Initial-Benchmarks/CUTEst-Newton-CG-test/problems.txt .   # lista de problemas (uma vez)

PRESET=amd ./run_benchmark.sh            # ou PRESET=nvidia (export CMAKE_PREFIX_PATH=... antes)
```

O script compila a parte fixa uma vez e, para cada problema, reconfigura com
`cmake --preset <p> -DPROBLEM=<NOME>` (isso roda o `sifdecoder` e gera `ELFUN.f`, `GROUP.f`,
`RANGE.f` e `EXTER.f` em `build/<preset>/`), compila e executa com limite de 600 s de CPU.
Um problema so e decodificado de novo se `PROBLEM` mudar.

Um problema individual:
```bash
cmake --preset amd -DPROBLEM=ROSENBR && cmake --build --preset amd --target benchmark_newton_cg
(cd build/amd && ./benchmark_newton_cg [nopre])
```

- Retomada: cada problema executado e anotado em `done_problems.txt` (nome, codigo de saida) e
  pulado nas proximas execucoes; apague a linha para repetir. Logs em `build/<preset>/logs/<problema>.log`.
- `NOPRE=1` (ou `nopre` no executavel) pula o `SpMV_preprocess`.
- Resultados em `results_newton_cg.csv` (variavel `RESULTS_CSV` muda o caminho). Colunas novas no fim:
  `t_GPU_prep_ms` (setup do SpMV; entra em `t_total_GPU_ms`, nao em `t_cg_GPU_ms`; `-1` se o preprocess
  falhou), `Gnorm_CPU` / `Gnorm_GPU` (norma do gradiente no ponto final; `-1` se a avaliacao falhou).

## Problemas comuns

- `Variavel de ambiente ... nao definida`: rode `source <raiz>/env_cutest.sh`.
- `sifdecoder falhou`: confira se `externals/mastsif/<PROBLEMA>.SIF` existe.
- `HIPSPARSE ERROR em <chamada>: codigo 3` em `preprocess`/`SpMV`: com `bufferSize = 0` o buffer externo
  tem de ser nulo, e so deve haver um preprocess por matriz (o modulo `rocm_context` ja cuida disso).
