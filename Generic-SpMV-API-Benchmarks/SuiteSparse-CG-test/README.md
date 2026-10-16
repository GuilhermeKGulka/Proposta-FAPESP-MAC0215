# Benchmark CG (SuiteSparse) — HIP/Fortran

Compara o gradiente conjugado na GPU (hipSPARSE generico + hipBLAS + kernels escalares HIP)
contra o CG na CPU, para matrizes Matrix Market (rhs `g = -1`, `x0 = 0`, `eps = 1e-8`).
Funciona nos backends AMD (ROCm) e NVIDIA (CUDA) via hipfort. O SpMV usa a API generica
do hipSPARSE (`bufferSize` + `preprocess` + `SpMV`), nao mais o `csrmv` legado.

## Requisitos

- CMake >= 3.21, gfortran e `hipcc` (para os kernels em `scalar_kernels.cpp`)
- hipfort (`libhipfort-amdgcn.a` / `libhipfort-nvptx.a` e os `.mod`), compilado pelo **mesmo** gfortran
- AMD: ROCm (hip, hipsparse, hipblas). NVIDIA: CUDA Toolkit (cudart, cusparse, cublas) + hipSPARSE/hipBLAS compilados para CUDA

## Compilar

```bash
# AMD (ROCm em /opt/rocm)
cmake --preset amd && cmake --build --preset amd

# NVIDIA (informe onde esta o hipfort/hipSPARSE/hipBLAS instalado)
export CMAKE_PREFIX_PATH=$HOME/hip-nv      # exemplo
cmake --preset nvidia && cmake --build --preset nvidia
```

Arquitetura dos kernels (opcional): `-DKERNEL_FLAGS="--offload-arch=gfx1201"` (AMD).
Se o CUDA nao for achado: `-DCUDAToolkit_ROOT=/usr/local/cuda`.

## Rodar

```bash
./build/amd/benchmark_cg caminho/matriz.mtx [nopre]     # uma matriz, anexa em results_cg.csv
BENCH=build/nvidia/benchmark_cg ./run_benchmark.sh      # varredura de ../../SuiteSparse/MM
```

`nopre` pula o `SpMV_preprocess`. Colunas novas no fim do CSV: `t_GPU_prep_ms` (setup do SpMV,
fora de `t_GPU_CG_ms`; `-1` se o preprocess falhou) e `Res_rel_GPU` / `Res_rel_CPU`
(residuo verdadeiro `||g + H d|| / max(1,||g||)` das duas solucoes).
Matrizes nao quadradas sao ignoradas (`IGNORADA`, sem linha no CSV).

## Problemas comuns

- `hipfort ... nao encontrado`: aponte `-DHIPFORT_ROOT=<prefixo>` ou passe `-DHIPFORT_MOD_DIR` e `-DHIPFORT_LIB`.
- Erro de versao em `.mod`: hipfort foi compilado com outra versao do gfortran.
- `HIPSPARSE ERROR em <chamada>: codigo 3`: se aparecer em `preprocess`/`SpMV`, lembre que com
  `bufferSize = 0` o buffer externo tem de ser **nulo** (o modulo `hipsparse_generic` ja cuida disso).
