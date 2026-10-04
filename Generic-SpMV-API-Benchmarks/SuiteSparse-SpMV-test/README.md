# Benchmark SpMV (SuiteSparse) — HIP/Fortran

Compara `hipsparse` na GPU contra um SpMV CSR na CPU para matrizes Matrix Market.
Funciona nos backends AMD (ROCm) e NVIDIA (CUDA) via hipfort.

## Requisitos

- CMake >= 3.21 e gfortran
- hipfort (`libhipfort-amdgcn.a` / `libhipfort-nvptx.a` e os `.mod`), compilado pelo **mesmo** gfortran
- AMD: ROCm (hip, hipsparse). NVIDIA: CUDA Toolkit (cudart, cusparse) + hipfort/hipSPARSE para o backend NVIDIA

## Compilar

```bash
# AMD (ROCm em /opt/rocm)
cmake --preset amd && cmake --build --preset amd

# NVIDIA (informe onde esta o hipfort/hipSPARSE instalado)
export CMAKE_PREFIX_PATH=$HOME/hip-nv      # exemplo
cmake --preset nvidia && cmake --build --preset nvidia
```

O prefixo padrao do preset `amd` e `/opt/rocm`; para outro local use
`cmake --preset amd -DCMAKE_PREFIX_PATH=/outro/prefixo`.
Se o CUDA nao for achado: `-DCUDAToolkit_ROOT=/usr/local/cuda`.

## Rodar

```bash
./build/amd/benchmark_suitesparse caminho/matriz.mtx      # uma matriz, anexa em results.csv
BENCH=build/nvidia/benchmark_suitesparse ./run_benchmark.sh   # varredura de ../../SuiteSparse/MM
python plot_results.py
```

## Problemas comuns

- `hipfort ... nao encontrado`: aponte `-DHIPFORT_ROOT=<prefixo>` (ex.: `/usr/local`) ou passe `-DHIPFORT_MOD_DIR=<pasta dos .mod>` e `-DHIPFORT_LIB=<libhipfort-*.a>`.
- Erro de versao em `.mod`: hipfort foi compilado com outra versao do gfortran.
- Undefined reference a `stdc++` no link NVIDIA: adicione `-DCMAKE_EXE_LINKER_FLAGS=-lstdc++`.
- NVIDIA: `hipsparseDcsrmv` (API legada) nao existe mais no cuSPARSE atual; ate migrar para a API generica (SpMV) o link do benchmark falha nesse simbolo.