# Gencan GPU — Aceleração do Gradiente Conjugado

Implementação do método de Lagrangianos Aumentados (Gencan) com o solver de Gradientes Conjugados (CG) executando na GPU via AMD HIP.

## Estrutura do Projeto
- `lssit.f90`, `gencan.f90`: Código original do Gencan (CPU).
- `cgm_gpu.f90`: O solver de Gradientes Conjugados implementado em GPU.
- `gencan_gpu_context.f90`: Gerenciamento de memória da GPU (vetores, handles, SpMV).
- `scalar_kernels.cpp`: Kernels HIP para reduções (ex: cálculo de `amax`).

## Como Compilar
```bash
# Requer o ambiente do CUTEst carregado
source ../../env_cutest.sh

# Compilar para AMD
cmake --preset amd && cmake --build --preset amd

# Compilar para NVIDIA
cmake --preset nvidia && cmake --build --preset nvidia