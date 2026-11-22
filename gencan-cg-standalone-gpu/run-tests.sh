#!/bin/bash
set -e
source ../env_cutest.sh          # opcional para os testes; mantem o PROBLEM do cache valido
PRESET="${PRESET:-amd}"
cmake --preset "$PRESET" > /dev/null
cmake --build --preset "$PRESET" --target test_hess_csr test_cgm_gpu
echo "######## etapa 1: Hessiana CSR ########"
./build/$PRESET/test_hess_csr
echo "######## etapas 3 e 4: cgm ########"
./build/$PRESET/test_cgm_gpu