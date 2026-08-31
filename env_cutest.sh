#!/bin/bash

# Pega o caminho absoluto do diretório atual
export ROOT_DIR=$(pwd)

# Variáveis exigidas pelo CUTEst
export ARCHDEFS="${ROOT_DIR}/externals/ARCHDefs"
export SIFDECODE="${ROOT_DIR}/externals/SIFDecode"
export CUTEST="${ROOT_DIR}/externals/CUTEst"

# Diretório onde os binários e bibliotecas compiladas vão ficar (isolado)
export MYCST="${ROOT_DIR}/externals/mycutest"

# Diretório onde você colocará os arquivos .SIF (problemas de teste)
export MASTSIF="${ROOT_DIR}/externals/mastsif"

# Atualiza PATH e MANPATH
export PATH="${MYCST}/bin:${PATH}"
export MANPATH="${MYCST}/man:${MANPATH}"

echo "Ambiente CUTEst configurado em: $MYCST"