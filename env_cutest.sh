#!/bin/bash
# env_cutest.sh — carregar com:  source ./env_cutest.sh

if [ -n "${BASH_SOURCE[0]}" ]; then
    _env_script="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION}" ]; then
    _env_script="${(%):-%x}"
else
    _env_script="$0"
fi
ROOT_DIR="$(cd "$(dirname "$_env_script")" && pwd)"

if [ "${BASH_SOURCE[0]}" = "$0" ] 2>/dev/null; then
    echo "AVISO: rode 'source $ROOT_DIR/env_cutest.sh', nao './env_cutest.sh'." >&2
fi

export ROOT_DIR
export ARCHDEFS="${ROOT_DIR}/externals/ARCHDefs"
export SIFDECODE="${ROOT_DIR}/externals/SIFDecode"
export CUTEST="${ROOT_DIR}/externals/CUTEst"
export MASTSIF="${ROOT_DIR}/externals/mastsif"
export MYARCH="pc64.lnx.gfa"

# Nao existe mycutest com prefixo: o build foi in-place.
# Definimos MYCST=CUTEST para compatibilidade com scripts que a citam.
export MYCST="${CUTEST}"

# sifdecoder esta em SIFDecode/bin (script wrapper)
export PATH="${SIFDECODE}/bin:${PATH}"

# Sanidade
[ -f "${CUTEST}/objects/${MYARCH}/double/libcutest.a" ] || \
    echo "AVISO: libcutest.a nao encontrada em ${CUTEST}/objects/${MYARCH}/double/" >&2
[ -x "${SIFDECODE}/bin/sifdecoder" ] || \
    echo "AVISO: sifdecoder nao executavel em ${SIFDECODE}/bin/" >&2

echo "Ambiente CUTEst:"
echo "  CUTEST    = $CUTEST"
echo "  SIFDECODE = $SIFDECODE"
echo "  MYARCH    = $MYARCH"