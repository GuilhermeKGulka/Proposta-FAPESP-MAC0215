#!/bin/bash
# Compara, bit a bit, a saida do build CMake "cpu" (fontes deste diretorio sem a
# macro GENCAN_GPU) com a da referencia gencan-cg-standalone/, compilada do jeito
# dela (gfortran -O3, como em compile.sh/run-cutest.sh). Exigido pelo CLAUDE.md
# depois de qualquer mudanca em lssit.f90, cgm.f90 ou gencan.f90.
#
#   ./compare-cpu-ref.sh               toyprob + problemas de pequena.txt
#   ./compare-cpu-ref.sh lista.txt     toyprob + problemas da lista
#
# Saidas e diffs em build/ref-cmp/. So a linha com "time" (tempo de CPU) e ignorada.
set -e
cd "$(dirname "$0")"
GPU=$(pwd)
REF=$GPU/../gencan-cg-standalone
LISTA=${1:-pequena.txt}
S=$GPU/build/ref-cmp
source ../env_cutest.sh > /dev/null
rm -rf "$S" && mkdir -p "$S/ref" "$S/out"
cp "$REF"/sources/*.f90 "$S/ref/"

# referencia: biblioteca + toyprob
cd "$S/ref"
gfortran -O3 -c lssit.f90 cgm.f90 gencan.f90
ar rcs libgencan.a lssit.o cgm.o gencan.o
gfortran -O3 -c toyprob.f90 && gfortran -O3 toyprob.o -L. -lgencan -o toyprob

falhas=0
./toyprob > "$S/out/toy_ref.txt" 2>&1
(cd "$GPU" && cmake --preset cpu > /dev/null && cmake --build --preset cpu --target toyprob_cpu > /dev/null)
"$GPU/build/cpu/toyprob_cpu" > "$S/out/toy_new.txt" 2>&1
if cmp -s "$S/out/toy_ref.txt" "$S/out/toy_new.txt"; then
  echo "toyprob: IDENTICO"
else
  echo "toyprob: DIFERENTE"; falhas=$((falhas+1))
fi

for p in $(grep -v '^#' "$GPU/$LISTA"); do
  # referencia
  mkdir -p "$S/ref/$p" && cd "$S/ref/$p"
  "$SIFDECODE/bin/sifdecoder" "$MASTSIF/$p.SIF" > /dev/null
  srcs="ELFUN.f GROUP.f RANGE.f"; [ -e EXTER.f ] && srcs="ELFUN.f EXTER.f GROUP.f RANGE.f"
  gfortran -c -O3 $srcs
  gfortran -c -O3 "$S/ref/cutestma.f90" -I"$S/ref"
  gfortran -O3 ${srcs//.f/.o} cutestma.o -L"$S/ref" -lgencan -L"$CUTEST/objects/$MYARCH/double" -lcutest -o cutestma
  ./cutestma > "$S/out/${p}_ref.txt" 2>&1 || true
  # build cpu deste diretorio
  (cd "$GPU" && cmake --preset cpu -DPROBLEM="$p" > /dev/null && cmake --build --preset cpu --target cutestma_cpu > /dev/null)
  (cd "$GPU/build/cpu" && ./cutestma_cpu > "$S/out/${p}_new.txt" 2>&1 || true)
  if diff <(grep -vi "time" "$S/out/${p}_ref.txt") <(grep -vi "time" "$S/out/${p}_new.txt") > "$S/out/${p}.diff"; then
    echo "$p: IDENTICO"
  else
    echo "$p: DIFERENTE (ver $S/out/${p}.diff)"; falhas=$((falhas+1))
  fi
done
echo "compare-cpu-ref: diferencas = $falhas"
