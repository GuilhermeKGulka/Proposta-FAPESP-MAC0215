#!/bin/bash

# Runs Gencan on CUTEst problems.
#
#   ./run-cutest.sh                 all the problems in problems.txt
#   ./run-cutest.sh ROSENBR         a single problem
#   ./run-cutest.sh -f mylist.txt   the problems in mylist.txt, one name per line

# If CUTEST, SIFDECODE, MASTSIF and MYARCH are not already set in your
# environment, uncomment the four lines below and write your paths.
#export CUTEST=$HOME/CUTEst/cutest
#export SIFDECODE=$HOME/CUTEst/sifdecode
#export MASTSIF=$HOME/CUTEst/mastsif
#export MYARCH=pc64.lnx.gfo

# CPU-time limit per problem, in seconds
tlim=600

# Set to 1 to also append one line per problem to table-gencan.txt
table=0

# ------------------------------------------------------------------

if [ "$1" = "-f" ]; then
    problems=$(grep -v '^#' "$2")
elif [ -n "$1" ]; then
    problems="$1"
else
    problems=$(grep -v '^#' problems.txt)
fi

cd sources
gfortran -O3 -c lssit.f90
gfortran -O3 -c cgm.f90
gfortran -O3 -c gencan.f90
ar rcs libgencan.a lssit.o cgm.o gencan.o
cd ..

if [ $table -eq 1 ] && [ ! -f table-gencan.txt ]; then
    echo "# Problem      SC IERR      CPU(s)         n                        f    gpsupn     iter     fcnt     gcnt     hcnt" > table-gencan.txt
fi

for i in $problems; do
    echo "Running $i ..."

    $SIFDECODE/bin/sifdecoder $MASTSIF/$i.SIF

    if [ -e EXTER.f ]; then
        gfortran -c -O3 ELFUN.f EXTER.f GROUP.f RANGE.f
        gfortran -c -O3 sources/cutestma.f90 -Isources
        gfortran -O3 ELFUN.o EXTER.o GROUP.o RANGE.o cutestma.o \
                 -Lsources -lgencan -L$CUTEST/objects/$MYARCH/double -lcutest -o cutestma
    else
        gfortran -c -O3 ELFUN.f GROUP.f RANGE.f
        gfortran -c -O3 sources/cutestma.f90 -Isources
        gfortran -O3 ELFUN.o GROUP.o RANGE.o cutestma.o \
                 -Lsources -lgencan -L$CUTEST/objects/$MYARCH/double -lcutest -o cutestma
    fi

    ulimit -St $tlim
    ./cutestma
    ulimit -St unlimited

    if [ $table -eq 1 ]; then
        if [ -e tabline.txt ]; then
            cat tabline.txt >> table-gencan.txt
        elif [ -e solver-interrupted-tabline.txt ]; then
            printf " %-14s 99   0 %10.3f" $i $tlim >> table-gencan.txt
            cat solver-interrupted-tabline.txt >> table-gencan.txt
        else
            printf " %-14s 99   0 %10.3f  (no information)\n" $i $tlim >> table-gencan.txt
        fi
    fi

    rm -f AUTOMAT.d OUTSDIF.d ELFUN.f ELFUN.o EXTER.f EXTER.o GROUP.f GROUP.o RANGE.f RANGE.o
    rm -f cutestma.o cutestma tabline.txt solver-interrupted-tabline.txt
done

rm -f sources/*.o sources/*.mod sources/libgencan.a
