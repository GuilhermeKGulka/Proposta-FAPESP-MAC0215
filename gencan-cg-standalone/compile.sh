#!/bin/bash
# Compiles the solver and the example program. Run ./toyprob afterwards.
set -e
cd sources
gfortran -O3 -c lssit.f90
gfortran -O3 -c cgm.f90
gfortran -O3 -c gencan.f90
gfortran -O3 -c toyprob.f90
gfortran -O3 lssit.o cgm.o gencan.o toyprob.o -o ../toyprob
rm -f *.o *.mod
cd ..
echo "Done. Run ./toyprob"
