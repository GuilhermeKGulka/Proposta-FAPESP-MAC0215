# =====================================================================
# Makefile — Benchmark Newton-CG (CPU vs GPU) em problemas CUTEst
# =====================================================================

PROBLEM ?= ROSENBR

# ---- toolchain ------------------------------------------------------
FC        := hipfc
CXX       := hipcc
FC_FLAGS  := -O3
CXX_FLAGS := -O3

SIFDECODE_BIN := $(SIFDECODE)/bin/sifdecoder
SIF_FILE      := $(MASTSIF)/$(PROBLEM).SIF

# ---- CUTEst (build in-place) ----------------------------------------
CUTEST_LIB := $(CUTEST)/objects/$(MYARCH)/double
# Modulos ficam fora de objects/; ajuste se necessario apos o `find`
CUTEST_MOD := $(CUTEST)/modules/$(MYARCH)/double

CUTEST_INC := -I$(CUTEST_MOD)
LDLIBS     := -L$(CUTEST_LIB) -lcutest -lhipsparse -lhipblas

# ---- fontes fixas ---------------------------------------------------
FIXED_F90 := rocm_context.f90 benchmark_cg_mod.f90 cutest_iface.f90 newton_cg_mod.f90
FIXED_OBJ := $(FIXED_F90:.f90=.o)
KERNEL    := scalar_kernels.o

PROG      := benchmark_newton_cg
PROG_OBJ  := $(PROG).o

DECODED_F := ELFUN.f GROUP.f RANGE.f
DECODED_O := $(DECODED_F:.f=.o)
EXTER_F   := $(wildcard EXTER.f)
EXTER_O   := $(EXTER_F:.f=.o)
PROBLEM_OBJ := $(DECODED_O) $(EXTER_O)

.PHONY: all run decode build clean distclean check-env

all: build

run:
	@$(MAKE) --no-print-directory decode
	@$(MAKE) --no-print-directory build
	@bash -c 'ulimit -St 600; ./$(PROG)'

build: check-env $(PROG)

decode: check-env
	@if [ ! -f ELFUN.f ] || \
	   [ "$$(cat .last_problem 2>/dev/null)" != "$(PROBLEM)" ]; then \
	    echo ">>> sifdecoder $(SIF_FILE)"; \
	    $(SIFDECODE_BIN) "$(SIF_FILE)" || exit 1; \
	    echo "$(PROBLEM)" > .last_problem; \
	fi

check-env:
	@if [ -z "$(CUTEST)" ] || [ -z "$(SIFDECODE)" ] || \
	    [ -z "$(MASTSIF)" ] || [ -z "$(MYARCH)" ]; then \
	    echo "Erro: defina CUTEST, SIFDECODE, MASTSIF e MYARCH (source env_cutest.sh)."; \
	    exit 1; \
	fi
	@if [ ! -f "$(CUTEST_LIB)/libcutest.a" ]; then \
	    echo "Erro: libcutest.a nao encontrada em $(CUTEST_LIB)/"; \
	    exit 1; \
	fi

rocm_context.o:      rocm_context.f90
cutest_iface.o:      cutest_iface.f90
benchmark_cg_mod.o:  benchmark_cg_mod.f90 rocm_context.o

newton_cg_mod.o: newton_cg_mod.f90 \
                 cutest_iface.o benchmark_cg_mod.o rocm_context.o

benchmark_newton_cg.o: benchmark_newton_cg.f90 \
                       cutest_iface.o benchmark_cg_mod.o rocm_context.o \
                       newton_cg_mod.o

%.o: %.f90
	$(FC) $(FC_FLAGS) $(CUTEST_INC) -c $< -o $@

%.o: %.f
	$(FC) $(FC_FLAGS) -c $< -o $@

scalar_kernels.o: scalar_kernels.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

$(PROG): $(PROG_OBJ) $(FIXED_OBJ) $(KERNEL) $(PROBLEM_OBJ)
	$(FC) $(FC_FLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f AUTOMAT.d OUTSDIF.d \
	      ELFUN.f ELFUN.o EXTER.f EXTER.o GROUP.f GROUP.o RANGE.f RANGE.o \
	      .last_problem $(PROG) $(PROG_OBJ)

distclean: clean
	rm -f $(FIXED_OBJ) $(KERNEL) *.mod