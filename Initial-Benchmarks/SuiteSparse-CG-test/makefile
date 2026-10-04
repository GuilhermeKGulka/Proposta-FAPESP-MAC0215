# Compiladores e Flags
FC = hipfc
CXX = hipcc
MODDIR = mod
OBJDIR = obj
FFLAGS = -O3 -J$(MODDIR) -I$(MODDIR)
CXXFLAGS = -O3
LIBS = -lhipsparse -lhipblas

# Nome do executavel (na raiz do projeto)
EXEC = benchmark_cg

# Objetos correspondentes no diretorio obj/
OBJS = $(OBJDIR)/mmio.o \
       $(OBJDIR)/rocm_context.o \
       $(OBJDIR)/scalar_kernels.o \
       $(OBJDIR)/benchmark_cg.o

# Alvos phony (nao sao arquivos reais)
.PHONY: all clean dirs

all: dirs $(EXEC)

# Criacao dos diretorios
dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(MODDIR)

# Linkagem do executavel final
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# Regras de dependencia
# benchmark_cg depende dos modulos mmio e rocm_context
$(OBJDIR)/benchmark_cg.o: $(OBJDIR)/mmio_2.o $(OBJDIR)/rocm_context.o

# Regra de compilacao para Fortran (.f90 -> .o)
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Regra de compilacao para C/C++ HIP (.cpp -> .o)
$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Limpeza do projeto
clean:
	rm -rf $(OBJDIR) $(MODDIR) $(EXEC) results_cg.csv