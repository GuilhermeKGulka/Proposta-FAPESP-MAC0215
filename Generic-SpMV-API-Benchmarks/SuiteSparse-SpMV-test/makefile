# Compilador e Flags
FC = hipfc
MODDIR = mod
OBJDIR = obj
FFLAGS = -O3 -J$(MODDIR) -I$(MODDIR)
LIBS = -lhipsparse

# Nome do executavel (na raiz do projeto)
EXEC = benchmark_suitesparse

# Arquivos fontes
SRCS = mmio.f90 benchmark_suitesparse.f90

# Objetos correspondentes no diretorio obj/
OBJS = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRCS))

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

# Regra de dependencia: benchmark depende do modulo mmio
# Isso garante que mmio.f90 seja compilado primeiro
$(OBJDIR)/benchmark_suitesparse.o: $(OBJDIR)/mmio.o

# Regra de compilacao generica para .f90 -> .o
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Limpeza do projeto
clean:
	rm -rf $(OBJDIR) $(MODDIR) $(EXEC) results.csv