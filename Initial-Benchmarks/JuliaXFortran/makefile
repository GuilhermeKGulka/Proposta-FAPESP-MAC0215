FC = hipfc
FFLAGS = -O3 -g -Jmod -Imod -lhipsparse
TARGET = test_perf
SOURCES = test_perf-2.f90

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(FC) $(FFLAGS) $(SOURCES) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

run-gpu-monitor: $(TARGET)
	@echo "Monitorando GPU durante execução..."
	@(watch -n 1 rocm-smi) &
	@./$(TARGET)
	@pkill watch

check-memory: $(TARGET)
	rocm-smi --showmeminfo vram

clean:
	rm -f $(TARGET) *.o *.mod

.PHONY: all run run-gpu-monitor profile clean