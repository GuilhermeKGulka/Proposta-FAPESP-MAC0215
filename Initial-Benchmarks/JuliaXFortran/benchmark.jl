using SparseArrays
using LinearAlgebra
using AMDGPU
using AMDGPU.rocSPARSE
using Plots

# 1. Implementação Própria (Ingênua) - Formato CSR
function spmv_naive(rowptr, colvals, nzvals, x, y)
    for i in 1:(length(rowptr)-1)
        sum = 0.0
        for j in rowptr[i]:(rowptr[i+1]-1)
            sum += nzvals[j] * x[colvals[j]]
        end
        y[i] = sum
    end
end

# 2. Configuração do Teste
# Mais pontos (tamanhos) para gerar um gráfico log-log com uma curva visível
sizes = [1000, 10000, 100000, 1000000, 10000000] 
samples = 100

# Arrays para armazenar os resultados de tempo (em milissegundos)
times_naive = Float64[]
times_cpu_opt = Float64[]
times_gpu_mul = Float64[]
times_gpu_mv = Float64[]

has_gpu = has_rocm_gpu()

println("Iniciando benchmarks...")

for n in sizes
    println("Avaliando dimensão: $n")
    
    # Gerar matriz esparsa aleatória e vetores
    d = fill(2.0, n); dl = fill(-1.0, n-1); du = fill(-1.0, n-1)
    A = sprand(n,n,10/n)#spdiagm(-1 => dl, 0 => d, 1 => du)
    x = rand(n)
    y = zeros(n)
    
    # Prepara o formato CSR: O formato CSC da matriz transposta 
    # contém exatamente os arrays CSR da matriz original.
    A_csr = SparseMatrixCSC(transpose(A))
    
    # --- Benchmark CPU Ingênuo ---
    # colptr age como rowptr, rowval age como colvals
    t_naive = minimum(@elapsed spmv_naive(A_csr.colptr, A_csr.rowval, A_csr.nzval, x, y) for _ in 1:samples)
    push!(times_naive, t_naive * 1000)
    
    # --- Benchmark CPU Otimizado (SuiteSparse/BLAS) ---
    t_cpu = minimum(@elapsed A * x for _ in 1:samples)
    push!(times_cpu_opt, t_cpu * 1000)
    
    # --- Benchmark GPU (AMD) ---
    if has_gpu
        A_gpu = AMDGPU.rocSPARSE.ROCSparseMatrixCSR(A)
        x_gpu = ROCArray(x)
        y_gpu = ROCArray(zeros(n))
        t_gpu_mul = minimum(@elapsed mul!(y_gpu, A_gpu, x_gpu) for _ in 1:samples)
        push!(times_gpu_mul, t_gpu_mul * 1000)

        y_gpu = ROCArray(zeros(n))
        t_gpu_mv = minimum(@elapsed mv!('N', Float64(1.0), A_gpu, x_gpu, Float64(0.0), y_gpu, 'O') for _ in 1:samples)
        push!(times_gpu_mv, t_gpu_mv * 1000)
    end
end

# 3. Geração do Gráfico (Replicando o comportamento do Matplotlib)
println("\nGerando gráfico de comparação...")

# Configuração do canvas com as mesmas características do Python
plt = plot(
    title="SpMV Wall Time Comparison",
    xaxis=:log10,  # x_scale='log'
    yaxis=:log10,  # y_scale='log'
    xticks = 10.0 .^ (3:7),
    yticks = 10.0 .^ (-3:4),
    xlabel="Problem size n (log scale)",
    ylabel="Wall time (ms) (log scale)",
    legend=:topleft,
    grid=true,
    gridalpha=0.5,     # alpha=0.5
    minorgrid=true,    # which='both'
    size=(800, 600),   # figsize=(8, 6) aproximado
    dpi=150
)

# Adicionando as linhas (equivalente ao 'o-' do matplotlib)
plot!(plt, sizes, times_naive, label="CPU Ingênuo", marker=:circle, linewidth=2)
plot!(plt, sizes, times_cpu_opt, label="CPU Otimizado", marker=:circle, linewidth=2)

if has_gpu
    plot!(plt, sizes, times_gpu_mul, label="GPU AMD-mul", marker=:circle, linewidth=2)
    plot!(plt, sizes, times_gpu_mv, label="GPU AMD-mv", marker=:circle, linewidth=2)
end

# Salvando a imagem
output_file = "wall_time_comparison-logxlog.png"
savefig(plt, output_file)
println("Gráfico salvo com sucesso em: $output_file")