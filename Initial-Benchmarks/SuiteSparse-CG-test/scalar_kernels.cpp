#include <hip/hip_runtime.h>

// Calcula alpha = r2 / dot_pHp
__global__ void compute_alpha_kernel(const double* r2, const double* dot_pHp, double* alpha) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        *alpha = (*r2) / (*dot_pHp);
    }
}

// Calcula beta e faz o "lastr2 = r2" na GPU
__global__ void compute_beta_and_save_kernel(double* r2, double* lastr2, double* beta) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        *beta = (*r2) / (*lastr2);
        *lastr2 = *r2; 
    }
}

// Kernel de Vetor Fundido: p = -r + beta * p (Elimina o -1 do host e poupa banda de memória)
__global__ void update_p_vector_kernel(int p_n, double* p, const double* r, const double* beta) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < p_n) {
        p[idx] = -r[idx] + (*beta) * p[idx];
    }
}

// Interfaces de vinculação C/Fortran
extern "C" {
    void launch_compute_alpha(const double* d_r2, const double* d_dot_pHp, double* d_alpha) {
        hipLaunchKernelGGL(compute_alpha_kernel, dim3(1), dim3(1), 0, 0, d_r2, d_dot_pHp, d_alpha);
    }

    void launch_compute_beta(double* d_r2, double* d_lastr2, double* d_beta) {
        hipLaunchKernelGGL(compute_beta_and_save_kernel, dim3(1), dim3(1), 0, 0, d_r2, d_lastr2, d_beta);
    }

    void launch_update_p_vector(int p_n, double* d_p, const double* d_r, const double* d_beta) {
        int threads = 256;
        int blocks = (p_n + threads - 1) / threads;
        hipLaunchKernelGGL(update_p_vector_kernel, dim3(blocks), dim3(threads), 0, 0, p_n, d_p, d_r, d_beta);
    }
}