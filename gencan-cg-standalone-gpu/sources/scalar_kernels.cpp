#include <hip/hip_runtime.h>

__global__ void compute_amax_abox_kernel(int nfr, double delta, const double* p, const double* d, 
                                         const double* x, const double* lbnd, const double* ubnd, 
                                         const int* lind, const int* uind, const int* ifr, 
                                         double* amax_out, double* abox_out) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nfr) {
        double amax_val = 1e16;
        double abox_val = 1e16;
        
        if (p[idx] > 0.0) amax_val = (delta - d[idx]) / p[idx];
        else if (p[idx] < 0.0) amax_val = (-delta - d[idx]) / p[idx];
        
        int global_i = ifr[idx] - 1; 
        if (d[idx] > 0.0 && uind[global_i]) abox_val = (ubnd[global_i] - x[global_i]) / d[idx];
        else if (d[idx] < 0.0 && lind[global_i]) abox_val = (lbnd[global_i] - x[global_i]) / d[idx];
        
        atomicMin((unsigned long long int*)amax_out, __double_as_longlong(amax_val));
        atomicMin((unsigned long long int*)abox_out, __double_as_longlong(abox_val));
    }
}

__global__ void update_d_r_kernel(int nfr, double alpha, const double* p, const double* hp, double* d, double* r) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nfr) {
        d[idx] += alpha * p[idx];
        r[idx] += alpha * hp[idx];
    }
}

__global__ void update_p_kernel(int nfr, double bbeta, const double* z, double* p) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nfr) {
        p[idx] = -z[idx] + bbeta * p[idx];
    }
}

extern "C" {
    void launch_compute_amax_abox(int nfr, double* delta_ptr, double* p, double* d, double* x, 
                                  double* lbnd, double* ubnd, int* lind, int* uind, int* ifr, 
                                  double* amax, double* abox) {
        double init_val = 1e16;
        (void) hipMemcpy(amax, &init_val, sizeof(double), hipMemcpyHostToDevice);
        (void) hipMemcpy(abox, &init_val, sizeof(double), hipMemcpyHostToDevice);
        
        double delta_host;
        (void) hipMemcpy(&delta_host, delta_ptr, sizeof(double), hipMemcpyDeviceToHost);
        
        int threads = 256;
        int blocks = (nfr + threads - 1) / threads;
        hipLaunchKernelGGL(compute_amax_abox_kernel, dim3(blocks), dim3(threads), 0, 0, 
                           nfr, delta_host, p, d, x, lbnd, ubnd, lind, uind, ifr, amax, abox);
    }
    
    void launch_update_d_r(int nfr, double* alpha_ptr, double* p, double* hp, double* d, double* r) {
        double alpha_host;
        (void) hipMemcpy(&alpha_host, alpha_ptr, sizeof(double), hipMemcpyDeviceToHost);
        int threads = 256;
        int blocks = (nfr + threads - 1) / threads;
        hipLaunchKernelGGL(update_d_r_kernel, dim3(blocks), dim3(threads), 0, 0, nfr, alpha_host, p, hp, d, r);
    }

    void launch_update_p_vector(int nfr, double* p, const double* z, double bbeta) {
        int threads = 256;
        int blocks = (nfr + threads - 1) / threads;
        hipLaunchKernelGGL(update_p_kernel, dim3(blocks), dim3(threads), 0, 0, nfr, bbeta, z, p);
    }
}