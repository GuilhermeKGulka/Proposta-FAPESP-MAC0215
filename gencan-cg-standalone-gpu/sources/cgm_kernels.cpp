#include <hip/hip_runtime.h>
#include <cfloat>
#include <cstdio>
#include <cstdlib>

#define BS   256
#define MAXB 1024
#define HIPCHK() do { hipError_t e_ = hipGetLastError(); if (e_ != hipSuccess) { \
  fprintf(stderr, "cgm_kernels: %s (%s:%d)\n", hipGetErrorString(e_), __FILE__, __LINE__); exit(1); } } while (0)

// Reducao de bloco: NS somas e NM minimos por thread -> resultado na thread 0.
template <int NS, int NM>
__device__ __forceinline__ void block_reduce(double* s, double* m, double* sh) {
  const int t = threadIdx.x;
#pragma unroll
  for (int k = 0; k < NS; k++) sh[k * BS + t] = s[k];
#pragma unroll
  for (int k = 0; k < NM; k++) sh[(NS + k) * BS + t] = m[k];
  __syncthreads();
  for (int off = BS / 2; off > 0; off >>= 1) {
    if (t < off) {
#pragma unroll
      for (int k = 0; k < NS; k++) sh[k * BS + t] += sh[k * BS + t + off];
#pragma unroll
      for (int k = 0; k < NM; k++)
        sh[(NS + k) * BS + t] = fmin(sh[(NS + k) * BS + t], sh[(NS + k) * BS + t + off]);
    }
    __syncthreads();
  }
  if (t == 0) {
#pragma unroll
    for (int k = 0; k < NS; k++) s[k] = sh[k * BS];
#pragma unroll
    for (int k = 0; k < NM; k++) m[k] = sh[(NS + k) * BS];
  }
  __syncthreads();
}

// Reducao de grade em um unico lancamento: o ultimo bloco a terminar faz a reducao final.
// Saida: out[0..NS-1] = somas, out[NS..NS+NM-1] = minimos.
// O contador (unsigned, zerado na criacao) se reinicia sozinho via atomicInc.
template <int NS, int NM>
__device__ void grid_finish(double* s, double* m, double* sh,
                            double* part, unsigned int* counter, double* out) {
  constexpr int K = NS + NM;
  __shared__ int last;
  block_reduce<NS, NM>(s, m, sh);
  if (threadIdx.x == 0) {
    for (int k = 0; k < NS; k++) part[blockIdx.x * K + k] = s[k];
    for (int k = 0; k < NM; k++) part[blockIdx.x * K + NS + k] = m[k];
    __threadfence();
    unsigned int ticket = atomicInc(counter, gridDim.x - 1);
    last = (ticket == gridDim.x - 1);
  }
  __syncthreads();
  if (!last) return;

  volatile const double* vp = part;   // evita ler linhas de cache antigas
  for (int k = 0; k < NS; k++) s[k] = 0.0;
  for (int k = 0; k < NM; k++) m[k] = DBL_MAX;
  for (int b = threadIdx.x; b < (int)gridDim.x; b += BS) {
    for (int k = 0; k < NS; k++) s[k] += vp[b * K + k];
    for (int k = 0; k < NM; k++) m[k] = fmin(m[k], vp[b * K + NS + k]);
  }
  block_reduce<NS, NM>(s, m, sh);
  if (threadIdx.x == 0) {
    for (int k = 0; k < NS; k++) out[k] = s[k];
    for (int k = 0; k < NM; k++) out[NS + k] = m[k];
  }
}

// p = -r + coef * p     (z = r, sem precondicionador)
__global__ __launch_bounds__(BS) void k_update_p(int n, double coef, const double* r, double* p) {
  for (int i = blockIdx.x * BS + threadIdx.x; i < n; i += gridDim.x * BS)
    p[i] = -r[i] + coef * p[i];
}

// Reducoes que dependem so de (p, hp, d, g, r):
//  out[0]=p.hp  out[1]=p.d  out[2]=g.p  out[3]=p.r  out[4]=p.p
//  out[5]=amax para +p   out[6]=amax para -p   (o host escolhe conforme o sinal de p.r)
__global__ __launch_bounds__(BS) void k_dots1(int n, double delta,
    const double* p, const double* hp, const double* d, const double* g, const double* r,
    double* part, unsigned int* counter, double* out) {
  __shared__ double sh[7 * BS];
  double s[5] = {0, 0, 0, 0, 0};
  double m[2] = {DBL_MAX, DBL_MAX};
  for (int i = blockIdx.x * BS + threadIdx.x; i < n; i += gridDim.x * BS) {
    const double pi = p[i], di = d[i];
    s[0] += pi * hp[i];
    s[1] += pi * di;
    s[2] += g[i] * pi;
    s[3] += pi * r[i];
    s[4] += pi * pi;
    if (pi > 0.0) {
      m[0] = fmin(m[0], (delta - di) / pi);
      m[1] = fmin(m[1], (-delta - di) / (-pi));
    } else if (pi < 0.0) {
      m[0] = fmin(m[0], (-delta - di) / pi);
      m[1] = fmin(m[1], (delta - di) / (-pi));
    }
  }
  grid_finish<5, 2>(s, m, sh, part, counter, out);
}

// d += c p ; r += c hp ; e reducoes sobre o NOVO d e o NOVO r:
//  out[0]=r.r   out[1]=numero de i que violam "samep"   out[2]=abox
__global__ __launch_bounds__(BS) void k_update(int n, double coef, double macheps,
    const double* p, const double* hp, double* d, double* r,
    const double* xf, const double* lbf, const double* ubf,
    const int* lindf, const int* uindf,
    double* part, unsigned int* counter, double* out) {
  __shared__ double sh[3 * BS];
  double s[2] = {0, 0};
  double m[1] = {DBL_MAX};
  for (int i = blockIdx.x * BS + threadIdx.x; i < n; i += gridDim.x * BS) {
    const double step = coef * p[i];
    const double di = d[i] + step;
    const double ri = r[i] + coef * hp[i];
    d[i] = di;
    r[i] = ri;
    s[0] += ri * ri;
    if (fabs(step) > macheps * fmax(1.0, fabs(di))) s[1] += 1.0;
    if (di > 0.0 && uindf[i])      m[0] = fmin(m[0], (ubf[i] - xf[i]) / di);
    else if (di < 0.0 && lindf[i]) m[0] = fmin(m[0], (lbf[i] - xf[i]) / di);
  }
  grid_finish<2, 1>(s, m, sh, part, counter, out);
}

static inline int nblocks(int n) {
  int b = (n + BS - 1) / BS;
  if (b < 1) b = 1;
  if (b > MAXB) b = MAXB;
  return b;
}

extern "C" {
void launch_cgm_update_p(int n, double coef, const double* r, double* p) {
  hipLaunchKernelGGL(k_update_p, dim3(nblocks(n)), dim3(BS), 0, 0, n, coef, r, p);
  HIPCHK();
}
void launch_cgm_dots1(int n, double delta, const double* p, const double* hp, const double* d,
                      const double* g, const double* r,
                      double* part, unsigned int* counter, double* out) {
  hipLaunchKernelGGL(k_dots1, dim3(nblocks(n)), dim3(BS), 0, 0,
                     n, delta, p, hp, d, g, r, part, counter, out);
  HIPCHK();
}
void launch_cgm_update(int n, double coef, double macheps, const double* p, const double* hp,
                       double* d, double* r, const double* xf, const double* lbf, const double* ubf,
                       const int* lindf, const int* uindf,
                       double* part, unsigned int* counter, double* out) {
  hipLaunchKernelGGL(k_update, dim3(nblocks(n)), dim3(BS), 0, 0,
                     n, coef, macheps, p, hp, d, r, xf, lbf, ubf, lindf, uindf, part, counter, out);
  HIPCHK();
}
}