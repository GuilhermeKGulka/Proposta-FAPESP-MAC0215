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

// Reducao de grade em UM lancamento: o ultimo bloco a terminar faz a reducao final.
// Saida: out[0..NS-1] = somas, out[NS..NS+NM-1] = minimos.
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

  volatile const double* vp = part;
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

// out[2]=sum r^2/D  out[3]=s.r  out[4]=u.r      (mesmo layout de k_update)
__global__ __launch_bounds__(BS) void k_dots0(int n, const double* r, const double* D,
    const double* s, const double* u, double* part, unsigned int* counter, double* out) {
  __shared__ double sh[5 * BS];
  double sm[5] = {0, 0, 0, 0, 0};
  double m[1] = {DBL_MAX};
  for (int i = blockIdx.x * BS + threadIdx.x; i < n; i += gridDim.x * BS) {
    const double ri = r[i];
    sm[2] += ri * ri / D[i];
    sm[3] += s[i] * ri;
    sm[4] += u[i] * ri;
  }
  grid_finish<5, 0>(sm, m, sh, part, counter, out);
}

// p = -z + coef p.  Escalares lidos do device (red), escritos por k_update / k_dots0:
//   red[0]=r.r  red[2]=sum r^2/D  red[3]=s.r  red[4]=u.r
// z = r (sem precond)  ou  z = r/D + c1 u + c2 s
__global__ __launch_bounds__(BS) void k_update_p_dev(int n, int use_pre, int corr,
    double sty, double psmdyty, double sgnprev, double ztrprev, int iter0,
    const double* red, const double* r, const double* D, const double* u, const double* s,
    double* p) {
  double c1 = 0.0, c2 = 0.0, ztr;
  if (use_pre) {
    const double rDr = red[2], sr = red[3], ur = red[4];
    if (corr) {
      c1 = sr / sty;
      c2 = ur / sty - psmdyty * sr / (sty * sty);
    }
    ztr = rDr + c1 * ur + c2 * sr;
  } else {
    ztr = red[0];
  }
  const double coef = iter0 ? 0.0 : sgnprev * (ztr / ztrprev);
  for (int i = blockIdx.x * BS + threadIdx.x; i < n; i += gridDim.x * BS) {
    if (use_pre) p[i] = -(r[i] / D[i] + c1 * u[i] + c2 * s[i]) + coef * p[i];
    else         p[i] = -r[i] + coef * p[i];
  }
}

// out[0]=p.hp out[1]=p.d out[2]=g.p out[3]=p.r out[4]=p.p out[5]=amax(+p) out[6]=amax(-p)
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

// d += c p ; r += c hp ; sobre o NOVO d e NOVO r:
// out[0]=r.r out[1]=#violacoes de samep out[2]=sum r^2/D out[3]=s.r out[4]=u.r out[5]=abox
// (out[2..4] so se dg != nullptr)
__global__ __launch_bounds__(BS) void k_update(int n, double coef, double macheps,
    const double* p, const double* hp, double* d, double* r,
    const double* xf, const double* lbf, const double* ubf,
    const int* lindf, const int* uindf,
    const double* dg, const double* sv, const double* uv,
    double* part, unsigned int* counter, double* out) {
  __shared__ double sh[6 * BS];
  double s[5] = {0, 0, 0, 0, 0};
  double m[1] = {DBL_MAX};
  for (int i = blockIdx.x * BS + threadIdx.x; i < n; i += gridDim.x * BS) {
    const double step = coef * p[i];
    const double di = d[i] + step;
    const double ri = r[i] + coef * hp[i];
    d[i] = di;
    r[i] = ri;
    s[0] += ri * ri;
    if (fabs(step) > macheps * fmax(1.0, fabs(di))) s[1] += 1.0;
    if (dg) {
      s[2] += ri * ri / dg[i];
      s[3] += sv[i] * ri;
      s[4] += uv[i] * ri;
    }
    if (di > 0.0 && uindf[i])      m[0] = fmin(m[0], (ubf[i] - xf[i]) / di);
    else if (di < 0.0 && lindf[i]) m[0] = fmin(m[0], (lbf[i] - xf[i]) / di);
  }
  grid_finish<5, 1>(s, m, sh, part, counter, out);
}

// val[pos-1] = v[i]   (scatter COO -> CSR; posicoes distintas, sem corrida)
__global__ void k_scatter(int hne, const double* v, const int* p1, const int* p2, double* val) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < hne) {
    int a = p1[i];
    if (a > 0) {
      double x = v[i];
      val[a - 1] = x;
      int b = p2[i];
      if (b > 0) val[b - 1] = x;
    }
  }
}

static inline int nblocks(int n) {
  int b = (n + BS - 1) / BS;
  if (b < 1) b = 1;
  if (b > MAXB) b = MAXB;
  return b;
}

extern "C" {
void launch_cgm_update_p_dev(int n, int use_pre, int corr, double sty, double psmdyty,
                             double sgnprev, double ztrprev, int iter0, const double* red,
                             const double* r, const double* D, const double* u, const double* s,
                             double* p) {
  hipLaunchKernelGGL(k_update_p_dev, dim3(nblocks(n)), dim3(BS), 0, 0,
                     n, use_pre, corr, sty, psmdyty, sgnprev, ztrprev, iter0, red, r, D, u, s, p);
  HIPCHK();
}
void launch_cgm_dots0(int n, const double* r, const double* D, const double* s, const double* u,
                      double* part, unsigned int* counter, double* out) {
  hipLaunchKernelGGL(k_dots0, dim3(nblocks(n)), dim3(BS), 0, 0, n, r, D, s, u, part, counter, out);
  HIPCHK();
}
// dots1 escreve em out[8..14] (out[0..5] ficam para k_update / k_dots0)
void launch_cgm_dots1(int n, double delta, const double* p, const double* hp, const double* d,
                      const double* g, const double* r,
                      double* part, unsigned int* counter, double* out) {
  hipLaunchKernelGGL(k_dots1, dim3(nblocks(n)), dim3(BS), 0, 0,
                     n, delta, p, hp, d, g, r, part, counter, out + 8);
  HIPCHK();
}
void launch_cgm_update(int n, double coef, double macheps, const double* p, const double* hp,
                       double* d, double* r, const double* xf, const double* lbf, const double* ubf,
                       const int* lindf, const int* uindf,
                       const double* dg, const double* sv, const double* uv,
                       double* part, unsigned int* counter, double* out) {
  hipLaunchKernelGGL(k_update, dim3(nblocks(n)), dim3(BS), 0, 0,
                     n, coef, macheps, p, hp, d, r, xf, lbf, ubf, lindf, uindf, dg, sv, uv,
                     part, counter, out);
  HIPCHK();
}
void launch_scatter(int hne, const double* v, const int* p1, const int* p2, double* val) {
  hipLaunchKernelGGL(k_scatter, dim3((hne + 255) / 256), dim3(256), 0, 0, hne, v, p1, p2, val);
  HIPCHK();
}
}