// flux.cu
#include <cuda_runtime.h>
#include <cstdio>
#include "array.h"


__device__ __forceinline__ void keep(Array2D<double> Q, const double* N, const double gamma, const int i, Array1DView<double> F) {
  double rho1 = Q(0,0), rho2 = Q(1,0);
  double V1[3] = {Q(0,1), Q(0,2), Q(0,3)};
  double V2[3] = {Q(1,1), Q(1,2), Q(1,3)};
  double p1 = Q(0,4), p2 = Q(1,4);
  F[0] = 0.25e0 * (rho1 + rho2) * (V1[i] + V2[i]);
  F[1] = F[0] * 0.5e0 * (V1[0] + V2[0]) + 0.5e0 * (p1 + p2) * N[0];
  F[2] = F[0] * 0.5e0 * (V1[1] + V2[1]) + 0.5e0 * (p1 + p2) * N[1];
  F[3] = F[0] * 0.5e0 * (V1[2] + V2[2]) + 0.5e0 * (p1 + p2) * N[2];
  F[4] = F[0] * 0.5e0 * (p1 / rho1 + p2 / rho2) / (gamma - 1.e0) 
         + 0.5e0 * (V1[i] * p2 + V2[i] * p1)
         + F[0] * 0.5e0 * (V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]);
}


__device__ __forceinline__ void set_q(Array1DView<const double> Q, Array1DView<double> q, const double gamma) {
  q[0] = Q[0];
  q[1] = Q[1] / Q[0];
  q[2] = Q[2] / Q[0];
  q[3] = Q[3] / Q[0];
  q[4] = (gamma - 1.e0) * (Q[4] - 0.5e0 * (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) / Q[0]);
}


__global__ void calc_e_kernel(Array4D<const double> Q, Array4D<double> E, 
                              const double gamma, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  double q_d[2*5];
  double N[3] = {1.e0, 0.e0, 0.e0};
  if (0 <= i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1) {
    Array2D<double> q(q_d, 2, 5);
    for (int ac = 0; ac < 2; ac++) {
      set_q(Q.slice(k,j,i+ac), q.slice(ac), gamma);
    }
    keep(q, N, gamma, 0, E.slice(k-1,j-1,i));
  }
}


__global__ void calc_f_kernel(Array4D<const double> Q, Array4D<double> F, 
                              const double gamma, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  double q_d[2*5];
  double N[3] = {0.e0, 1.e0, 0.e0};
  if (0 < i && i < nx-1 && 0 <= j && j < ny-1 && 0 < k && k < nz-1) {
    Array2D<double> q(q_d, 2, 5);
    for (int ac = 0; ac < 2; ac++) {
      set_q(Q.slice(k,j+ac,i), q.slice(ac), gamma);
    }
    keep(q, N, gamma, 1, F.slice(k-1,j,i-1));
  }
}


__global__ void calc_g_kernel(Array4D<const double> Q, Array4D<double> G, 
                              const double gamma, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  double q_d[2*5];
  double N[3] = {0.e0, 0.e0, 1.e0};
  if (0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 <= k && k < nz-1) {
    Array2D<double> q(q_d, 2, 5);
    for (int ac = 0; ac < 2; ac++) {
      set_q(Q.slice(k+ac,j,i), q.slice(ac), gamma);
    }
    keep(q, N, gamma, 2, G.slice(k,j-1,i-1));
  }
}


__device__ __forceinline__ double mu(double T) {
  double mu0 = 1.716e-5;
  double T0  = 273.2e0;
  double S   = 111.e0;
  return mu0 * ((T0 + S) / (T + S)) * (T / T0) * sqrt(T / T0);
}


__device__ __forceinline__ double mu2(double* T) {
  return 0.5e0 * (mu(T[0]) + mu(T[1]));
}


__device__ __forceinline__ double* mu23(Array2D<double> T) {
  double m[2];
  m[0] = 0.25e0 * (mu(T(0,0)) + mu(T(0,1)) + mu(T(1,0)) + mu(T(1,1)));
  m[1] = 0.25e0 * (mu(T(0,1)) + mu(T(0,2)) + mu(T(1,1)) + mu(T(1,2)));
  return m;
}


__device__ __forceinline__ double* mu32(Array2D<double> T) {
  double m[2];
  m[0] = 0.25e0 * (mu(T(0,0)) + mu(T(1,0)) + mu(T(0,1)) + mu(T(1,1)));
  m[1] = 0.25e0 * (mu(T(1,0)) + mu(T(2,0)) + mu(T(1,1)) + mu(T(2,1)));
  return m;
}


__device__ __forceinline__ double dy23(double* m, Array2D<double> a, const double d) {
  return 0.25e0 * (m[0] * (-a(0,0) + a(0,1) - a(1,0) + a(1,1)) 
                 + m[1] * (-a(0,1) + a(0,2) - a(1,1) + a(1,2)));
}


__device__ __forceinline__ double dy32(double* m, Array2D<double> a, const double d) {
  return 0.25e0 * (m[0] * (-a(0,0) + a(1,0) - a(0,1) + a(1,1)) 
                 + m[1] * (-a(1,0) + a(2,0) - a(1,1) + a(2,1)));
}


__global__ void calc_ev_kernel(Array4D<const double> Q, Array4D<double> E, const double gamma, 
                               const double Rgas, const double Pr, const int nx, const int ny, const int nz, 
                               const double dx, const double dy, const double dz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  double Cp = gamma * Rgas / (gamma - 1.e0);
  double uy[], uz[], v[], w[], mx, my[], mz[];
  double mx, mux, mvy, mwz, muy, mvx, mwx, muz, txx, txy, txz, utxx, vtxy, wtxz, kTx;
  if () {
    uy(0,0), uy(1,0), uy(2,0) = Q(k,j-1,i,1)   / Q(k,j-1,i,0),   Q(k,j,i,1)   / Q(k,j,i,0),   Q(k,j+1,i,1)   / Q(k,j+1,i,0);
    uy(0,1), uy(1,1), uy(2,1) = Q(k,j-1,i+1,1) / Q(k,j-1,i+1,0), Q(k,j,i+1,1) / Q(k,j,i+1,0), Q(k,j+1,i+1,1) / Q(k,j+1,i+1,0);
    uz(0,0), uz(1,0), uz(2,0) = Q(k-1,j,i,1)   / Q(k-1,j,i,0),   Q(k,j,i,1)   / Q(k,j,i,0),   Q(k+1,j,i,1)   / Q(k+1,j,i,0);
    uz(0,1), uz(1,1), uz(2,1) = Q(k-1,j,i+1,1) / Q(k-1,j,i+1,0), Q(k,j,i+1,1) / Q(k,j,i+1,0). Q(k+1,j,i+1,1) / Q(k+1,j,i+1,0);
    v(0,0), v(1,0), v(2,0)    = Q(k,j-1,i,2)   / Q(k,j-1,i,0),   Q(k,j,i,2)   / Q(k,j,i,0),   Q(k,j+1,i,2)   / Q(k,j+1,i,0);
    v(0,1), v(1,1), v(2,1)    = Q(k,j-1,i+1,2) / Q(k,j-1,i+1,0), Q(k,j,i+1,2) / Q(k,j,i+1,0), Q(k,j+1,i+1,2) / Q(k,j+1,i+1,0);
    w(0,0), w(1,0), w(2,0)    = Q(k-1,j,i,3)   / Q(k-1,j,i,0),   Q(k,j,i,3)   / Q(k,j,i,0),   Q(k+1,j,i,3)   / Q(k+1,j,i,0);
    w(0,1), w(1,1), w(2,1)    = Q(k-1,j,i+1,3) / Q(k-1,j,i+1,0), Q(k,j,i+1,3) / Q(k,j,i+1,0). Q(k+1,j,i+1,3) / Q(k+1,j,i+1,0);
    mx   = mu2();
    my   = mu32();
    mz   = mu32();
    mux  = mx * (-uy(1,0) + uy(1,1)) * dx;
    mvx  = mx * ( -v(1,0) +  v(1,1)) * dx;
    mwx  = mx * ( -w(1,0) +  w(1,1)) * dx;
    muy  = dy32(my, uy, dy);
    mvy  = dy32(my,  v, dy);
    muz  = dy32(mz, uz, dz);
    mwz  = dy32(mz,  w, dz);
    txx  = 2.e0 * (2.e0 * mux - mvy - mwz) / 3.e0;
    txy  = muy + mvx;
    txz  = mwx + muz;
    utxx = 0.5e0 * (uy(1,0) + uy(1,1)) * txx;
    vtxy = 0.5e0 * ( v(1,0) +  v(1,1)) * txy;
    wtxz = 0.5e0 * ( w(1,0) +  w(1,1)) * txz;
    kTx  = Cp * mx * () * dx / Pr;
    E(k-1,j-1,i,1) -= txx;
    E(k-1,j-1,i,2) -= txy;
    E(k-1,j-1,i,3) -= txz;
    E(k-1,j-1,i,4) -= (utxx + vtxy + wtxz + kTx);
  }
}


__global__ void calc_fv_kernel(Array4D<const double> Q, Array4D<double> F, const double gamma, 
                               const double Rgas, const double Pr, const int nx, const int ny, const int nz,
                               const double dx, const double dy, const double dz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
  double Cp = gamma * Rgas / (gamma - 1.e0);
  double vx[], vz[], u[], w[], mx[], my, mz[];
  double my, muy, mvx, mvy, mwz, mux, mvz, mwy, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy;
  if () {
    vx(0,0), vx(0,1), vx(0,2) = Q(k,j,i-1,2)   / Q(k,j,i-1,0),   Q(k,j,i,2)   / Q(k,j,i,0),   Q(k,j,i+1,2)   / Q(k,j,i+1,0);
    vx(1,0), vx(1,1), vx(1,2) = Q(k,j+1,i-1,2) / Q(k,j+1,i-1,0), Q(k,j+1,i,2) / Q(k,j+1,i,0), Q(k,j+1,i+1,2) / Q(k,j+1,i+1,0);
    vz(0,0), vz(1,0), vz(2,0) = Q(k-1,j,i,2)   / Q(k-1,j,i,0),   Q(k,j,i,2)   / Q(k,j,i,0),   Q(k+1,j,i,2)   / Q(k+1,j,i,0);
    vz(0,1), vz(1,1), vz(2,1) = Q(k-1,j+1,i,2) / Q(k-1,j+1,i,0), Q(k,j+1,i,2) / Q(k,j+1,i,0), Q(k+1,j+1,i,2) / Q(k+1,j+1,i,0);
    u(0,0), u(0,1), u(0,2)    = Q(k,j,i-1,1)   / Q(k,j,i-1,0),   Q(k,j,i,1)   / Q(k,j,i,0),   Q(k,j,i+1,1)   / Q(k,j,i+1,0);
    u(1,0), u(1,1), u(1,2)    = Q(k,j+1,i-1,1) / Q(k,j+1,i-1,0), Q(k,j+1,i,1) / Q(k,j+1,i,0), Q(k,j+1,i+1,1) / Q(k,j+1,i+1,0);
    w(0,0), w(1,0), w(2,0)    = Q(k-1,j,i,3)   / Q(k-1,j,i,0),   Q(k,j,i,3)   / Q(k,j,i,0),   Q(k+1,j,i,3)   / Q(k+1,j,i,0);
    w(0,1), w(1,1), w(2,1)    = Q(k-1,j+1,i,3) / Q(k-1,j+1,i,0), Q(k,j+1,i,3) / Q(k,j+1,i,0), Q(k+1,j+1,i,3) / Q(k+1,j+1,i,0);
    mx   = mu23();
    my   = mu2();
    mz   = mu32();
    muy  = my * ( -u(0,1) +  u(1,1)) * dy;
    mvy  = my * (-vx(0,1) + vx(1,1)) * dy;
    mwy  = my * ( -w(1,0) +  w(1,1)) * dy;
    mux  = dy23(mx,  u, dx);
    mvx  = dy23(mx, vx, dx);
    mvw  = dy32(mz, vz, dz);
    mww  = dy32(mz,  w, dz);
    tyx  = muy + mvx;
    tyy  = 2.e0 * (2.e0 * mvy - mwz - mux) / 3.e0;
    tyz  = mvz + mwy;
    utyx = 0.5e0 * ( u(0,1) +  u(1,1)) * tyx;
    vtyy = 0.5e0 * (vx(0,1) + vx(1,1)) * tyy;
    wtyz = 0.5e0 * ( w(1,0) +  w(1,1)) * tyz;
    kTy  = Cp * my * () * dy / Pr;;
    F(k-1,j,i-1,1) -= tyx;
    F(k-1,j,i-1,2) -= tyy;
    F(k-1,j,i-1,3) -= tyz;
    F(k-1,j,i-1,4) -= (utyx + vtyy + wtyz + kTy);
  }
}


__global__ void calc_gv_kernel(Array4D<const double> Q, Array4D<double> G, const double gamma, 
                               const double Rgas, const double Pr, const int nx, const int ny, const int nz, 
                               const double dx, const double dy, const double dz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  double Cp = gamma * Rgas / (gamma - 1.e0);
  double wx[], wy[], u[], v[], mx[], my[], mz;
  double mz, mwx, muz, mvz, mwy, mwz, mux, mvy, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz;
  if () {
    wx(0,0), wx(0,1), wx(0,2) = Q(k,j,i-1,3)   / Q(k,j,i-1,0),   Q(k,j,i,3)   / Q(k,j,i,0),   Q(k,j,i+1,3)   / Q(k,j,i+1,0);
    wx(1,0), wx(1,1), wx(1,2) = Q(k+1,j,i-1,3) / Q(k+1,j,i-1,0), Q(k+1,j,i,3) / Q(k+1,j,i,0), Q(k+1,j,i+1,3) / Q(k+1,j,i+1,0);
    wy(0,0), wy(0,1), wy(0,2) = Q(k,j-1,i,3)   / Q(k,j-1,i,0),   Q(k,j,i,3)   / Q(k,j,i,0),   Q(k,j+1,i,3)   / Q(k,j+1,i,0);
    wy(1,0), wy(1,1), wy(1,2) = Q(k+1,j-1,i,3) / Q(k+1,j-1,i,0), Q(k+1,j,i,3) / Q(k+1,j,i,0), Q(k+1,j+1,i,3) / Q(k+1,j+1,i,0);
    u(0,0), u(0,1), u(0,2)    = Q(k,j,i-1,1)   / Q(k,j,i-1,0),   Q(k,j,i,1)   / Q(k,j,i,0),   Q(k,j,i+1,1)   / Q(k,j,i+1,0);
    u(1,0), u(1,1), u(1,2)    = Q(k+1,j,i-1,1) / Q(k+1,j,i-1,0), Q(k+1,j,i,1) / Q(k+1,j,i,0), Q(k+1,j,i+1,1) / Q(k+1,j,i+1,0);
    v(0,0), v(0,1), v(0,2)    = Q(k,j-1,i,2)   / Q(k,j-1,i,0),   Q(k,j,i,2)   / Q(k,j,i,0),   Q(k,j+1,i,2)   / Q(k,j+1,i,0);
    v(1,0), w(1,1), w(1,2)    = Q(k+1,j-1,i,2) / Q(k+1,j-1,i,0), Q(k+1,j,i,2) / Q(k+1,j,i,0), Q(k+1,j+1,i,2) / Q(k+1,j+1,i,0);
    mx   = mu23();
    my   = mu23();
    mz   = mu2();
    muz  = mz * ( -u(0,1) +  u(1,1)) * dz;
    mvz  = mz * ( -v(0,1) +  v(1,1)) * dz;
    mwz  = mz * (-wx(0,1) + wx(1,1)) * dz;
    mux  = dy23(mx,  u, dx);
    mwx  = dy23(mx, wx, dx);
    mvy  = dy23(my,  v, dy);
    mwy  = dy23(my, wy, dy);
    tzx  = mwx + muz;
    tzy  = mvz + mwy;
    tzz  = 2.e0 * (2.e0 * mwz - mux - mvy) / 3.e0;
    utzx = 0.5e0 * ( u(0,1) +  u(1,1)) * tzx;
    vtzy = 0.5e0 * ( v(0,1) +  v(1,1)) * tzy;
    wtzz = 0.5e0 * (wx(0,1) + wx(1,1)) * tzz;
    kTz  = Cp * mz * () * dz / Pr;
    G(k,j-1,i-1,1) -= tzx;
    G(k,j-1,i-1,2) -= tzy;
    G(k,j-1,i-1,3) -= tzz;
    G(k,j-1,i-1,4) -= (utzx + vtzy + wtzz + kTz);
  }
}


void calc_flux(const double* q, double* e, double* f, double* g, 
               const double gamma, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q(q, nz, ny, nx, 5);
  Array4D<double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<double> G(g, nz-1, ny-2, nx-2, 5);
  calc_e_kernel<<<grid, block>>>(Q, E, gamma, nx, ny, nz);
  calc_f_kernel<<<grid, block>>>(Q, F, gamma, nx, ny, nz);
  calc_g_kernel<<<grid, block>>>(Q, G, gamma, nx, ny, nz);
  calc_ev_kernel<<<grid, block>>>(Q, E, gamma, Rgas, nx, ny, nz);
  calc_fv_kernel<<<grid, block>>>(Q, F, gamma, Rgas, nx, ny, nz);
  calc_gv_kernel<<<grid, block>>>(Q, G, gamma, Rgas, nx, ny, nz);
}

