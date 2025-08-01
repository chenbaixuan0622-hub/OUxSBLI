// rk.cu
#include <cuda_runtime.h>
#include <cstdio>
#include "array.h"


__device__ __forceinline__ double R(Array4D<const double> e, Array4D<const double> f, Array4D<const double> g,
                                    const double dx, const double dy, const double dz, int i, int j, int k, int itr) {
  return (-e(k-1,j-1,i-1,itr) + e(k-1,j-1,i,itr)) / dx + 
         (-f(k-1,j-1,i-1,itr) + f(k-1,j,i-1,itr)) / dy + 
         (-g(k-1,j-1,i-1,itr) + g(k,j-1,i-1,itr)) / dz;
}


//3rd-order TVD Runge-kutta


__global__ void RK33_1_kernel(Array4D<const double> q, Array4D<const double> e, Array4D<const double> f, Array4D<const double> g, Array4D<double> q1,
                              const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1) {
    for (int itr = 0; itr < 5; itr++) {
      q1(k,j,i,itr) = q(k,j,i,itr) - dt * R(e, f, g, dx, dy, dz, i, j, k, itr);
    }
  }
}


void RK33_1(const double* q, const double* e, const double* f, const double* g, double* q1, 
            const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q(q, nz, ny, nx, 5);
  Array4D<const double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<const double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<const double> G(g, nz-1, ny-2, nx-2, 5);
  Array4D<double> Q1(q1, nz, ny, nx, 5);
  RK33_1_kernel<<<grid, block>>>(Q, E, F, G, Q1, dt, dx, dy, dz, nx, ny, nz);
}


__global__ void RK33_2_kernel(Array4D<const double> q, Array4D<const double> e, Array4D<const double> f, Array4D<const double> g, Array4D<double> q1, 
                              const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1) {
    for (int itr = 0; itr < 5; itr++) {
      q1(k,j,i,itr) = 0.25e0 * (3.e0 * q(k,j,i,itr) + q1(k,j,i,itr) - dt * R(e, f, g, dx, dy, dz, i, j, k, itr));
    }
  }
}


void RK33_2(const double* q, const double* e, const double* f, const double* g, double* q1, 
            const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q(q, nz, ny, nx, 5);
  Array4D<const double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<const double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<const double> G(g, nz-1, ny-2, nx-2, 5);
  Array4D<double> Q1(q1, nz, ny, nx, 5);
  RK33_2_kernel<<<grid, block>>>(Q, E, F, G, Q1, dt, dx, dy, dz, nx, ny, nz);
}


__global__ void RK33_3_kernel(Array4D<const double> q1, Array4D<const double> e, Array4D<const double> f, Array4D<const double> g, Array4D<double> q, 
                              const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1) {
    for (int itr = 0; itr < 5; itr++) {
      q(k,j,i,itr) = (q(k,j,i,itr) + 2.e0 * q1(k,j,i,itr) - 2.e0 * dt * R(e, f, g, dx, dy, dz, i, j, k, itr)) / 3.e0;
    }
  }
}


void RK33_3(const double* q1, const double* e, const double* f, const double* g, double* q, 
            const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q1(q1, nz, ny, nx, 5);
  Array4D<const double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<const double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<const double> G(g, nz-1, ny-2, nx-2, 5);
  Array4D<double> Q(q, nz, ny, nx, 5);
  RK33_3_kernel<<<grid, block>>>(Q1, E, F, G, Q, dt, dx, dy, dz, nx, ny, nz);
}


//4th-order Runge-kutta


__global__ void RK44_1_3_kernel(Array4D<const double> q, Array4D<const double> e, Array4D<const double> f, Array4D<const double> g, Array4D<double> q1,
                              Array4D<double> r1, const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  double r_tmp;
  if (0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1) {
    for (int itr = 0; itr < 5; itr++) {
      r_tmp = R(e, f, g, dx, dy, dz, i, j, k, itr);
      q1(k,j,i,itr)  = q(k,j,i,itr) - dt * r_tmp;
      r1(k,j,i,itr) += r_tmp;
    }
  }
}


__global__ void RK44_4_kernel(Array4D<const double> q, Array4D<const double> e, Array4D<const double> f, Array4D<const double> g, Array4D<double> q1,
                              Array4D<double> r1, const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1) {
    for (int itr = 0; itr < 5; itr++) {
      r1(k,j,i,itr) += R(e, f, g, dx, dy, dz, i, j, k, itr);
      q1(k,j,i,itr) = q(k,j,i,itr) - dt * r1(k,j,i,itr) / 6.e0;
      r1(k,j,i,itr) = 0.e0;
    }
  }
}


void RK44_1(const double* q, const double* e, const double* f, const double* g, double* q1, double* r1,
            const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q(q, nz, ny, nx, 5);
  Array4D<const double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<const double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<const double> G(g, nz-1, ny-2, nx-2, 5);
  Array4D<double> Q1(q1, nz, ny, nx, 5);
  Array4D<double> R1(r1, nz-2, ny-2, nx-2, 5);
  RK44_1_3_kernel<<<grid, block>>>(Q, E, F, G, Q1, R1, dt, dx, dy, dz, nx, ny, nz);
}


void RK44_2(const double* q, const double* e, const double* f, const double* g, double* q1, double* r1,
            const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q(q, nz, ny, nx, 5);
  Array4D<const double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<const double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<const double> G(g, nz-1, ny-2, nx-2, 5);
  Array4D<double> Q1(q1, nz, ny, nx, 5);
  Array4D<double> R1(r1, nz-2, ny-2, nx-2, 5);
  RK44_1_3_kernel<<<grid, block>>>(Q, E, F, G, Q1, R1, dt, dx, dy, dz, nx, ny, nz);
}


void RK44_3(const double* q, const double* e, const double* f, const double* g, double* q1, double* r1,
            const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q(q, nz, ny, nx, 5);
  Array4D<const double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<const double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<const double> G(g, nz-1, ny-2, nx-2, 5);
  Array4D<double> Q1(q1, nz, ny, nx, 5);
  Array4D<double> R1(r1, nz-2, ny-2, nx-2, 5);
  RK44_1_3_kernel<<<grid, block>>>(Q, E, F, G, Q1, R1, dt, dx, dy, dz, nx, ny, nz);
}


void RK44_4(const double* q1, const double* e, const double* f, const double* g, double* q, double* r1,
            const double dt, const double dx, const double dy, const double dz, const int nx, const int ny, const int nz) {
  dim3 block(8, 8, 8);
  dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y, (nz + block.z - 1) / block.z);
  Array4D<const double> Q1(q1, nz, ny, nx, 5);
  Array4D<const double> E(e, nz-2, ny-2, nx-1, 5);
  Array4D<const double> F(f, nz-2, ny-1, nx-2, 5);
  Array4D<const double> G(g, nz-1, ny-2, nx-2, 5);
  Array4D<double> Q(q, nz, ny, nx, 5);
  Array4D<double> R1(r1, nz-2, ny-2, nx-2, 5);
  RK44_4_kernel<<<grid, block>>>(Q1, E, F, G, Q, R1, dt, dx, dy, dz, nx, ny, nz);
}

