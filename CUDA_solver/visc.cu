// visc.cu
#include <cuda_runtime.h>
#include <cstdio>
#include "array.h"


__device__ __forceinline__ void Sutherland(T) {
  double mu0 = 1.716e-5;
  double T0  = 273.2e0;
  double S   = 111.e0;
  return mu0 * ((T0 + S) / (T + S)) * (T / T0) * sqrt(T / T0);
}

