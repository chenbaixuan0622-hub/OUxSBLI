//array.d
#pragma once


template <typename T>
struct Array1DView {
  T* data;
  int size;
  int stride;

  __device__ __host__
  Array1DView(T* d, int s, int st)
    : data(d), size(s), stride(st) {}

  __device__ __host__ __forceinline__
  T& operator[](int i) {
    return data[i * stride];
  }
};


template <typename T>
struct Array2D {
  T* data;
  int nx, ny;
  __device__ __host__
  Array2D(T* d, int y, int x)
    : data(d), ny(y), nx(x) {}

  __device__ __host__ __forceinline__
  T& operator()(int j, int i) {
    return data[i + nx * j];
  }

  __device__ __host__ __forceinline__
  Array1DView<T> slice(int j) {
    return Array1DView<T>(&data[nx * j], nx, 1);
  }
};


template <typename T>
struct Array3D {
  T* data;
  int nx, ny, nz;
  __device__ __host__
  Array3D(T* d, int z, int y, int x)
    : data(d), nz(z), ny(y), nx(x) {}

  __device__ __host__ __forceinline__
  T& operator()(int k, int j, int i) {
    return data[i + nx * (j + ny * k)];
  }

  __device__ __host__ __forceinline__
  Array1DView<T> slice(int k, int j) {
    return Array1DView<T>(&data[nx * (j + ny * k)], nx, 1);
  }
};


template <typename T>
struct Array4D {
  T* data;
  int nw, nx, ny, nz;
  __device__ __host__
  Array4D(T* d, int z, int y, int x, int w)
    : data(d), nz(z), ny(y), nx(x), nw(w) {}

  __device__ __host__ __forceinline__
  T& operator()(int k, int j, int i, int h) {
    return data[h + nw * (i + nx * (j + ny * k))];
  }

  __device__ __host__ __forceinline__
  Array1DView<T> slice(int k, int j, int i) {
    return Array1DView<T>(&data[nw * (i + nx * (j + ny * k))], nw, 1);
  }
};

