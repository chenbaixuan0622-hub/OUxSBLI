#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/cast.h>
#include <pybind11/stl.h>
#include <cuda_runtime.h>


extern void calc_flux(const double* q, double* e, double* f, double* g, const double c, 
                      const int nx, const int ny, const int nz);
// 3rd-order TVD Runge-kutta
extern void RK33_1(const double* q, const double* e, const double* f, const double* g, 
                   double* q1, double dt, double dx, double dy, double dz, int nx, int ny, int nz);
extern void RK33_2(const double* q, const double* e, const double* f, const double* g, 
                   double* q1, double dt, double dx, double dy, double dz, int nx, int ny, int nz);
extern void RK33_3(const double* q1, const double*e, const double* f, const double* g, 
                   double* q,  double dt, double dx, double dy, double dz, int nx, int ny, int nz);
// 4th-order Tunge-kutta
extern void RK44_1(const double* q, const double* e, const double* f, const double* g, 
                   double* q1, double* r1, double dt, double dx, double dy, double dz, int nx, int ny, int nz);
extern void RK44_2(const double* q, const double* e, const double* f, const double* g, 
                   double* q1, double* r1, double dt, double dx, double dy, double dz, int nx, int ny, int nz);
extern void RK44_3(const double* q, const double* e, const double* f, const double* g, 
                   double* q1, double* r1, double dt, double dx, double dy, double dz, int nx, int ny, int nz);
extern void RK44_4(const double* q, const double* e, const double* f, const double* g, 
                   double* q1, double* r1, double dt, double dx, double dy, double dz, int nx, int ny, int nz);


namespace py = pybind11;


template <typename T>
T* get_cuda_pointer(py::object arr) {
    py::tuple data_tuple = arr.attr("__cuda_array_interface__")["data"].cast<py::tuple>();
    uintptr_t ptr = data_tuple[0].cast<uintptr_t>();
    return reinterpret_cast<T*>(ptr);
}


extern "C" void calc_flux_wrapper(py::object q_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                                  double c, int nx, int ny, int nz) {
  double* q_ptr = get_cuda_pointer<double>(q_obj);
  double* e_ptr = get_cuda_pointer<double>(e_obj);
  double* f_ptr = get_cuda_pointer<double>(f_obj);
  double* g_ptr = get_cuda_pointer<double>(g_obj);
  calc_flux(q_ptr, e_ptr, f_ptr, g_ptr, c, nx, ny, nz);
}


// 3rd-order TVD Runge-kutta


extern "C" void RK33_1_wrapper(py::object q_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                               py::object q1_obj, double dt, double dx, double dy, double dz, 
                               int nx, int ny, int nz) {
  double* q_ptr  = get_cuda_pointer<double>(q_obj);
  double* e_ptr  = get_cuda_pointer<double>(e_obj);
  double* f_ptr  = get_cuda_pointer<double>(f_obj);
  double* g_ptr  = get_cuda_pointer<double>(g_obj);
  double* q1_ptr = get_cuda_pointer<double>(q1_obj);
  RK33_1(q_ptr, e_ptr, f_ptr, g_ptr, q1_ptr, dt, dx, dy, dz, nx, ny, nz);
}


extern "C" void RK33_2_wrapper(py::object q_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                               py::object q1_obj, double dt, double dx, double dy, double dz, 
                               int nx, int ny, int nz) {
  double* q_ptr  = get_cuda_pointer<double>(q_obj);
  double* e_ptr  = get_cuda_pointer<double>(e_obj);
  double* f_ptr  = get_cuda_pointer<double>(f_obj);
  double* g_ptr  = get_cuda_pointer<double>(g_obj);
  double* q1_ptr = get_cuda_pointer<double>(q1_obj);
  RK33_2(q_ptr, e_ptr, f_ptr, g_ptr, q1_ptr, dt, dx, dy, dz, nx, ny, nz);
}


extern "C" void RK33_3_wrapper(py::object q1_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                               py::object q_obj, double dt, double dx, double dy, double dz, 
                               int nx, int ny, int nz) {
  double* q_ptr  = get_cuda_pointer<double>(q_obj);
  double* e_ptr  = get_cuda_pointer<double>(e_obj);
  double* f_ptr  = get_cuda_pointer<double>(f_obj);
  double* g_ptr  = get_cuda_pointer<double>(g_obj);
  double* q1_ptr = get_cuda_pointer<double>(q1_obj);
  RK33_3(q1_ptr, e_ptr, f_ptr, g_ptr, q_ptr, dt, dx, dy, dz, nx, ny, nz);
}


// 4th-order Tunge-kutta


extern "C" void RK44_1_wrapper(py::object q_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                               py::object q1_obj, py::object r1_obj, double dt, double dx, double dy, double dz, 
                               int nx, int ny, int nz) {
  double* q_ptr  = get_cuda_pointer<double>(q_obj);
  double* e_ptr  = get_cuda_pointer<double>(e_obj);
  double* f_ptr  = get_cuda_pointer<double>(f_obj);
  double* g_ptr  = get_cuda_pointer<double>(g_obj);
  double* q1_ptr = get_cuda_pointer<double>(q1_obj);
  double* r1_ptr = get_cuda_pointer<double>(r1_obj);
  RK44_1(q_ptr, e_ptr, f_ptr, g_ptr, q1_ptr, r1_ptr, dt, dx, dy, dz, nx, ny, nz);
}


extern "C" void RK44_2_wrapper(py::object q_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                               py::object q1_obj, py::object r1_obj, double dt, double dx, double dy, double dz, 
                               int nx, int ny, int nz) {
  double* q_ptr  = get_cuda_pointer<double>(q_obj);
  double* e_ptr  = get_cuda_pointer<double>(e_obj);
  double* f_ptr  = get_cuda_pointer<double>(f_obj);
  double* g_ptr  = get_cuda_pointer<double>(g_obj);
  double* q1_ptr = get_cuda_pointer<double>(q1_obj);
  double* r1_ptr = get_cuda_pointer<double>(r1_obj);
  RK44_2(q_ptr, e_ptr, f_ptr, g_ptr, q1_ptr, r1_ptr, dt, dx, dy, dz, nx, ny, nz);
}


extern "C" void RK44_3_wrapper(py::object q_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                               py::object q1_obj, py::object r1_obj, double dt, double dx, double dy, double dz, 
                               int nx, int ny, int nz) {
  double* q_ptr  = get_cuda_pointer<double>(q_obj);
  double* e_ptr  = get_cuda_pointer<double>(e_obj);
  double* f_ptr  = get_cuda_pointer<double>(f_obj);
  double* g_ptr  = get_cuda_pointer<double>(g_obj);
  double* q1_ptr = get_cuda_pointer<double>(q1_obj);
  double* r1_ptr = get_cuda_pointer<double>(r1_obj);
  RK44_3(q_ptr, e_ptr, f_ptr, g_ptr, q1_ptr, r1_ptr, dt, dx, dy, dz, nx, ny, nz);
}


extern "C" void RK44_4_wrapper(py::object q_obj, py::object e_obj, py::object f_obj, py::object g_obj, 
                               py::object q1_obj, py::object r1_obj, double dt, double dx, double dy, double dz, 
                               int nx, int ny, int nz) {
  double* q_ptr  = get_cuda_pointer<double>(q_obj);
  double* e_ptr  = get_cuda_pointer<double>(e_obj);
  double* f_ptr  = get_cuda_pointer<double>(f_obj);
  double* g_ptr  = get_cuda_pointer<double>(g_obj);
  double* q1_ptr = get_cuda_pointer<double>(q1_obj);
  double* r1_ptr = get_cuda_pointer<double>(r1_obj);
  RK44_4(q_ptr, e_ptr, f_ptr, g_ptr, q1_ptr, r1_ptr, dt, dx, dy, dz, nx, ny, nz);
}


PYBIND11_MODULE(cufd, m) {
  m.def("calc_flux", &calc_flux_wrapper);
  m.def("RK33_1", &RK33_1_wrapper);
  m.def("RK33_2", &RK33_2_wrapper);
  m.def("RK33_3", &RK33_3_wrapper);
  m.def("RK44_1", &RK44_1_wrapper);
  m.def("RK44_2", &RK44_2_wrapper);
  m.def("RK44_3", &RK44_3_wrapper);
  m.def("RK44_4", &RK44_4_wrapper);
}

