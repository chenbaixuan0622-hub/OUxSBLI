#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/cast.h>
#include <pybind11/stl.h>
#include <cuda_runtime.h>
#include "fortran_interface.h"


namespace py = pybind11;


template <typename T>
T* get_cuda_pointer(py::object arr) {
  py::tuple data_tuple = arr.attr("__cuda_array_interface__")["data"].cast<py::tuple>();
  uintptr_t ptr = data_tuple[0].cast<uintptr_t>();
  return reinterpret_cast<T*>(ptr);
}


extern "C" void calc_EFG_Euler_wrapper(int nx, int ny, int nz, py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                       py::object Jacobian_obj, py::object QJ_obj, \
                                       py::object E_obj, py::object F_obj, py::object G_obj, \
                                       py::object fx_obj = py::none(), \
                                       py::object fy_obj = py::none(), \
                                       py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *Jacobian = get_cuda_pointer<double>(Jacobian_obj);
  double *QJ       = get_cuda_pointer<double>(QJ_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_EFG_Euler_forcing_c(nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz);
  } else {
    calc_EFG_Euler_c(nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G);
  }
}


extern "C" void calc_EFG_NS_wrapper(int nx, int ny, int nz, py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                    py::object Jacobian_obj, py::object QJ_obj, \
                                    py::object E_obj, py::object F_obj, py::object G_obj, \
                                    py::object fx_obj = py::none(), \
                                    py::object fy_obj = py::none(), \
                                    py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *Jacobian = get_cuda_pointer<double>(Jacobian_obj);
  double *QJ       = get_cuda_pointer<double>(QJ_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_EFG_NS_forcing_c(nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz);
  } else {
    calc_EFG_NS_c(nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G);
  }
}


extern "C" void calc_EFG_LES_wrapper(int nx, int ny, int nz, py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                     py::object Jacobian_obj, py::object QJ_obj, \
                                     py::object E_obj, py::object F_obj, py::object G_obj, \
                                     py::object fx_obj = py::none(), \
                                     py::object fy_obj = py::none(), \
                                     py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *Jacobian = get_cuda_pointer<double>(Jacobian_obj);
  double *QJ       = get_cuda_pointer<double>(QJ_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_EFG_LES_forcing_c(nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz);
  } else {
    calc_EFG_LES_c(nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G);
  }
}


extern "C" void calc_R_wrapper(int nx, int ny, int nz, py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                               py::object E_obj, py::object F_obj, py::object G_obj, py::object R_obj, \
                               py::object fx_obj = py::none(), \
                               py::object fy_obj = py::none(), \
                               py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *R = get_cuda_pointer<double>(R_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_R_forcing_c(nx, ny, nz, dx, dy, dz, E, F, G, fx, fy, fz, R);
  } else {
    calc_R_c(nx, ny, nz, dx, dy, dz, E, F, G, R);
  }
}


extern "C" void calc_step1_wrapper(int nx, int ny, int nz, double coef, py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                    py::object E_obj, py::object F_obj, py::object G_obj, \
                                    py::object Q_obj, py::object Q2_obj, \
                                    py::object fx_obj = py::none(), \
                                    py::object fy_obj = py::none(), \
                                    py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *Q  = get_cuda_pointer<double>(Q_obj);
  double *Q2 = get_cuda_pointer<double>(Q2_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_step1_forcing_c(nx, ny, nz, coef, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2);
  } else {
    calc_step1_c(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2);
  }
}


extern "C" void calc_step_wrapper(int nx, int ny, int nz, double coef1, double coef2, \
                                  py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                  py::object E_obj, py::object F_obj, py::object G_obj, \
                                  py::object Q_obj, py::object Q2_obj, py::object Rs_obj, \
                                  py::object fx_obj = py::none(), \
                                  py::object fy_obj = py::none(), \
                                  py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *Q  = get_cuda_pointer<double>(Q_obj);
  double *Q2 = get_cuda_pointer<double>(Q2_obj);
  double *Rs = get_cuda_pointer<double>(Rs_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_step_forcing_c(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2, Rs);
  } else {
    calc_step_c(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs);
  }
}


extern "C" void calc_step2_3_wrapper(int nx, int ny, int nz, double coef1, double coef2, \
                                      double coef3, double coef4, \
                                      py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                      py::object E_obj, py::object F_obj, py::object G_obj, \
                                      py::object Qin_obj, py::object Qinout_obj, \
                                      py::object fx_obj = py::none(), \
                                      py::object fy_obj = py::none(), \
                                      py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *Qin    = get_cuda_pointer<double>(Qin_obj);
  double *Qinout = get_cuda_pointer<double>(Qinout_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_step2_3_forcing_c(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, fx, fy, fz, Qin, Qinout);
  } else {
    calc_step2_3_c(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Qin, Qinout);
  }
}


extern "C" void calc_step4_wrapper(int nx, int ny, int nz, \
                                    py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                    py::object E_obj, py::object F_obj, py::object G_obj, \
                                    py::object Rs_obj, py::object Q_obj, \
                                    py::object fx_obj = py::none(), \
                                    py::object fy_obj = py::none(), \
                                    py::object fz_obj = py::none()) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *E = get_cuda_pointer<double>(E_obj);
  double *F = get_cuda_pointer<double>(F_obj);
  double *G = get_cuda_pointer<double>(G_obj);
  double *Rs = get_cuda_pointer<double>(Rs_obj);
  double *Q  = get_cuda_pointer<double>(Q_obj);
  if (!fx_obj.is_none() && !fy_obj.is_none() && !fz_obj.is_none()){
    double *fx = get_cuda_pointer<double>(fx_obj);
    double *fy = get_cuda_pointer<double>(fy_obj);
    double *fz = get_cuda_pointer<double>(fz_obj);
    calc_step4_forcing_c(nx, ny, nz, dx, dy, dz, E, F, G, fx, fy, fz, Rs, Q);
  } else {
    calc_step4_c(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q);
  }
}


extern "C" void calc_error_wrapper(int nx, int ny, int nz, \
                                   py::object R1_obj, py::object R2_obj, \
                                   py::object R1_new_obj, py::object R2_new_obj) {
  double *R1     = get_cuda_pointer<double>(R1_obj);
  double *R2     = get_cuda_pointer<double>(R2_obj);
  double *R1_new = get_cuda_pointer<double>(R1_new_obj);
  double *R2_new = get_cuda_pointer<double>(R2_new_obj);
  calc_error_c(nx, ny, nz, R1, R2, R1_new, R2_new);
}


extern "C" void calc_Gauss_step_wrapper(int nx, int ny, int nz, double a1, double a2, \
                                        py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                        py::object R1_obj, py::object R2_obj, \
                                        py::object Q_obj, py::object Q2_obj) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *R1 = get_cuda_pointer<double>(R1_obj);
  double *R2 = get_cuda_pointer<double>(R2_obj);
  double *Q  = get_cuda_pointer<double>(Q_obj);
  double *Q2 = get_cuda_pointer<double>(Q2_obj);
  calc_Gauss_step_c(nx, ny, nz, a1, a2, dx, dy, dz, R1, R2, Q, Q2);
}


extern "C" void calc_Gauss_step_Q_wrapper(int nx, int ny, int nz, double a1, double a2, \
                                          py::object dx_obj, py::object dy_obj, py::object dz_obj, \
                                          py::object R1_obj, py::object R2_obj, \
                                          py::object Q_obj) {
  double *dx = get_cuda_pointer<double>(dx_obj);
  double *dy = get_cuda_pointer<double>(dy_obj);
  double *dz = get_cuda_pointer<double>(dz_obj);
  double *R1 = get_cuda_pointer<double>(R1_obj);
  double *R2 = get_cuda_pointer<double>(R2_obj);
  double *Q  = get_cuda_pointer<double>(Q_obj);
  calc_Gauss_step_Q_c(nx, ny, nz, a1, a2, dx, dy, dz, R1, R2, Q);
}


PYBIND11_MODULE(cufd, m) {
  m.def("calc_EFG_Euler", &calc_EFG_Euler_wrapper, 
        py::arg("nx"), py::arg("ny"), py::arg("nz"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("Jacobian"), py::arg("QJ"),
        py::arg("E"),  py::arg("F"),  py::arg("G"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_EFG_NS", &calc_EFG_NS_wrapper,
        py::arg("nx"), py::arg("ny"), py::arg("nz"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("Jacobian"), py::arg("QJ"),
        py::arg("E"),  py::arg("F"),  py::arg("G"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_EFG_LES", &calc_EFG_LES_wrapper,
        py::arg("nx"), py::arg("ny"), py::arg("nz"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("Jacobian"), py::arg("QJ"),
        py::arg("E"),  py::arg("F"),  py::arg("G"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_R", &calc_R_wrapper,
        py::arg("nx"), py::arg("ny"), py::arg("nz"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("E"),  py::arg("F"),  py::arg("G"), py::arg("R"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_step1", &calc_step1_wrapper,
        py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("coef"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("E"),  py::arg("F"),  py::arg("G"), py::arg("Q"), py::arg("Q2"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_step", &calc_step_wrapper,
        py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("coef1"), py::arg("coef2"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("E"),  py::arg("F"),  py::arg("G"),
        py::arg("Q"),  py::arg("Q2"), py::arg("Rs"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_step2_3", &calc_step2_3_wrapper,
        py::arg("nx"), py::arg("ny"), py::arg("nz"),
        py::arg("coef1"), py::arg("coef2"), py::arg("coef3"), py::arg("coef4"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("E"),  py::arg("F"),  py::arg("G"),
        py::arg("Qin"),  py::arg("Qinout"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_step4", &calc_step4_wrapper,
        py::arg("nx"), py::arg("ny"), py::arg("nz"),
        py::arg("dx"), py::arg("dy"), py::arg("dz"),
        py::arg("E"),  py::arg("F"),  py::arg("G"),
        py::arg("Rs"),  py::arg("Q"),
        py::arg("fx") = py::none(),
        py::arg("fy") = py::none(),
        py::arg("fz") = py::none());
  m.def("calc_error", &calc_error_wrapper);
  m.def("calc_Gauss_step", &calc_Gauss_step_wrapper);
  m.def("calc_Gauss_step_Q", &calc_Gauss_step_Q_wrapper);
}

