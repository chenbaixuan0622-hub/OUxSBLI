#include <stdio.h>
#include <cuda_runtime.h>


extern "C" void calc_EFG_Euler_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                                 double* Jacobian, double* QJ, double* E, double* F, double* G);
extern "C" void calc_EFG_Euler_forcing_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                                         double* Jacobian, double* QJ, double* E, double* F, double* G, \
                                         double* fx, double* fy, double* fz);
extern "C" void calc_EFG_NS_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                              double* Jacobian, double* QJ, double* E, double* F, double* G);
extern "C" void calc_EFG_NS_forcing_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                                      double* Jacobian, double* QJ, double* E, double* F, double* G, \
                                      double* fx, double* fy, double* fz);
extern "C" void calc_EFG_LES_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                               double* Jacobian, double* QJ, double* E, double* F, double* G);
extern "C" void calc_EFG_LES_forcing_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                                       double* Jacobian, double* QJ, double* E, double* F, double* G, \
                                       double* fx, double* fy, double* fz);

extern "C" void calc_R_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                         double* E, double* F, double* G, double* R);
extern "C" void calc_R_forcing_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                                 double* E, double* F, double* G, \
                                 double* fx, double* fy, double* fz, double* R);
extern "C" void calc_step1_c(int nx, int ny, int nz, double coef, double* dx, double* dy, double* dz, \
                             double* E, double* F, double* G, double* Q, double* Q2);
extern "C" void calc_step1_forcing_c(int nx, int ny, int nz, double coef, double* dx, double* dy, double* dz, \
                                     double* E, double* F, double* G, \
                                     double* fx, double* fy, double* fz, double* Q, double* Q2);
extern "C" void calc_step_c(int nx, int ny, int nz, double coef1, double coef2, double* dx, double* dy, double* dz, \
                            double* E, double* F, double* G, double* Q, double* Q2, double* Rs);
extern "C" void calc_step_forcing_c(int nx, int ny, int nz, double coef1, double coef2, double* dx, double* dy, double* dz, \
                                    double* E, double* F, double* G, \
                                    double* fx, double* fy, double* fz, double* Q, double* Q2, double* Rs);
extern "C" void calc_step2_3_c(int nx, int ny, int nz, double coef1, double coef2, double coef3, double coef4, \
                               double* dx, double* dy, double* dz, double* E, double* F, double* G, \
                               double* Qin, double* Qinout);
extern "C" void calc_step2_3_forcing_c(int nx, int ny, int nz, double coef1, double coef2, double coef3, double coef4, \
                                       double* dx, double* dy, double* dz, double* E, double* F, double* G, \
                                       double* fx, double* fy, double* fz, double* Qin, double* Qinout);
extern "C" void calc_step4_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                             double* E, double* F, double* G, double* Rs, double* Q);
extern "C" void calc_step4_forcing_c(int nx, int ny, int nz, double* dx, double* dy, double* dz, \
                                     double* E, double* F, double* G, \
                                     double* fx, double* fy, double* fz, double* Rs, double* Q);
extern "C" void calc_error_c(int nx, int ny, int nz, double* R1_ptr, double* R2, double* R1_new, double* R2_new);
extern "C" void calc_Gauss_step_c(int nx, int ny, int nz, double a1, double a2, \
                                  double* dx, double* dy, double* dz, double* R1, double* R2, double* Q, double* Q2);
extern "C" void calc_Gauss_step_Q_c(int nx, int ny, int nz, double a1, double a2, \
                                    double* dx, double* dy, double* dz, double* R1, double* R2, double* Q);

