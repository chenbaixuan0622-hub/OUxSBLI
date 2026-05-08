FC = mpif90
OUTPUT_OPTIN = -cuda -acc -fast -Mpreprocess -Mchkptr -Minfo=accel,inline -gpu=ptxinfo,rdc,lto -lnvhpcwrapnvtx

.SUFFIXES : .f90
%.o: %.f90
	$(FC) -c $(OUTPUT_OPTIN) $<

TARGET = a.out
vpath %f90 ../src:../../src:../../3D_solver/src

OBJ = mod_globals.o \
      mod_constant.o \
      cpu_gpu_mpi.o \
      set_coordinate.o \
      calc_physical_quantities.o \
      calc_muscl.o \
      calc_hybrid.o \
      calc_hybrid_curv.o \
      calc_keep_kernel_curv.o \
      calc_slau_kernel_curv.o \
      calc_hybrid_kernel_curv.o \
      calc_steps_curv.o \
      load_smem_visc2_curv.o \
      calc_visc2_curv.o \
      set_bc_common.o \
      set.o \
      calc_flux_base_curv.o \
      print.o \
      print_curv.o \
      preprocess_curv.o \
      calc_time_dev_curv.o \
      main_curv.o

$(TARGET): $(OBJ)
	$(FC) -cuda -acc -fast -Minfo -gpu=rdc,lto -lnvhpcwrapnvtx -I./ $^

mod_constant.o: mod_globals.mod
set_coordinate.o: mod_globals.mod
set_bc_common.o: mod_globals.mod
set.o: mod_globals.mod mod_constant.mod set_coordinate.mod set_bc_common.mod
calc_physical_quantities.o: mod_globals.mod mod_constant.mod
calc_keep_kernel_curv.o: mod_globals.mod calc_keep_3d.f90
calc_slau_kernel_curv.o: mod_globals.mod mod_constant.mod calc_hybrid.mod calc_slau_3d.f90
calc_hybrid_kernel_curv.o: mod_globals.mod calc_keep_3d.f90 calc_hybrid_curv.mod calc_slau_3d.f90 calc_muscl.mod calc_hybrid.mod
calc_steps_curv.o: mod_globals.mod mod_constant.mod
calc_flux_base_curv.o: mod_globals.mod calc_physical_quantities.mod calc_hybrid_curv.mod calc_keep_kernel_curv.mod calc_slau_kernel_curv.mod calc_hybrid_kernel_curv.mod
preprocess_curv.o: mod_globals.mod print.mod print_curv.mod calc_flux_base_curv.mod
calc_time_dev_curv.o: mod_globals.mod preprocess_curv.mod set.mod calc_flux_base_curv.mod calc_steps_curv.mod print.mod print_curv.mod
set.o: mod_globals.mod mod_constant.mod set_bc_common.mod set_coordinate.mod
print.o: mod_globals.mod mod_constant.mod
print_curv.o: mod_globals.mod print.mod
main_curv.o: mod_globals.mod set_coordinate.mod set.mod calc_time_dev_curv.mod

clean:
	rm -f *.mod *.o *.out a.out
