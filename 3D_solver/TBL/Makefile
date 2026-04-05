FC = mpif90

subroutines = name:calc_tau_straight,calc_tau_cross,flux4

OUTPUT_OPTIN = -cuda -acc -fast -Mchkptr -Minline=reshape,$(subroutines) -Minfo=accel,inline -Mchkptr -gpu=ptxinfo,rdc,lto -lnvhpcwrapnvtx 

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = shock_tube
vpath %f90 ../src:../../src

OBJ =	mod_globals.o \
      mod_constant.o \
      cpu_gpu_mpi.o \
      set_coordinate.o \
      calc_physical_quantities.o \
      calc_visc2.o \
      calc_visc4.o \
      calc_les.o \
      calc_hybrid.o \
      calc_muscl.o \
      calc_keep_kernel.o \
      calc_keep_kernel_internal.o \
      calc_slau_kernel.o \
      calc_slau_kernel_internal.o \
      calc_roe_kernel.o \
      calc_hybrid_kernel.o \
      calc_flux_base.o \
      calc_steps.o \
      calc_rescale.o \
      calc_para.o \
      set_bc_common.o \
      set_compressible_bl.o \
      set_init_common.o \
      set.o \
      print.o \
      preprocess.o \
      calc_time_dev.o \
      main.o

a.out: $(OBJ)
		$(FC) -cuda -acc -fast -Minfo -gpu=rdc,lto -lnvhpcwrapnvtx -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

mod_constant.o: mod_globals.mod
set_bc_common.o: mod_globals.mod
set_init_common.o: mod_globals.mod mod_constant.mod set_compressible_bl.mod
set_compressible_bl.o: mod_globals.mod mod_constant.mod
set.o: mod_globals.mod mod_constant.mod calc_para.mod set_bc_common.mod set_init_common.mod set_coordinate.mod
calc_physical_quantities.o: mod_globals.mod mod_constant.mod
calc_visc2.o: mod_globals.mod mod_constant.mod
calc_visc4.o: mod_globals.mod mod_constant.mod
calc_muscl.o: mod_globals.mod
calc_keep_kernel.o: mod_globals.mod calc_keep_3d.f90
calc_keep_kernel_internal.o: mod_globals.mod calc_keep_3d.f90
calc_slau_kernel.o: mod_globals.mod mod_constant.mod calc_muscl.mod calc_hybrid.mod calc_slau_3d.f90
calc_slau_kernel_internal.o: mod_globals.mod mod_constant.mod calc_muscl.mod calc_hybrid.mod calc_slau_3d.f90
calc_roe_kernel.o: mod_globals.mod mod_constant.mod calc_muscl.mod calc_hybrid.mod calc_roe_3d.f90
calc_hybrid_kernel.o: mod_globals.mod mod_constant.mod calc_keep_3d.f90 calc_muscl.mod calc_hybrid.mod calc_slau_3d.f90
calc_flux_base.o: mod_globals.mod calc_physical_quantities.mod calc_hybrid.mod calc_visc2.mod calc_visc4.mod calc_les.mod set.mod calc_keep_kernel.mod calc_slau_kernel.mod calc_roe_kernel.mod calc_hybrid.mod
calc_hybrid.o: mod_globals.mod
calc_steps.o: mod_globals.mod
calc_rescale.o: mod_globals.mod mod_constant.mod cpu_gpu_mpi.mod
preprocess.o: mod_globals.mod print.mod
calc_time_dev.o: mod_globals.mod preprocess.mod set.mod calc_flux_base.mod calc_steps.mod calc_rescale.mod calc_para.mod print.mod
print.o: mod_globals.mod
main.o: main.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

clean:
		rm -f *.mod *.o *.out

