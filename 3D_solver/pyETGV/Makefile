# Compilers
FC             = nvfortran
CC             = nvcc
# CUDA
CUDA_PATH      = /work/opt/local/aarch64/cores/nvidia/24.9/Linux_aarch64/24.9
CUDA_INC       = -I$(CUDA_PATH)/include -I$(CUDA_PATH)/compilers/include
CUDA_LIB       = -L$(CUDA_PATH)/lib64 -L$(CUDA_PATH)/compilers/lib64 -lcudart
# Python and pybind11
PYTHON_VERSION = 3.10
PYTHON_INC     = $(shell python$(PYTHON_VERSION)-config --includes)
PYBIND11_INC   = $(shell python3 -m pybind11 --includes)
# Flags
FFLAGS         = -fPIE -cuda -acc -fast
CXXFLAGS       = -Xcompiler -fPIC -O3 -std=c++17 $(PYTHON_INC) $(PYBIND11_INC) $(CUDA_INC)
LDFLAGS        = -Mnomain -cuda -acc -shared -fPIC $(CUDA_LIB) -lstdc++

vpath %f90 ../src:../../src
vpath %cu  ../src
vpath %h 	 ../src

# Source files
FSRC =mod_globals.o \
      mod_constant.o \
      set_coordinate.o \
      calc_common_dim.o \
      calc_mat.o \
      calc_physical_quantities.o \
      calc_sutherland.o \
      calc_visc.o \
      calc_les.o \
      calc_common.o \
      calc_hybrid.o \
      calc_slau.o \
      calc_term.o \
      calc_keep.o \
			calc_roe.o \
      calc_muscl.o \
			set_bc_common.o \
			set.o \
      calc_flux.o \
      calc_flux_base.o \
      calc_flux_base_c.o \
      calc_steps.o \
      calc_steps_c.o
CSRC 	= wrapper.cu

FOBJ  = $(FSRC:.f90=.o)
COBJ  = $(CSRC:.cu=.o)
OBJ   = $(FOBJ) $(COBJ)

TARGET = cufd.so
all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(LDFLAGS) -o $@ $(OBJ) -L$(CUDA_PATH)/lib64 -lcudart

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.cu
	$(CC) $(CXXFLAGS) -c $< -o $@


wrapper.o: $(FSRC:.f90=.o)
cufd.so: wrapper.o

mod_constant.o: mod_globals.mod
set_bc_common.o: mod_globals.mod mod_constant.mod
set.o: mod_globals.mod mod_constant.mod set_bc_common.mod set_coordinate.mod
calc_physical_quantities.o: mod_globals.mod mod_constant.mod
calc_common.o: mod_globals.mod
calc_term.o: calc_common_dim.mod
calc_keep.o: mod_globals.mod mod_constant.mod calc_common_dim.mod calc_term.mod calc_mat.mod
calc_sutherland.o: mod_globals.mod mod_constant.mod
calc_visc.o: mod_globals.mod mod_constant.mod calc_sutherland.mod
calc_mat.o: mod_globals.mod calc_common_dim.mod
calc_slau.o: mod_globals.mod mod_constant.mod calc_common.mod calc_hybrid.mod
calc_roe.o: mod_globals.mod calc_common.mod calc_common_dim.mod calc_physical_quantities.mod calc_mat.mod
calc_flux.o: mod_globals.mod calc_keep.mod calc_slau.mod calc_roe.mod calc_hybrid.mod calc_muscl.mod
calc_flux_base.o: mod_globals.mod calc_physical_quantities.mod calc_hybrid.mod calc_flux.mod calc_visc.mod calc_les.mod set.mod
calc_hybrid.o: mod_globals.mod
calc_steps.o: mod_globals.mod

clean:
	rm -f *.o *.so

