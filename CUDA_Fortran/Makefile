FC = nvfortran

OUTPUT_OPTIN = -cuda -acc -gpu=cc75 -fast -Minfo

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = CUDA_Fortran
OBJ =	mod_variables.o \
		set_bc.o \
		set_parallel.o \
		set_init.o \
		calc_term.o \
		calc_flux.o \
		calc_steps.o \
		print.o \
		calc_time_dev.o \
		main.o

a.out: $(OBJ)
		$(FC) -cuda -acc -gpu=cc75 -fast -Minfo -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set_init.o: set_parallel.mod set_bc.mod
print.o: set_parallel.mod
calc_flux.o: mod_variables.mod calc_term.mod
calc_steps.o: calc_flux.mod
calc_time_dev.o: set_bc.mod calc_steps.mod print.mod
main.o: set_parallel.mod set_bc.mod set_init.mod \
		calc_term.mod calc_flux.mod calc_steps.mod calc_time_dev.mod