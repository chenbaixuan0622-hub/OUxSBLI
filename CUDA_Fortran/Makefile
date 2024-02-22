FC = nvfortran

OUTPUT_OPTIN = -cuda -acc -fast -Minfo

.SUFFIXES : .cuf .f90

%.o: %.cuf
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.cuf %.f90 %.o
		@:

TARGET = CUDA_Fortran
OBJ =	set_bc.o \
		set_parallel.o \
		set_init.o \
		calc_term.o \
		calc_flux.o \
		calc_steps.o \
		set_bc_d.o \
		print.o \
		calc_time_dev.o \
		main.o

a.out: $(OBJ)
		$(FC) -cuda -acc -fast -Minfo -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set_init.o: set_parallel.mod set_bc.mod
print.o: set_parallel.mod
calc_flux.o: calc_term.mod
calc_steps.o: calc_flux.mod
calc_time_dev.o: set_bc_d.mod calc_steps.mod print.mod
main.o: set_parallel.mod set_bc.mod set_init.mod \
		set_bc_d.mod calc_term.mod calc_flux.mod calc_steps.mod calc_time_dev.mod