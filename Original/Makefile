FC = nvfortran

OUTPUT_OPTIN = -fast

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = FORTRAN_MPI
OBJ =	mod_globals.o \
		set_bc.o \
		set_init.o \
		calc_physical_quantities.o \
		calc_term.o \
		calc_KEEP.o \
		calc_SLAU.o \
		calc_visc.o \
		calc_steps.o \
		calc_time_dev.o \
		print.o \
		main.o

a.out: $(OBJ)
		$(FC) -fast -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set_init.o: mod_globals.mod set_bc.mod
calc_flux.o: calc_term.mod
calc_time_dev.o: mod_globals.mod set_bc.mod calc_physical_quantities.mod calc_KEEP.mod calc_SLAU.mod calc_visc.mod calc_steps.mod
main.o: print.mod set_bc.mod set_init.mod calc_time_dev.mod
