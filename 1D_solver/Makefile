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
	calc_Sutherland.o \
	calc_visc.o \
	calc_flux.o \
	calc_steps.o \
	print.o \
	calc_time_dev.o \
	main.o

a.out: $(OBJ)
		$(FC) -fast -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set_init.o: mod_globals.mod
calc_flux.o: calc_term.mod
calc_visc.o: calc_Sutherland.mod
calc_time_dev.o: mod_globals.mod calc_physical_quantities.mod set_bc.mod calc_flux.mod calc_visc.mod calc_steps.mod
main.o: main.f90 mod_globals.mod print.mod set_bc.mod set_init.mod \
		calc_flux.mod calc_steps.mod calc_time_dev.mod
