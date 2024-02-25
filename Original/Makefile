FC = gfortran

OUTPUT_OPTIN = -O2

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = FORTRAN_MPI
OBJ =	mod_globals.o \
		set_bc.o \
		set_init.o \
		calc_term.o \
		calc_flux.o \
		calc_steps.o \
		calc_time_dev.o \
		print.o \
		main.o

a.out: $(OBJ)
		$(FC) -O2 -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set_init.o: mod_globals.mod set_bc.mod
calc_flux.o: calc_term.mod
calc_time_dev.o: mod_globals.mod set_bc.mod calc_flux.mod calc_steps.mod
main.o: print.mod set_bc.mod set_init.mod \
		calc_flux.mod calc_steps.mod calc_time_dev.mod