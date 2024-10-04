FC = gfortran

OUTPUT_OPTIN = -Ofast

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = TGV
vpath %f90 ../src

OBJ =	mod_globals.o \
      set_coordinate.o \
      set.o \
      calc_physical_quantities.o \
      calc_term.o \
      calc_keep.o \
      calc_flux.o \
      calc_steps.o \
      print.o \
      calc_time_dev.o \
      main.o

a.out: $(OBJ)
		$(FC) -Ofast -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set.o: mod_globals.mod
calc_physical_quantities.o: mod_globals.mod
calc_keep.o: mod_globals.mod calc_term.mod
calc_flux.o: mod_globals.mod calc_keep.mod
calc_steps.o: mod_globals.mod
calc_time_dev.o: mod_globals.mod calc_physical_quantities.mod set.mod calc_flux.mod calc_steps.mod print.mod
print.o: mod_globals.mod
main.o: main.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

