FC = mpif90

OUTPUT_OPTIN = -cuda -fast -Minfo=accel

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = shock_tube
vpath %f90 ../src:../../src

OBJ =	mod_globals.o \
      set.o \
      calc_common_dim.o \
      calc_sutherland.o \
      calc_visc.o \
      calc_common.o \
      calc_term.o \
      calc_kep1d.o \
      calc_slau1d.o \
      calc_flux.o \
      calc_steps.o \
      print1d.o \
      calc_time_dev.o \
      main1d.o

a.out: $(OBJ)
		$(FC) -cuda -fast -Minfo -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

set.o: mod_globals.mod
calc_common.o: mod_globals.mod
calc_term.o: calc_common_dim.mod
calc_kep1d.o: mod_globals.mod calc_common_dim.mod calc_term.mod
calc_slau1d.o: mod_globals.mod
calc_visc.o: mod_globals.mod calc_sutherland.mod
calc_flux.o: mod_globals.mod calc_kep1d.mod calc_slau1d.mod
calc_steps.o: mod_globals.mod
calc_time_dev.o: mod_globals.mod set.mod calc_flux.mod calc_visc.mod calc_steps.mod print1d.mod
print1d.o: mod_globals.mod
main1d.o: main1d.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

