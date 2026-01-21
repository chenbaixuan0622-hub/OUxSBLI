FC = mpif90
STUB_PATH    = /opt/nvidia/hpc_sdk/Linux_x86_64/24.7/cuda/12.5/lib64/stubs
OUTPUT_OPTIN = -cuda -fast -Minfo=accel -cudalib=cusparse -L$(STUB_PATH) -lcuda

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = ST
vpath %f90 ../src:../../src

OBJ =	mod_globals.o \
      mod_constant.o \
      calc_flux.o \
      calc_visc.o \
      calc_steps.o \
      calc_igr.o \
      set.o \
      print.o \
      calc_time_dev.o \
      main1d.o

a.out: $(OBJ)
		$(FC) -cuda -fast -Minfo -cudalib=cusparse -L$(STUB_PATH) -lcuda -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

mod_constant.o: mod_globals.mod
calc_flux.o: mod_constant.mod
calc_visc.o: mod_globals.mod mod_constant.mod
calc_steps.o: mod_constant.mod
calc_igr.o: mod_globals.mod
set.o: mod_globals.mod
calc_time_dev.o: mod_globals.mod mod_constant.mod calc_flux.mod calc_visc.mod calc_steps.mod calc_igr.mod set.mod print.mod
print.o: mod_globals.mod
main1d.o: main1d.f90 mod_globals.mod set.mod calc_time_dev.mod

clean:
		rm -f *.mod *.o *.out

