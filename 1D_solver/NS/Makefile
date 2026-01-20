FC = mpif90

OUTPUT_OPTIN = -cuda -fast -Minfo=accel

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = ST
vpath %f90 ../src:../../src

OBJ =	mod_globals.o \
      mod_constant.o \
      set.o \
      print.o \
      calc_time_dev.o \
      main1d.o

a.out: $(OBJ)
		$(FC) -cuda -fast -Minfo -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

mod_constant.o: mod_globals.mod
set.o: mod_globals.mod
calc_time_dev.o: mod_globals.mod mod_constant.mod set.mod print.mod
print.o: mod_globals.mod
main1d.o: main1d.f90 mod_globals.mod set.mod calc_time_dev.mod

clean:
		rm -f *.mod *.o *.out

