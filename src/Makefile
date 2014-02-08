include ../make.sys

objects = input.o calc_spectrum.o output_spectrum.o

LIBS = ../modules/mod.a

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

unfold: unfold.f90 $(objects) $(LIBS)
	$(F90) $(F90FLAGS) -fpp -o unfold.x unfold.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

clean:
	rm -rf *.o *.mod *.x
