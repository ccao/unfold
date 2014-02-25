include ../make.sys

objects = input.o calc_spectrum.o output_spectrum.o

LIBS = ../modules/mod.a

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

unfold: unfold.f90 $(objects) $(LIBS)
	$(F90) $(F90FLAGS) -cpp -o unfold.x unfold.f90 $(objects) $(LIBS) $(INCLUDES) $(LAPACK95LIBS) $(LAPACKLIBS)

clean:
	rm -rf *.o *.mod *.x
