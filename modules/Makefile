include ../make.sys

objects = constants.o para.o wanndata.o mapdata.o specdata.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a
