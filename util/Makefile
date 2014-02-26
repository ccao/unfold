include ../make.sys

objects = vector.o wanndata.o poscar.o


all: extend_ham preprocess reorder

extend_ham: extend_ham.c $(objects)
	$(CC) $(CFLAGS) -o extend_ham.x extend_ham.c $(objects) -lm
preprocess: preprocess.c $(objects)
	$(CC) $(CFLAGS) -o preprocess.x preprocess.c $(objects) -lm
reorder: reorder.c $(objects)
	$(CC) $(CFLAGS) -o reorder.x reorder.c $(objects) -lm

clean:
	rm -rf *.o *.mod *.x
