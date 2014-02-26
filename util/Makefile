include ../make.sys

objects = vector.o wanndata.o poscar.o mapping.o


all: extend_ham preprocess reorder dope_ham

extend_ham: extend_ham.c $(objects)
	$(CC) $(CFLAGS) -o extend_ham.x extend_ham.c $(objects) -lm
preprocess: preprocess.c $(objects)
	$(CC) $(CFLAGS) -o preprocess.x preprocess.c $(objects) -lm
reorder: reorder.c $(objects)
	$(CC) $(CFLAGS) -o reorder.x reorder.c $(objects) -lm
dope_ham: dope_ham.c $(objects)
	$(CC) $(CFLAGS) -o dope_ham.x dope_ham.c $(objects) -lm

clean:
	rm -rf *.o *.mod *.x
