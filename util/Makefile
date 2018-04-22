include ../make.sys

objects = vector.o wanndata.o poscar.o mapping.o read_orbdef.o


all: extend_ham preprocess unfold2plot dope_ham mixham

extend_ham: extend_ham.c $(objects)
	$(CC) $(CFLAGS) -o extend_ham.x extend_ham.c $(objects) -lm
preprocess: preprocess.c $(objects)
	$(CC) $(CFLAGS) -o preprocess.x preprocess.c $(objects) -lm
unfold2plot: unfold2plot.c $(objects)
	$(CC) $(CFLAGS) -o unfold2plot.x unfold2plot.c $(objects) -lm
dope_ham: dope_ham.c $(objects)
	$(CC) $(CFLAGS) -o dope_ham.x dope_ham.c $(objects) -lm

mixham: mixham.c $(objects)
	$(CC) $(CFLAGS) -o mixham.x mixham.c $(objects) -lm

clean:
	rm -rf *.o *.mod *.x
