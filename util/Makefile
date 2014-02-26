.SUFFIXES :
.SUFFIXES : .o .c

CC=gcc

.c.o:
	$(CC) -c $(CFLAGS) $<

objects = vector.o wanndata.o poscar.o mapping.o


all: extend_ham preprocess reorder dope_ham

extend_ham: extend_ham.c $(objects)
	$(CC) $(CFLAGS) -o extend_ham.x extend_ham.c $(objects)
preprocess: preprocess.c $(objects)
	$(CC) $(CFLAGS) -o preprocess.x preprocess.c $(objects)
reorder: reorder.c $(objects)
	$(CC) $(CFLAGS) -o reorder.x reorder.c $(objects)
dope_ham: dope_ham.c $(objects)
	$(CC) $(CFLAGS) -o dope_ham.x dope_ham.c $(objects)

clean:
	rm -rf *.o *.mod *.x
