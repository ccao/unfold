#ifdef WANNDATA_H
#else
#define WANNDATA_H

#include <complex.h>

#include "vector.h"

typedef struct __wanndata {
  int norb;
  int nrpt;
  double complex * ham;
  vector * rvec;
  int * weight;
} wanndata;

void init_wanndata(wanndata * wann);
void read_ham(wanndata * wann, char * seed);
void finalize_wanndata(wanndata wann);
void write_ham(wanndata * wann);

#endif
