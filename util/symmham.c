#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "constants.h"
#include "wanndata.h"
#include "vector.h"

int main(int argc, char ** argv) {

  wanndata in;
  wanndata * out;
  wanndata final;
  vector * symm;
  vector * site;
  int * info;
  int nwann, nsymm;
  int mode;

  FILE * fin;
  char line[MAXLEN];
  char seed[MAXLEN];
  int ii, jj;

  if (argc<3) {
    printf("Usage: %s <INPUTFILE> MODE\n", argv[0]);
    printf("   MODE: 0: output Hamiltonian after symmetry operation.\n");
    printf("         1: symmetrize over the operation.\n");
    exit(0);
  }

  fin=fopen(argv[1], "r");
  fgets(line, MAXLEN, fin);
  sscanf(line, " %s %d", seed, &nwann);

  site=(vector *) malloc(sizeof(vector)*nwann);
  info=(int *) malloc(sizeof(int)*nwann);

  for (ii=0; ii<nwann; ii++) {
    fgets(line, MAXLEN, fin);
    sscanf(line, " %lf %lf %lf %d", &(site[ii].x, site[ii].x+1, site[ii].x+2, info+ii);
  }

  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &nsymm);

  symm=(vector *)malloc(sizeof(vector)*4*nsymm);
  for (ii=0; ii<nsymm; ii++) {
    for (jj=0; jj<4; jj++) {
      fgets(line, MAXLEN, fin);
      sscanf(line, " %lf %lf %lf ", (symm+ii*4+jj)->x, (symm+ii*4+jj)->x+1, (symm+ii*4+jj)->x+2);
    }
    fgets(line, MAXLEN, fin);
  }

  fclose(fin);

  sscanf(argv[2], " %d", &mode);
  if (mode==0) {
    printf(" Will generate %5d hamiltonians.\n", nsymm);
  }
  else {
    printf(" Will symmetrize the Hailtonian.\n");
  }

  read_ham(&in, seed);
  if (nwann!=in.norb) {
    printf(" !!! ERROR: Wrong number of wannier sites\n");
    exit(0);
  }

  out=(wanndata *)malloc(sizeof(wanndata)*nsymm);

  for (ii=0; ii<nsymm; ii++) {
    out[ii].nrpt=in.nrpt;
    out[ii].norb=in.norb;
    init_wanndata(out+ii);
    rotate_ham(out+ii, &in, symm+ii*4, site, info);
  }

  if (mode==0) {
    for (ii=0; ii<nsymm; ii++) {
      sprintf(seed, "output%d", ii+1);
      write_ham(out+ii, seed);
    }
  }
  else {
    final.nrpt=in.nrpt;
    final.norb=in.norb;
    init_wanndata(&final);
    memcpy(final.rvec, in.rvec, sizeof(vector)*in.nrpt);
    memcpy(final.weight, in.weight, sizeof(int)*in.nrpt);
    for (ii=0; ii<in.nrpt*in.norb*in.norb; ii++) {
      final.ham[ii]=0.0;
      for (jj=0; jj<nsymm; jj++) {
        final.ham[ii]+=out[jj].ham[ii];
      }
      final.ham[ii]/=nsymm;
    }
    write_ham(&final, "output");
    finalize_wanndata(final);
  }

  for (ii=0; ii<nsymm; ii++) finalize_wanndata(out[ii]);
  finalize_wanndata(in);
  free(symm);
  free(site);

  return 0;
}
