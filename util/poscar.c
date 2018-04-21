#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "vector.h"
#include "poscar.h"

void init_poscar(poscar * psc) {
  psc->nsp=0;
  psc->nat_per_sp=NULL;
  psc->tau=NULL;
}

void read_poscar_header(poscar * psc, FILE * fin) {
  char line[MAXLEN];
  char * p;
  double lat;
  int ii;

  fgets(line, MAXLEN, fin);
  fgets(line, MAXLEN, fin);
  sscanf(line, " %lf", &lat);

  for(ii=0; ii<3; ii++) {
    fgets(line, MAXLEN, fin);
    sscanf(line, " %lf %lf %lf", (psc->cell+ii)->x, (psc->cell+ii)->x+1, (psc->cell+ii)->x+2);
    (psc->cell[ii]).x[0]*=lat;
    (psc->cell[ii]).x[1]*=lat;
    (psc->cell[ii]).x[2]*=lat;
  }

  fgets(line, MAXLEN, fin);
  ii=0;
  p=strtok(line, " ");
  while(p!=NULL) {
    ii++;
    p=strtok(NULL, " ");
  }
  psc->nsp=ii;

  psc->nat_per_sp=(int *) malloc(sizeof(int)*psc->nsp);

  fgets(line, MAXLEN, fin);

  psc->nat=0;
  p=strtok(line, " ");
  for(ii=0; ii<psc->nsp; ii++) {
    sscanf(p, " %d", psc->nat_per_sp+ii);
    psc->nat+=psc->nat_per_sp[ii];
    p=strtok(NULL, " ");
  }

  psc->tau=NULL;

}

void read_poscar(poscar * psc, char * fn) {
  FILE * fin;
  char line[MAXLEN];
  char * p;
  int ii;

  fin=fopen(fn, "r");

  read_poscar_header(psc, fin);
  printf("Done after read header\n");
  printf("Num of species: %d\n", psc->nsp);
  printf("Num of atoms: %d\n", psc->nat);

  psc->tau=(vector *)malloc(sizeof(vector)*psc->nat);
  fgets(line, MAXLEN, fin);
  for(ii=0; ii<psc->nat; ii++) {
    fgets(line, MAXLEN, fin);
    sscanf(line, " %lf %lf %lf ", (psc->tau+ii)->x, (psc->tau+ii)->x+1, (psc->tau+ii)->x+2);
  }
  fclose(fin);
}

void finalize_poscar(poscar psc) {
  if(psc.nat_per_sp)
    free(psc.nat_per_sp);
  if(psc.tau)
    free(psc.tau);
}
