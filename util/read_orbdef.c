#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"

int read_orbdef(int ** norb_per_sp, char * fuc, char * fsc, FILE * fin) {
  int ii, nsp;
  char line[MAXLEN];
  char * p;

  fgets(line, MAXLEN, fin);
  p=line;
  while((*p)!='\r' && (*p)!='\0' && (*p)!='\n') {
    fuc[p-line]=(*p);
    p++;
  }
  fuc[p-line]=0;

  fgets(line, MAXLEN, fin);
  p=line;
  while((*p)!='\r' && (*p)!='\0' && (*p)!='\n') {
    fsc[p-line]=(*p);
    p++;
  }
  fsc[p-line]=0;

  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &nsp);
  (* norb_per_sp)=(int *) malloc(sizeof(int)*nsp);
  fgets(line, MAXLEN, fin);
  p=strtok(line, " ");
  for(ii=0; ii<nsp; ii++) {
    sscanf(p, "%d", (*norb_per_sp)+ii);
    p=strtok(NULL, " ");
  }
  return nsp;
}
