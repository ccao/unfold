#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "vector.h"
#include "wanndata.h"
#include "mapping.h"

void read_mapping(mapping * map, FILE * fin) {
  char line[MAXLEN];
  int ii, norb;

  for(ii=0; ii<3; ii++) 
    fgets(line, MAXLEN, fin);

  fgets(line, MAXLEN, fin);
  sscanf(line, " %*d %d", &norb);

  for(ii=0; ii<norb; ii++) {
    fgets(line, MAXLEN, fin);
    sscanf(line, " %d %lf %lf %lf", &((map+ii)->nat), ((map+ii)->rvec).x, ((map+ii)->rvec).x+1, ((map+ii)->rvec).x+2);
    (map+ii)->nat--;
  }
}


void extend_wann(wanndata * sc, wanndata * uc, mapping * map, int n[3]) {
  int irpt, iorb, jorb;
  int iirpt, iiorb, jjorb;
  int nnr[3];
  int ii, jj;
  vector vr, nr;

  sc->norb=uc->norb*n[0]*n[1]*n[2];       /*  Set norb */

  for(ii=0; ii<3; ii++) nnr[ii]=0;        /*  set nrpt */
  for(irpt=0; irpt<uc->nrpt; irpt++) {
    for(ii=0; ii<3; ii++)
      if ((uc->rvec+irpt)->x[ii]>nnr[ii]) nnr[ii]=(uc->rvec+irpt)->x[ii];
  }

  for(ii=0; ii<3; ii++)
    nr.x[ii]=(int)(nnr[ii]/n[ii]);

  sc->nrpt=0;
  for(irpt=0; irpt<uc->nrpt; irpt++) {
    if (isenclosed(uc->rvec[irpt], nr))
      sc->nrpt++;
  }

  init_wanndata(sc);

  ii=0;                                  /*  Set weight & rvec */
  for(irpt=0; irpt<uc->nrpt; irpt++) {
    if (isenclosed(uc->rvec[irpt],nr)) {
      (sc->rvec+ii)->x[0]=(uc->rvec+irpt)->x[0];
      (sc->rvec+ii)->x[1]=(uc->rvec+irpt)->x[1];
      (sc->rvec+ii)->x[2]=(uc->rvec+irpt)->x[2];
      sc->weight[ii]=uc->weight[irpt];
      ii++;
    }
  }

  for(irpt=0; irpt<sc->nrpt; irpt++) {
    vector_multiply(&nr, sc->rvec[irpt], n);

    for(iorb=0; iorb<sc->norb; iorb++) {
      iiorb=(map+iorb)->nat;

      for(jorb=0; jorb<sc->norb; jorb++) {
        jjorb=(map+jorb)->nat;

        vector_add(&vr, nr, (map+iorb)->rvec);
        vector_sub(&vr, vr, (map+jorb)->rvec);

        if ((fabs(vr.x[0])>nnr[0]) ||
            (fabs(vr.x[1])>nnr[1]) ||
            (fabs(vr.x[2])>nnr[2]))
          sc->ham[irpt*sc->norb*sc->norb+iorb*sc->norb+jorb]=0.0;
        else {
          iirpt=locate_rpt(uc, vr);
          if(iirpt==-1) {
            printf("!!!ERROR: Cannot locate rpt for:\n");
            printf(" iorb: %5d, jorb: %5d, vec:(%5d,%5d,%5d) ==> (%5d,%5d,%5d)\n", iorb, jorb, 
              (int)(sc->rvec+irpt)->x[0], (int)(sc->rvec+irpt)->x[1], (int)(sc->rvec+irpt)->x[2], 
              (int)vr.x[0], (int)vr.x[1], (int)vr.x[2]);
            exit(0);
          }
          else
            sc->ham[irpt*sc->norb*sc->norb+iorb*sc->norb+jorb]=uc->ham[iirpt*uc->norb*uc->norb+iiorb*uc->norb+jjorb];
        }
      }
    }
  }
}

int main(int argc, char ** argv) {
  wanndata wann_uc, wann_sc;
  mapping * map;
  FILE * fin;
  int n[3];

  if(argc<5) {
    printf(" Usage:%s <SEED> nx ny nz\n", argv[0]);
    exit(0);
  }

  sscanf(argv[2], " %d", n);
  sscanf(argv[3], " %d", n+1);
  sscanf(argv[4], " %d", n+2);

  read_ham(&wann_uc, argv[1]);
  map=(mapping *) malloc(sizeof(mapping)*wann_uc.norb*n[0]*n[1]*n[2]);

  fin=fopen("unfold.map", "r");
  read_mapping(map, fin);
  fclose(fin);

  extend_wann(&wann_sc, &wann_uc, map, n);
/*
  write_ham(&wann_sc);
*/
  write_reduced_ham(&wann_sc);

  finalize_wanndata(wann_uc);
  finalize_wanndata(wann_sc);
}
