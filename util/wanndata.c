#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "constants.h"
#include "vector.h"
#include "wanndata.h"

void init_wanndata(wanndata * wann) {
  wann->ham=(double complex *) malloc(sizeof(double complex)*wann->norb*wann->norb*wann->nrpt);
  wann->rvec=(vector *) malloc(sizeof(vector)*wann->nrpt);
  wann->weight=(int *) malloc(sizeof(int)*wann->nrpt);
}

void read_ham(wanndata * wann, char * seed) {
  FILE * fin;
  int ii, iorb, jorb;
  char line[MAXLEN];
  char * pch;
  int t1, t2, t3;
  double a, b;

  strcpy(line, seed);
  strcat(line, "_hr.dat");

  fin=fopen(line, "r");

  fgets(line, MAXLEN, fin);
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &(wann->norb));
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &(wann->nrpt));

  init_wanndata(wann);

  for(ii=0; ii<wann->nrpt; ii++) {
    if(ii%15==0) {
      fgets(line, MAXLEN, fin);
      pch=strtok(line, " ");
    }
    else {
      pch=strtok(NULL, " ");
    }
    sscanf(pch, "%d", wann->weight+ii);
  }

  for(ii=0; ii<wann->nrpt; ii++) {
    for(iorb=0; iorb<wann->norb; iorb++) {
      for(jorb=0; jorb<wann->norb; jorb++) {
        fgets(line, MAXLEN, fin);
        sscanf(line, " %d %d %d %*d %*d %lf %lf ", &t1, &t2, &t3, &a, &b);
        if (iorb==0 && jorb==0) {
          (wann->rvec+ii)->x[0]=t1;
          (wann->rvec+ii)->x[1]=t2;
          (wann->rvec+ii)->x[2]=t3;
        }
        wann->ham[ii*wann->norb*wann->norb+iorb*wann->norb+jorb]=a+_Complex_I*b;
      }
    }
  }

  fclose(fin);
}

void finalize_wanndata(wanndata wann) {
  free(wann.rvec);
  free(wann.weight);
  free(wann.ham);
}

void write_ham(wanndata * wann) {
  int irpt, iorb, jorb;
  printf("# Wannier Hamiltonian extended to supercell\n");
  printf("%5d\n%5d", wann->norb, wann->nrpt);
  for(irpt=0; irpt<wann->nrpt; irpt++) {
    if(irpt%15==0) printf("\n");
    printf("%5d", wann->weight[irpt]);
  }
  printf("\n");

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    for(iorb=0; iorb<wann->norb; iorb++) {
      for(jorb=0; jorb<wann->norb; jorb++) {
        printf("%5d%5d%5d%5d%5d%12.6f%12.6f\n", 
          (int)(wann->rvec+irpt)->x[0], (int)(wann->rvec+irpt)->x[1], (int)(wann->rvec+irpt)->x[2], 
          iorb+1, jorb+1,
          creal(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]),
          cimag(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]));
      }
    }
  }
}

int locate_rpt(wanndata * wann, vector vr) {
  int ii;
  vector dx;
  
  for(ii=0; ii<wann->nrpt; ii++) {
    vector_sub(&dx, wann->rvec[ii], vr);
    if( fabs(dx.x[0])<eps &&
        fabs(dx.x[1])<eps &&
        fabs(dx.x[2])<eps )
      return ii;
  }
  return -1;
}

