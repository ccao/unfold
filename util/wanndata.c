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

void write_ham(wanndata * wann, char * seed) {
  int irpt, iorb, jorb;
  FILE * fout;
  char fn[MAXLEN];

  if (fn==NULL) {
    fout=stdout;
  }
  else {
    sprintf(fn, "%s_hr.dat", seed);
    fout=fopen(fn, "w");
  }
  fprintf(fout, "# Wannier Hamiltonian extended to supercell\n");
  fprintf(fout, "%5d\n%5d", wann->norb, wann->nrpt);
  for(irpt=0; irpt<wann->nrpt; irpt++) {
    if(irpt%15==0) fprintf(fout, "\n");
    fprintf(fout, "%5d", wann->weight[irpt]);
  }
  fprintf(fout, "\n");

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    for(iorb=0; iorb<wann->norb; iorb++) {
      for(jorb=0; jorb<wann->norb; jorb++) {
        fprintf(fout, "%5d%5d%5d%5d%5d%22.16f%22.16f\n", 
          (int)(wann->rvec+irpt)->x[0], (int)(wann->rvec+irpt)->x[1], (int)(wann->rvec+irpt)->x[2], 
          jorb+1, iorb+1,
          creal(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]),
          cimag(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]));
      }
    }
  }
  if (fn)  fclose(fout);
}

void write_reduced_ham(wanndata * wann, char * seed) {
  int irpt, iorb, jorb;
  int iirpt;
  long long int nwann;
  vector rv;
  char fn[MAXLEN];

  FILE * fout;

  if (fn==NULL) {
    fout=stdout;
  }
  else {
    sprintf(fn, "%s_hr.dat", seed);
    fout=fopen(fn, "w");
  }

  fprintf(fout, "# Wannier Hamiltonian extended to supercell\n");
  fprintf(fout, "%5d\n%5d", wann->norb, wann->nrpt);
  for(irpt=0; irpt<wann->nrpt; irpt++) {
    if(irpt%15==0) fprintf(fout, "\n");
    fprintf(fout, "%5d", wann->weight[irpt]);
  }
  fprintf(fout, "\n");

  nwann=0;

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    fprintf(fout, "%5d%5d%5d\n", (int)(wann->rvec+irpt)->x[0], (int)(wann->rvec+irpt)->x[1], (int)(wann->rvec+irpt)->x[2]);
    for(iorb=0; iorb<wann->norb; iorb++)
      for(jorb=0; jorb<=iorb; jorb++)
        if(cabs(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])>eps8)
          nwann++;
  }

  fprintf(fout, "%ld\n", nwann);

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    rv.x[0]=-(wann->rvec+irpt)->x[0];
    rv.x[1]=-(wann->rvec+irpt)->x[1];
    rv.x[2]=-(wann->rvec+irpt)->x[2];
    iirpt=locate_rpt(wann, rv);
    for(iorb=0; iorb<wann->norb; iorb++) {
      for(jorb=0; jorb<=iorb; jorb++) {
        if(cabs(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])>eps8) {
          fprintf(fout, "%5d%5d%5d%5d%22.16f%22.16f\n",
            iorb+1, jorb+1, irpt+1, iirpt+1,
            creal(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]),
            cimag(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]));
          if(fabs(creal(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])-
                  creal(wann->ham[iirpt*wann->norb*wann->norb+jorb*wann->norb+iorb]))>eps8 ||
             fabs(cimag(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])+
                  cimag(wann->ham[iirpt*wann->norb*wann->norb+jorb*wann->norb+iorb]))>eps8) 
            fprintf(fout, "!!!!WARNING: Unhermitian Hamiltonian:%5d%5d%5d%5d%5d\n", iorb+1, jorb+1, (int)(wann->rvec+irpt)->x[0], (int)(wann->rvec+irpt)->x[1], (int)(wann->rvec+irpt)->x[2]);
        }
      }
    }
  }
  if (fn) fclose(fout);
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

