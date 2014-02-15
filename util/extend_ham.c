#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#define MAXLEN 512

typedef struct __vector {
  double x[3];
} vector;

typedef struct __wanndata {
  int norb;
  int nrpt;
  double complex * ham;
  vector * rvec;
  int * weight;
} wanndata;

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

int locate_rpt(wanndata * wann, int nr[3]) {
  int ii;
  for(ii=0; ii<wann->nrpt; ii++) {
    if( (wann->rvec+ii)->x[0]==nr[0] &&
        (wann->rvec+ii)->x[1]==nr[1] &&
        (wann->rvec+ii)->x[2]==nr[2] )
      return ii;
  }
  return -1;
}

void extend_wann(wanndata * sc, wanndata * uc, int nx, int ny, int nz) {
  int irpt, iorb, jorb;
  int iirpt, iiorb, jjorb;
  int ix[3], jx[3], nr[3];
  int nnr[3];
  int ii, jj;

  sc->norb=uc->norb*nx*ny*nz;            /*  Set norb */

  for(ii=0; ii<3; ii++) nnr[ii]=0;        /*  set nrpt */
  for(irpt=0; irpt<uc->nrpt; irpt++) {
    for(ii=0; ii<3; ii++)
      if ((uc->rvec+irpt)->x[ii]>nnr[ii]) nnr[ii]=(uc->rvec+irpt)->x[ii];
  }

  nr[0]=nnr[0]/nx;
  nr[1]=nnr[1]/ny;
  nr[2]=nnr[2]/nz;

  sc->nrpt=0;
  for(irpt=0; irpt<uc->nrpt; irpt++) {
    if ((fabs((uc->rvec+irpt)->x[0])<=nr[0]) &&
        (fabs((uc->rvec+irpt)->x[1])<=nr[1]) &&
        (fabs((uc->rvec+irpt)->x[2])<=nr[2]))
      sc->nrpt++;
  }

  init_wanndata(sc);

  printf(" ##: Supercell hamiltonian dimensions determined.\n");

  ii=0;                                  /*  Set weight & rvec */
  for(irpt=0; irpt<uc->nrpt; irpt++) {
    if ((fabs((uc->rvec+irpt)->x[0])<=nr[0]) &&
        (fabs((uc->rvec+irpt)->x[1])<=nr[1]) &&
        (fabs((uc->rvec+irpt)->x[2])<=nr[2])) {
      (sc->rvec+ii)->x[0]=(uc->rvec+irpt)->x[0];
      (sc->rvec+ii)->x[1]=(uc->rvec+irpt)->x[1];
      (sc->rvec+ii)->x[2]=(uc->rvec+irpt)->x[2];
      sc->weight[ii]=uc->weight[irpt];
      ii++;
    }
  }

  printf(" ##: weight & rvecs determined.\n");
  printf(" ##:  sc->nrpt: %5d\n", sc->nrpt);

  for(irpt=0; irpt<sc->nrpt; irpt++) {
    for(iorb=0; iorb<sc->norb; iorb++) {
      iiorb=iorb%uc->norb;
      ix[2]=(iorb/uc->norb)%nz;
      ix[1]=(iorb/(uc->norb*nz))%ny;
      ix[0]=(iorb/(uc->norb*nz*ny))%nx;
      for(jorb=0; jorb<sc->norb; jorb++) {
        jjorb=jorb%uc->norb;
        jx[2]=(jorb/uc->norb)%nz;
        jx[1]=(jorb/(uc->norb*nz))%ny;
        jx[0]=(jorb/(uc->norb*nz*ny))%nx;


        nr[0]=((sc->rvec+irpt)->x[0])*nx+(ix[0]-jx[0]);
        nr[1]=((sc->rvec+irpt)->x[1])*ny+(ix[1]-jx[1]);
        nr[2]=((sc->rvec+irpt)->x[2])*nz+(ix[2]-jx[2]);

        if ((fabs(nr[0])>nnr[0]) ||
            (fabs(nr[1])>nnr[1]) ||
            (fabs(nr[2])>nnr[2]))
          sc->ham[irpt*sc->norb*sc->norb+iorb*sc->norb+jorb]=0.0;
        else {
          iirpt=locate_rpt(uc, nr);
          if(iirpt==-1) {
            printf("!!!ERROR: Cannot locate rpt for:\n");
            printf(" iorb: %5d, jorb: %5d, vec:(%5d,%5d,%5d) ==> (%5d,%5d,%5d)\n", iorb, jorb, (int)(sc->rvec+irpt)->x[0], (int)(sc->rvec+irpt)->x[1], (int)(sc->rvec+irpt)->x[2], nr[0], nr[1], nr[2]);
            exit(0);
          }
          else
            sc->ham[irpt*sc->norb*sc->norb+iorb*sc->norb+jorb]=uc->ham[iirpt*uc->norb*uc->norb+iiorb*uc->norb+jjorb];
        }
      }
    }
  }
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

int main(int argc, char ** argv) {
  wanndata wann_uc, wann_sc;
  int nkx_uc, nky_uc, nkz_uc;
  int nkx_sc, nky_sc, nkz_sc;
  int nx, ny, nz;

  if(argc<5) {
    printf(" Usage:%s <SEED> nx ny nz\n", argv[0]);
    exit(0);
  }

  sscanf(argv[2], " %d", &nx);
  sscanf(argv[3], " %d", &ny);
  sscanf(argv[4], " %d", &nz);

  read_ham(&wann_uc, argv[1]);
  printf(" %5d x %5d Hamiltonian read for %6d rpts.\n", wann_uc.norb, wann_uc.norb, wann_uc.nrpt);
  printf("   Extending it to %5d x %5d x %5d supercell.\n", nx, ny, nz);
  extend_wann(&wann_sc, &wann_uc, nx, ny, nz);
  write_ham(&wann_sc);
}
