#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "vector.h"
#include "wanndata.h"

void extend_wann(wanndata * sc, wanndata * uc, int nx, int ny, int nz) {
  int irpt, iorb, jorb;
  int iirpt, iiorb, jjorb;
  int ix[3], jx[3], nr[3];
  int nnr[3];
  int ii, jj;
  vector vr;

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
          vr.x[0]=nr[0];
          vr.x[1]=nr[1];
          vr.x[2]=nr[2];
          iirpt=locate_rpt(uc, vr);
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
  extend_wann(&wann_sc, &wann_uc, nx, ny, nz);
  write_ham(&wann_sc);

  finalize_wanndata(wann_uc);
  finalize_wanndata(wann_sc);
}
