#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "wanndata.h"

void shift_ham(wanndata * wann, double ef) {
  int ir, io;

  for(ir=0; ir<wann->nrpt; ir++)
    for(io=0; io<wann->norb; io++)
      wann->ham[ir*wann->norb*wann->norb+io*wann->norb+io]-=ef;

}

void mix_ham(wanndata * wann, wanndata * wann1, wanndata * wann2, double x) {
  int ir, io, jo;

  wann->nrpt=wann1->nrpt;
  wann->norb=wann1->norb;

  init_wanndata(wann);

  for(ir=0; ir<wann->nrpt; ir++) {
    wann->weight[ir]=wann1->weight[ir];
    (wann->rvec+ir)->x[0]=(wann1->rvec+ir)->x[0];
    (wann->rvec+ir)->x[1]=(wann1->rvec+ir)->x[1];
    (wann->rvec+ir)->x[2]=(wann1->rvec+ir)->x[2];

    for(io=0; io<wann->norb; io++) {
      for(jo=0; jo<wann->norb; jo++) {
        wann->ham[ir*wann->norb*wann->norb+io*wann->norb+jo] = 
           x*wann1->ham[ir*wann->norb*wann->norb+io*wann->norb+jo] +
           (1.0-x)*wann2->ham[ir*wann->norb*wann->norb+io*wann->norb+jo];
      }
    }
  }
}

main(int argc, char ** argv) {
  double x;  /* Percentage of first hamiltonian */
  double ef1, ef2;

  wanndata  wann1, wann2, wann;

  if (argc < 4) {
    printf("  Usage: %s <SEED1> <SEED2> percentage [Ef1] [Ef2].\n", argv[0]);
    exit(0);
  }

  ef1=0.0;
  ef2=0.0;

  if (argc==6) {
    sscanf(argv[4], " %lf", &ef1);
    sscanf(argv[5], " %lf", &ef2);
  }

  sscanf(argv[3], " %lf", &x);

  read_ham(&wann1, argv[1]);
  read_ham(&wann2, argv[2]);

  if ((wann1.norb!=wann2.norb) || (wann1.nrpt!=wann2.nrpt)) {
    printf("   Dimensions of wanniltonians doesn't match!!! Exiting...\n");
    finalize_wanndata(wann1);
    finalize_wanndata(wann2);
    exit(0);
  }
  else {
    shift_ham(&wann1, ef1);
    shift_ham(&wann2, ef2);
    mix_ham(&wann, &wann1, &wann2, x);
    write_ham(&wann);
  }

  finalize_wanndata(wann1);
  finalize_wanndata(wann2);
  finalize_wanndata(wann);

}
