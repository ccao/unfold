#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "wanndata.h"

main(int argc, char ** argv) {
  double x;  /* Percentage of first hamiltonian */

  wanndata  ham1, ham2, ham;

  if (argc < 4) {
    printf("  Usage: %s <SEED1> <SEED2> percentage.\n", argv[0]);
    exit(0);
  }

  sscanf(argv[3], " %lf", &x);

  read_ham(&ham1, argv[1]);
  read_ham(&ham2, argv[2]);

  if ((ham1.norb!=ham2.norb) || (ham1.nrpt!=ham2.nrpt)) {
    printf("   Dimensions of hamiltonians doesn't match!!! Exiting...\n");
    finalize_wanndata(ham1);
    finalize_wanndata(ham2);
    exit(0);
  }
  else {
    mix_ham(&ham, &ham1, &ham2, x);
    write_ham(&ham);
  }

  finalize_wanndata(ham1);
  finalize_wanndata(ham2);
  finalize_wanndata(ham);

}
