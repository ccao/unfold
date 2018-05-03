#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "constants.h"
#include "vector.h"
#include "wanndata.h"
#include "wannorb.h"
#include "mapping.h"

void rotate_ham(wanndata * out, wanndata * in, vector * symm, wannorb * wann) {
/*
 * This function generates a new Hamiltonian by performing symmetry operation
 *
 * output:
 *   out   : output Hamiltonian
 * input:
 *   in    : input Hamiltonian
 *   symm  : symmetry operator (in lattice vector unit)
 *   wann  : Wannier orbitals
 */
  mapping * map0;    /* Mapping of orbitals at (0, 0, 0) cell */
  mapping * mapr;    /* Mapping of oribtals at (Rx, Ry, Rz) cell */

  int ir, ii, jj, nn;
  int iir;
  
  vector v;

  nn=in->norb;
  map0=(mapping *)malloc(sizeof(mapping)*nn);
  mapr=(mapping *)malloc(sizeof(mapping)*nn);

  /* The symmetry operation should not change Rvec and weights */
  memcpy(out->rvec, in->rvec, sizeof(vector)*in->nrpt);
  memcpy(out->weight, in->weight, sizeof(int)*in->nrpt);

  iir=setup_symm_mapping(map0, symm, NULL, wann, nn);
  if (iir<0) {
    printf("ERR: Cannot find mapping for orbital %5d\n", (-iir));
    printf("ERR: shift: (0, 0, 0)\n");
  }

  for (ir=0; ir<in->nrpt; ir++) {
    iir=setup_symm_mapping(mapr, symm, in->rvec+ir, wann, nn);
    if (iir<0) {
      printf("ERR: Cannot find mapping for orbital %5d\n", -iir);
       printf("ERR:  shift: (%5d, %5d, %5d)\n", (int)in->rvec[ir].x, (int)in->rvec[ir].y, (int)in->rvec[ir].z);
    }
    for (ii=0; ii<nn; ii++) {
      for (jj=0; jj<nn; jj++) {
        /* We need the relative distance between orbitals */
        v=vector_sub(mapr[jj].rvec, map0[ii].rvec);
        iir=find_vector(v, in->rvec, in->nrpt);
        if (iir>=0) {
          out->ham[ir*nn*nn+ii*nn+jj]=in->ham[iir*nn*nn+map0[ii].nat*nn+mapr[jj].nat];
        }
        else {
          printf(" # Cannot find mapping for (%5d, %5d, %5d)\n", (int)v.x, (int)v.y, (int)v.z);
          printf(" # Setting it to be 0\n");
          out->ham[ir*nn*nn+ii*nn+jj]=0.0;
        }
      }
    }
  }
  free(map0);
  free(mapr);
}
