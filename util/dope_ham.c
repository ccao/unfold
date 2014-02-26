#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "vector.h"
#include "wanndata.h"
#include "poscar.h"
#include "mapping.h"

void setup_orblst(vector * orblst, int * orbsub, poscar psc) {
  int ii, jj, kk, iorb, iatm;
  int nsp;
  int * norb_per_sp;
  char line[MAXLEN];
  char * p;
  FILE * fin;

  fin=fopen("orbital.def", "r");
  fgets(line, MAXLEN, fin);
  fgets(line, MAXLEN, fin);
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &nsp);
  if (nsp!=psc.nsp) {
    printf("!!! WARNING: nsp in orbital definition incompatible with POSCAR.\n");
  }
  norb_per_sp=(int *) malloc(sizeof(int)*nsp);
  fgets(line, MAXLEN, fin);
  p=strtok(line, " ");
  for(ii=0; ii<nsp; ii++) {
    sscanf(p, "%d", norb_per_sp+ii);
    p=strtok(NULL, " ");
  }

  fclose(fin);

  iorb=0;
  iatm=0;
  for(ii=0; ii<nsp; ii++) {
    for(jj=0; jj<psc.nat_per_sp[ii]; jj++) {
      for(kk=0; kk<norb_per_sp[ii]; kk++) {
        (orblst+iorb)->x[0]=(psc.tau+iatm)->x[0];
        (orblst+iorb)->x[1]=(psc.tau+iatm)->x[1];
        (orblst+iorb)->x[2]=(psc.tau+iatm)->x[2];
        orbsub[iorb]=kk;
        iorb++;
      }
      iatm++;
    }
  }

  free(norb_per_sp);
}

void apply_doping(wanndata ham0, wanndata dham, vector * orblst, int * orbsub, vector shft){
  mapping * orbmap;
  int nrpt, norb;

  int irpt, iorb, jorb;
  /* indices in new (translated) hamiltonian) */
  int iirpt, iiorb, jjorb;
  /* indices in old (original) hamiltonian), obtained by orbmap */
  int ierr;
  vector rv;

  nrpt=ham0.nrpt;
  norb=ham0.norb;

  orbmap=(mapping *)malloc(sizeof(mapping)*norb);
  ierr=setup_mapping(orbmap, orblst, orblst, &shft, orbsub, orbsub, norb, norb);

  if (ierr) {
    printf("!!!! ERROR: Cannot map orbital # %4d.\n", ierr);
    printf("!!!!   Please check if dopants are equivalent.\n");
    exit(0);
  }

  for(irpt=0; irpt<nrpt; irpt++) {
    for(iorb=0; iorb<norb; iorb++) {
      iiorb=(orbmap+iorb)->nat;
      for(jorb=0; jorb<norb; jorb++) {
        jjorb=(orbmap+jorb)->nat;
        vector_sub(&rv, (orbmap+iorb)->rvec, (orbmap+jorb)->rvec);
        vector_add(&rv, rv, ham0.rvec[irpt]);
        iirpt=locate_rpt(&ham0, rv);

        if (iirpt>=0) {
          ham0.ham[irpt*norb*norb+iorb*norb+jorb]+=dham.ham[iirpt*norb*norb+iiorb*norb+jjorb];
        }
      }
    }
  }

  free(orbmap);

}

void read_input(int * ndpnt, vector ** shftvec, poscar psc) {
  int origin, target;
  char line[MAXLEN];
  char *p;
  int i;
  FILE * fin;

  fin=fopen("dopants.in", "r");
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", ndpnt);
  (*shftvec)=(vector *) malloc(sizeof(vector)*(*ndpnt));

  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &origin);

  fgets(line, MAXLEN, fin);
  p=strtok(line, " ");
  for(i=0; i<(*ndpnt); i++) {
    sscanf(p, "%d", &target);
    ((*shftvec)+i)->x[0]=(psc.tau+target)->x[0]-(psc.tau+origin)->x[0];
    ((*shftvec)+i)->x[1]=(psc.tau+target)->x[1]-(psc.tau+origin)->x[1];
    ((*shftvec)+i)->x[2]=(psc.tau+target)->x[2]-(psc.tau+origin)->x[2];
    p=strtok(NULL, " ");
  }
}

int main(int argc, char ** argv) {
  int ndpnt;

  vector * orblst;
  int * orbsub;
  /* orbital vector list */
  vector * shftvec;
  /* doping vector list */

  wanndata ham0, dham;
  /* ham0: hamiltonian for the prestine system (in) 
           hamiltonian for the doped system (out)
     dham: hamiltonian difference for single dopant */

  poscar psc;
  /* psc0: crystal structure for the system */

  int idp;

  if (argc<3) {
    printf("Usage: %s <PRESTINE_SEED> <DELTAH_SEED>\n", argv[0]);
    printf("  Input file (dopants.in) details:\n");
    printf("  line #1: # of dopants.\n");
    printf("  line #2: the original dopant (in DELTAH_SEED)\n");
    printf("  line #3: new dopant positions (indices in POSCAR)\n");
    exit(0);
  }

  read_poscar(&psc, "POSCAR.sc");
  read_input(&ndpnt, &shftvec, psc);

  read_ham(&ham0, argv[1]);
  read_ham(&dham, argv[2]);

  orblst=(vector *) malloc(sizeof(vector)*ham0.norb);
  orbsub=(int *) malloc(sizeof(int)*ham0.norb);
  setup_orblst(orblst, orbsub, psc);

  for(idp=0; idp<ndpnt; idp++) {
    apply_doping(ham0, dham, orblst, orbsub, shftvec[idp]);
  }

  write_ham(&ham0);

  free(orblst);
  free(shftvec);

  finalize_wanndata(ham0);
  finalize_wanndata(dham);

}
