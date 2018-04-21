#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "vector.h"
#include "poscar.h"

void real_to_reciprocal(vector * b, vector * a) {
  int i, j;
  double vol;
  vol=volume_product(a[0], a[1], a[2]);
  cross_product(&(b[0]), a[1], a[2]);
  cross_product(&(b[1]), a[2], a[0]);
  cross_product(&(b[2]), a[0], a[1]);
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      b[i].x[j]/=vol;
    }
  }
}

int main(int argc, char ** argv) {
  int nkpt, nen;
  char line[MAXLEN];
  char * p;
  FILE * fin;
  int ik, ii;
  double dos, kpos;
  poscar psc;
  vector a[3], b[3], k1, k2;

  init_poscar(&psc);

  fin=fopen("POSCAR.uc", "r");
  read_poscar_header(&psc, fin);
  fclose(fin);

  printf("#  Read unit cell lattices:\n");
  for(ii=0; ii<3; ii++) {
    printf("#   %12.8f%12.8f%12.8f\n", (psc.cell+ii)->x[0], (psc.cell+ii)->x[1], (psc.cell+ii)->x[2]);
  }

  real_to_reciprocal(b, psc.cell);

  printf("#  Reciprocal lattice:\n");
  for(ii=0; ii<3; ii++) {
    printf("#   %12.8f%12.8f%12.8f\n", b[ii].x[0], b[ii].x[1], b[ii].x[2]);
  }

  fin=fopen(argv[1], "r");
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d %d", &nkpt, &nen);

  kpos=0.0;
  for(ik=0; ik<nkpt; ik++) {
    printf("\n");
    fgets(line, MAXLEN, fin);
    sscanf(line, " %lf %lf %lf", &(k2.x[0]), &(k2.x[1]), &(k2.x[2]));
    if(ik>0) kpos+=distance(k2, k1, b);

    for(ii=0; ii<nen; ii++) {
      if(ii%10==0) {
        fgets(line, MAXLEN, fin);
        p=strtok(line, " ");
      }
      else {
        p=strtok(NULL, " ");
      }
      sscanf(p, " %lf", &dos);
      printf("%12.8f%12.8f%12.8f\n", kpos, (ii/1000.0), dos);
    }

    k1.x[0]=k2.x[0];
    k1.x[1]=k2.x[1];
    k1.x[2]=k2.x[2];
  }
  fclose(fin);
  finalize_poscar(psc);
  return 0;
}
