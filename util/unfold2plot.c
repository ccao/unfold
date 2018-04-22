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
  FILE * fin, * fout;
  double elow, ehigh;
  double * kpos, * dos;
  int ik, ii;
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
  sscanf(line, " %d %d %lf %lf", &nkpt, &nen, &elow, &ehigh);

  kpos=(double *)malloc(sizeof(double)*nkpt);
  dos=(double *)malloc(sizeof(double)*nkpt*nen);

  for(ik=0; ik<nkpt; ik++) {
    fgets(line, MAXLEN, fin);
    sscanf(line, " %lf %lf %lf", &(k2.x[0]), &(k2.x[1]), &(k2.x[2]));
    if(ik>0) {
      kpos[ik]=distance(k2, k1, b)+kpos[ik-1];
    }
    else {
      kpos[ik]=0.0;
    }

    for(ii=0; ii<nen; ii++) {
      if(ii%10==0) {
        fgets(line, MAXLEN, fin);
        p=strtok(line, " ");
      }
      else {
        p=strtok(NULL, " ");
      }
      sscanf(p, " %lf", dos+ik*nen+ii);
    }

    k1.x[0]=k2.x[0];
    k1.x[1]=k2.x[1];
    k1.x[2]=k2.x[2];
  }

  fout=fopen("plot.dat", "w");
  fprintf(fout, "%10d%10d", nkpt, nen);
  for(ik=0; ik<nkpt; ik++) {
    if (ik%10==0) fprintf(fout,"\n");
    fprintf(fout, "%14.9f ", kpos[ik]);
  }
  for(ii=0; ii<nen; ii++) {
    if (ii%10==0) fprintf(fout,"\n");
    fprintf(fout, "%14.9f ", elow+ii*(ehigh-elow)/(nen-1));
  }
  for(ik=0; ik<nkpt; ik++) {
    for(ii=0; ii<nen; ii++) {
      if (ii%10==0) fprintf(fout,"\n");
      fprintf(fout, "%14.9f ", dos[ik*nen+ii]);
    }
  }
  fclose(fout);

  free(kpos);
  free(dos);

  fclose(fin);
  finalize_poscar(psc);
  return 0;
}
