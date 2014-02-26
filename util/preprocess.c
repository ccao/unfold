#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "constants.h"
#include "vector.h"
#include "poscar.h"
#include "mapping.h"

void orbital_index(int *** orb_idx, poscar uc, int * norb_sp) {
  int isp, ii, iat, iorb, iiorb;

  (*orb_idx)=(int **) malloc(sizeof(int *)*uc.nat);

  iat=0;
  iiorb=1;
  for(isp=0; isp<uc.nsp; isp++) {
    for(ii=0; ii<uc.nat_per_sp[isp]; ii++) {
      if(norb_sp[isp]>0) {
        (*orb_idx)[iat]=(int *) malloc(sizeof(int)*norb_sp[isp]);
        for(iorb=0; iorb<norb_sp[isp]; iorb++) {
          (*orb_idx)[iat][iorb]=iiorb++;
        }
      }
      else {
        (*orb_idx)[iat]=(int *) malloc(sizeof(int));
        (*orb_idx)[iat][0]=-1;
      }
      iat++;
    }
  }
}

int search_transform(int * t, vector target, vector * source) {
  int jj;
  vector tc;

  for(t[0]=-MAXN; t[0]<MAXN+1; t[0]++)
   for(t[1]=-MAXN; t[1]<MAXN+1; t[1]++)
    for(t[2]=-MAXN; t[2]<MAXN+1; t[2]++) {
      for(jj=0; jj<3; jj++)
        tc.x[jj]=t[0]*source[0].x[jj]+t[1]*source[1].x[jj]+t[2]*source[2].x[jj];
      if (equal(target, tc))
        return 0;
    }

  return -1;
}

void retrieve_scell(double scell[3][3], vector * sc, vector * uc) {
  int ii, jj;
  int t[3];

  for(ii=0; ii<3; ii++) {
    if(search_transform(t, sc[ii], uc)==0) {
      for(jj=0; jj<3; jj++) scell[ii][jj]=t[jj];
    }
    else {
      printf("  !!! ERROR: cannot find transformation matrix for lattice %d\n", ii);
    }
  }
}

void transform_structure(poscar * sc, double scell[3][3]) {
  int ii, jj;
  vector tv;
  for(ii=0; ii<sc->nat; ii++) {
    for(jj=0; jj<3; jj++)
      tv.x[jj]=scell[0][jj]*(sc->tau+ii)->x[0]+scell[1][jj]*(sc->tau+ii)->x[1]+scell[2][jj]*(sc->tau+ii)->x[2];
    for(jj=0; jj<3; jj++)
      (sc->tau+ii)->x[jj]=tv.x[jj];
  }
}

void read_input(char * fsc, char * fuc, int * nsp, int ** norb_sp, FILE * fin) {
  char line[MAXLEN];
  char *p;
  int ii;

  fgets(line, MAXLEN, fin);
  sscanf(line, " %s", fuc);
  fgets(line, MAXLEN, fin);
  sscanf(line, " %s", fsc);
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", nsp);
  (*norb_sp)=(int *) malloc(sizeof(int)*(*nsp));
  fgets(line, MAXLEN, fin);
  p=strtok(line, " ");
  for (ii=0; ii<(*nsp); ii++) {
    sscanf(p, "%d", (*norb_sp)+ii);
    p=strtok(NULL, " ");
  }
}

void output_mapping(double scell[3][3], mapping *map, poscar uc, poscar sc, int * norb_sp, int ** orb_idx_at) {
  int ii, jj, kk, iat;
  int norb_sc, norb_uc;
  for(ii=0; ii<3; ii++) {
    printf("%8d%8d%8d\n", (int)(scell[ii][0]), (int)(scell[ii][1]), (int)(scell[ii][2]));
  }
  norb_sc=0;
  norb_uc=0;
  for(ii=0; ii<sc.nsp; ii++) {
    norb_sc+=sc.nat_per_sp[ii]*norb_sp[ii];
    norb_uc+=uc.nat_per_sp[ii]*norb_sp[ii];
  }
  printf("%10d%10d\n", norb_uc, norb_sc);

  iat=0;
  for(ii=0; ii<sc.nsp; ii++) {
    for(jj=0; jj<sc.nat_per_sp[ii]; jj++) {
      for(kk=0; kk<norb_sp[ii]; kk++) {
        printf("%5d  %5d%5d%5d\n", orb_idx_at[(map+iat)->nat][kk], (int)((map+iat)->rvec).x[0], (int)((map+iat)->rvec).x[1], (int)((map+iat)->rvec).x[2]);
      }
      iat++;
    }
  }
}
        

int main(int argc, char ** argv) {
  FILE * fin;
  poscar sc, uc;
  mapping * map;
  double scell[3][3];
  int nsp, ii;
  int * norb_sp, ** orb_idx_at;
  char fsc[20], fuc[20];

  if(argc<2) {
    printf("Usage: %s <INPUT>\n", argv[0]);
    exit(0);
  }
  fin=fopen(argv[1], "r");
  read_input(fsc, fuc, &nsp, &norb_sp, fin);
  fclose(fin);

  read_poscar(&sc, fsc);
  read_poscar(&uc, fuc);

  if(nsp!=sc.nsp || nsp!=uc.nsp) {
    printf("  !!! WARNING: Number of species doesn't match.\n");
  }

  orbital_index(&orb_idx_at, uc, norb_sp);

  retrieve_scell(scell, sc.cell, uc.cell);

  transform_structure(&sc, scell);

  map=(mapping *) malloc(sizeof(mapping)*sc.nat);

  setup_mapping(map, sc.tau, uc.tau, NULL, NULL, NULL, sc.nat, uc.nat);

  output_mapping(scell, map, uc, sc, norb_sp, orb_idx_at); 

  finalize_poscar(sc);
  finalize_poscar(uc);
  free(map);
  free(norb_sp);
  for(ii=0; ii<uc.nat; ii++)
    free(orb_idx_at[ii]);
  free(orb_idx_at);

  return 0;
}
