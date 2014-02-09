#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXLEN 512
#define MAXN   10
#define eps   1E-5

typedef struct __vector{
  double x[3];
} vector;

typedef struct __mapping{
  int nat;
  int rvec[3];
} mapping;

typedef struct __poscar{
  int nsp;
  int nat;
  int * nat_per_sp;
  vector cell[3];
  vector * tau;
} poscar;

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

int equal(vector v1, vector v2) {
  if ((fabs(v1.x[0]-v2.x[0])<eps) &&
      (fabs(v1.x[1]-v2.x[1])<eps) &&
      (fabs(v1.x[2]-v2.x[2])<eps))
    return 1;
  else
    return 0;
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

int match(int *rv, vector x1, vector x2) {
  int ii;
  double dx[3];
  for(ii=0; ii<3; ii++) {
    rv[ii]=(int)(rint(x1.x[ii]-x2.x[ii]));
    dx[ii]=x1.x[ii]-x2.x[ii]-rv[ii];
  }
  if( (fabs(dx[0])<eps) &&
      (fabs(dx[1])<eps) &&
      (fabs(dx[2])<eps) )
    return 1;
  else
    return 0;
}

void retrieve_mapping(mapping * map, poscar sc, poscar uc, double scell[3][3]) {
  int ii, jj;
  int rv[3];
  vector tv;

  for(ii=0; ii<sc.nat; ii++) {
    for(jj=0; jj<3; jj++) 
      tv.x[jj]=scell[0][jj]*(sc.tau+ii)->x[0]+scell[1][jj]*(sc.tau+ii)->x[1]+scell[2][jj]*(sc.tau+ii)->x[2];

    for(jj=0; jj<uc.nat; jj++) {
      if( match(rv, tv, uc.tau[jj]) ){
        (map+ii)->nat=jj;
        (map+ii)->rvec[0]=rv[0];
        (map+ii)->rvec[1]=rv[1];
        (map+ii)->rvec[2]=rv[2];
        break;
      }
    }
    if (jj==uc.nat) {
      printf("  !!! ERROR: cannot find match for atom %d in super cell\n", ii);
      exit(0);
    }
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
  p=strtok(line, " \t");
  for (ii=0; ii<(*nsp); ii++) {
    sscanf(p, "%d", (*norb_sp)+ii);
    p=strtok(NULL, " \t");
  }
}

void read_poscar(poscar * psc, char * fn) {
  FILE * fin;
  char line[MAXLEN];
  char * p;
  double lat;
  int ii;

  fin=fopen(fn, "r");
  fgets(line, MAXLEN, fin);
  fgets(line, MAXLEN, fin);
  sscanf(line, " %lf", &lat);

  for(ii=0; ii<3; ii++) {
    fgets(line, MAXLEN, fin);
    sscanf(line, " %lf %lf %lf", (psc->cell+ii)->x, (psc->cell+ii)->x+1, (psc->cell+ii)->x+2);
    (psc->cell[ii]).x[0]*=lat;
    (psc->cell[ii]).x[1]*=lat;
    (psc->cell[ii]).x[2]*=lat;
  }

  fgets(line, MAXLEN, fin);
  ii=0;
  p=strtok(line, " \t");
  while(p!=NULL) {
    ii++;
    p=strtok(NULL, " \t");
  }
  psc->nsp=ii;

  psc->nat_per_sp=(int *) malloc(sizeof(int)*psc->nsp);

  fgets(line, MAXLEN, fin);

  psc->nat=0;
  p=strtok(line, " \t");
  for(ii=0; ii<psc->nsp; ii++) {
    sscanf(p, " %d", psc->nat_per_sp+ii);
    psc->nat+=psc->nat_per_sp[ii];
    p=strtok(NULL, " \t");
  }
  psc->tau=(vector *)malloc(sizeof(vector)*psc->nat);
  fgets(line, MAXLEN, fin);
  for(ii=0; ii<psc->nat; ii++) {
    fgets(line, MAXLEN, fin);
    sscanf(line, " %lf %lf %lf ", (psc->tau+ii)->x, (psc->tau+ii)->x+1, (psc->tau+ii)->x+2);
  }
  fclose(fin);
}

void finalize(poscar psc) {
  free(psc.nat_per_sp);
  free(psc.tau);
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
        printf("%5d  %5d%5d%5d\n", orb_idx_at[(map+iat)->nat][kk], (map+iat)->rvec[0], (map+iat)->rvec[1], (map+iat)->rvec[2]);
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

  map=(mapping *) malloc(sizeof(mapping)*sc.nat);

  retrieve_mapping(map, sc, uc, scell);

  output_mapping(scell, map, uc, sc, norb_sp, orb_idx_at); 

  finalize(sc);
  finalize(uc);
  free(map);
  free(norb_sp);
  for(ii=0; ii<uc.nat; ii++)
    free(orb_idx_at[ii]);
  free(orb_idx_at);

  return 0;
}
