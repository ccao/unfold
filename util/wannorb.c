#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "constants.h"
#include "vector.h"

int parse_line(char * line, vector * site, int * lm, vector * axis) {
  char * pch, *ppc;
  int itok;
  int ii;
  itok=0;
  pch=strtok(line, ":");
  while (pch!=NULL) {
    switch(itok) {
      case 0:
        sscanf(pch, " %lf, %lf, %lf", &(site->x), &(site->y), &(site->z));
        break;
      case 1:
        identify(pch);
        break;
      default:
        ppc=strtok(pch, "=");
        switch(ppc[0]) {
          case 'x':
            ii=0;
            break;
          case 'z':
            ii=2;
            break;
          default:
            printf("!!! Wrong input while specifying the directions of Wannier orbitals.\n");
        }
        ppc=strtok(NULL, "=");
        sscanf(ppc, " %lf, %lf, %lf", &(axis[ii].x), &(axis[ii].y), &(axis[ii].z));
    }
  }
  else {
      
}

void setup_wannorb_dir(wannorb * wann, vector * axis, int ii, int lm) {
}

int read_wannorb(wannorb * wann, FILE * fin) {
/*
 * This function reads Wannier orbital definitions from a file
 * If wann is NULL, then determines the actual wannier orbital
 * Otherwise assumes wann is allocated and populates the 
 * wannier definition
 * INPUT:
 *   wann : if NULL returns number of wannier orbitals
 *          otherwise is an output
 *   fin  : input file handle
 * OUTPUT:
 *   wann : actual wannier orbitals
 * RETURNS the number of wannier orbitals
 */
  int ii, jj, nn;
  int iw, nw;
  int lm;
  /*
   * lm: specification of orbitals (index)
   *   a three digits integer number
   *   whose first digit determines l
   *   and second digit determines the type of hybridization
   *   and the last digit determines the actual orbital(s)
   * A List:
   *   lm | orbital(s)
   *   000 |  s
   *   100 |  p 
   *   101-103 |  pz, px, py
   *   110 |  sp
   *   111-112 |  sp-1, sp-2
   *   120 |  sp2
   *   121-123 |  sp2-1, sp2-2, sp2-3
   *   130 |  sp3
   *   131-134 |  sp3-1, sp3-2, sp3-3, sp3-4
   *   200 |  d
   *   201-205 | dz2, dzx, zy ...
   *   210 |  t2g/eg
   *   211-213 | t2g (dzy, dxz, dxy)
   *   214-215 | eg (...)
   *   300 |  f...
   */
  int ll;
  char line[MAXLEN];
  vector axis[3];
  vector site;

  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &nn);
  iwann=0;
  for (ii=0; ii<nn; ii++) {
    fgets(line, MAXLEN, fin);
    nw=parse_line(line, &site, &lm, axis);
    ll=lm/100;
    if (wann) {
      for (jj=0; jj<nw; jj++) {
        init_wannorb(wann+iwann, &site, ll);
        set_wannorb_dir(wann+iwann, axis, jj, lm);
        iwann++;
      }
    }
    else {
      iwann+=nw;
    }
  }

  return iwann;
}
  

void init_wannorb(wannorb * wann, vector * v, int l) {
  if (v!=NULL) {
    wann->site=(*v);
  }
  wann->l=l;
  if (l>0) {
    wann->dir=(vector *)malloc(sizeof(vector)*l);
  }
  else {
    wann->dir=NULL;
  }
}

void copy_wannorb(wannorb * target, wannorb source) {
  init_wann(target, source.site, source.l);
  memcpy(target->dir, source.dir, sizeof(vector)*source.l);
}

void finalize_wannorb(wannorb wann) {
  free(wann.dir);
}

int reappear(int * list, int length) {
  int i;
  for (i=0; i<length-1; i++)
    if (list[i]==list[length]) return -1;
  return 0;
}

void symmop_wannorb(wannorb * out, wannorb in, vector shift, vector * symm) {
/*
 * This function performs symmetry operator on a specific wannier orbital
 * OUTPUT:
 *   out: output wannier orbital
 * INPUT :
 *   in : input wannier orbital
 *   shift : additional shift vector (usually integer lattice constants)
 *   symm : symmetry operator
 */
  int ii;
  vector vt;
  vt=vector_add(in.site, shift);
  out->site=rotate_vector(vt, symm);
  out->l=in.l;
  for (ii=0; ii<in.l; ii++) {
    out->dir[ii]=rotate_vector(in.dir[ii], symm);
  }
}

int match_wannorb(vector * vt, wannorb orb1, wannorb orb2) {
  int i, j;
  int * m;

  if (orb1.l!=orb2.l) {
    return 0;
  }
  if (!translate_match(vt, orb1.site, orb2.site)) {
    return 0;
  }
  if (orb1.l>0) {
    m=(int *) malloc(sizeof(int)*orb1.l);
    for (i=0; i<orb1.l; i++) {
      for (m[i]=0; m[i]<orb2.l; m[i]++) {
        if (reappear(m, i)) continue;

        if (equal(orb1.dir[i], orb2.dir[m[i]]) break;
      }
      if (m[i]==orb2.l) {
        free(m);
        return 0;
      }
    }
  }
  free(m);
  return 1;
}
