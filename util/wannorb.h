#ifndef __WANNORB_H
#define __WANNORB_H

#include "constants.h"
#include "vector.h"

typedef struct __wannorb {
  vector site;	/* Orbital site */
  int l;	/* Maximal angular momentum for this orbital */
  vector * dir;	/* Orbital definition*/
/*
 * l & dir definition:
 * s : l=0; dir=NULL (s orbital do not have direction
 * p : l=1; px: dir=(1,0,0) etc...
 * sp2 : l=1; sp2-1: (1, 0, 0) etc...
 * d : l=2; dzx: dir=(1,0,0) and (0,0,1) ...
 */
} wannorb;

void init_wannorb(wannorb * wann, vector * v, int l);
void copy_wannorb(wannorb * target, wannorb source);
void finalize_wannorb(wannorb wann);
int reappear(int * list, int length);
void symmop_wannorb(wannorb * out, wannorb in, vector shift, vector * symm);
int match_wannorb(vector * vt, wannorb orb1, wannorb orb2);

#endif   /* __WANNORB_H*/
