#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "vector.h"

void init_vector(vector * v, double x, double y, double z) {
  v->x=x;
  v->y=y;
  v->z=z;
}

int translate_match(vector * rv, vector x1, vector x2) {
  int ii;
  vector tmp;
  rv->x=(int)(rint(x1.x-x2.x);
  rv->y=(int)(rint(x1.y-x2.y);
  rv->z=(int)(rint(x1.z-x2.z);
  tmp=vector_add(x2, (*rv));
  if (distance(x1, tmp)<eps)
    return 1;
  else
    return 0;
}

vector vector_scale(double a, vector v) {
  vector r;
  r.x=a*v.x;
  r.y=a*v.y;
  r.z=a*v.z;

  return r;
}

vector vector_multiply(vector v1, int * n) {
  vector r;
  r.x=v1.x*n[0];
  r.y=v1.y*n[1];
  r.z=v1.z*n[2];

  return r;
}

vector vector_sub(vector v1, vector v2) {
  vector vr;
  vr.x=v1.x-v2.x;
  vr.y=v1.y-v2.y;
  vr.z=v1.z-v2.z;
  return vr;
}

vector vector_add(vector v1, vector v2) {
  vector vr;
  vr.x=v1.x+v2.x;
  vr.y=v1.y+v2.y;
  vr.z=v1.z+v2.z;
  return vr;
}

int equal(vector v1, vector v2) {
  if ((fabs(v1.x[0]-v2.x[0])<eps) &&
      (fabs(v1.x[1]-v2.x[1])<eps) &&
      (fabs(v1.x[2]-v2.x[2])<eps))
    return 1;
  else
    return 0;
}

vector matrix_dot(vector * Tmat, vector v) {
  vector r;
  r.x=v.x*Tmat[0].x+v.y*Tmat[1].x+v.z*Tmat[2].x;
  r.y=v.x*Tmat[0].y+v.y*Tmat[1].y+v.z*Tmat[2].y;
  r.z=v.x*Tmat[0].z+v.y*Tmat[1].z+v.z*Tmat[2].z;
  return r;
}

double distance(vector v1, vector v2, vector * Tmat) {
  double r;
  vector x1, x2;
  if (Tmat!=NULL) {
    x1=matrix_dot(Tmat, v1);
    x2=matrix_dot(Tmat, v2);
  }
  else {
    x1=v1;
    x2=v2;
  }

  r=sqrt((x1.x-x2.x)*(x1.x-x2.x)+(x1.y-x2.y)*(x1.y-x2.y)+(x1.z-x2.z)*(x1.z-x2.z));
  return r;
}

vector cross_product(vector v1, vector v2) {
  vector r;
  r.x=v1.y*v2.z-v1.z*v2.y;
  r.y=v1.z*v2.x-v1.x*v2.z;
  r.z=v1.x*v2.y-v1.y*v2.x;

  return r;
}

double dot_product(vector v1, vector v2) {
  return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

double volume_product(vector v1, vector v2, vector v3) {
  double vol;
  vol=v1.x*v2.y*v3.z+v1.y*v2.z*v3.x+v1.z*v2.x*v3.y-
     (v1.x*v2.z*v3.y+v1.y*v2.x*v3.z+v1.z*v2.y*v3.x);
  return vol;
}

int isenclosed(vector v1, vector v2) {
  return ( fabs(v1.x[0])<=v2.x[0] &&
           fabs(v1.x[1])<=v2.x[1] &&
           fabs(v1.x[2])<=v2.x[2] );
}

int find_vector(vector v, vector * list, int nlist) {
  int ii;
  for (ii=0; ii<nlist; ii++) {
    if (equal(v, list[ii]))
      return ii;
  }
  return -1;
}

vector rotate_vector(vector in, vector * symm) {
  vector out;
  out.x=dot_product(symm[0], in);
  out.y=dot_product(symm[1], in);
  out.z=dot_product(symm[2], in);

  return out;
}
