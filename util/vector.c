#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "vector.h"

int equal(vector v1, vector v2) {
  if ((fabs(v1.x[0]-v2.x[0])<eps) &&
      (fabs(v1.x[1]-v2.x[1])<eps) &&
      (fabs(v1.x[2]-v2.x[2])<eps))
    return 1;
  else
    return 0;
}

double distance(vector v1, vector v2, vector * Tmat) {
  vector x1, x2;
  int i;
  if(Tmat!=NULL) {
    for(i=0; i<3; i++) {
      x1.x[i]=Tmat[0].x[i]*v1.x[0]+Tmat[1].x[i]*v1.x[1]+Tmat[2].x[i]*v1.x[2];
      x2.x[i]=Tmat[0].x[i]*v2.x[0]+Tmat[1].x[i]*v2.x[1]+Tmat[2].x[i]*v2.x[2];
    }
  }
  else {
    for(i=0; i<3; i++) {
      x1.x[i]=v1.x[i];
      x2.x[i]=v2.x[i];
    }
  }
    
  return sqrt((x1.x[0]-x2.x[0])*(x1.x[0]-x2.x[0])+
              (x1.x[1]-x2.x[1])*(x1.x[1]-x2.x[1])+
              (x1.x[2]-x2.x[2])*(x1.x[2]-x2.x[2]));
}

void cross_product(vector * res, vector v1, vector v2) {
  res->x[0]=v1.x[1]*v2.x[2];
  res->x[1]=v1.x[2]*v2.x[0];
  res->x[2]=v1.x[0]*v2.x[1];
}

double dot_product(vector v1, vector v2) {
  return v1.x[0]*v2.x[0]+v1.x[1]*v2.x[1]+v1.x[2]*v2.x[2];
}

double volume_product(vector v1, vector v2, vector v3) {
  double vol;
  vol=v1.x[0]*v2.x[1]*v3.x[2]+v1.x[1]*v2.x[2]*v3.x[0]+v1.x[2]*v2.x[0]*v3.x[1]-
     (v1.x[0]*v2.x[2]*v3.x[1]+v1.x[1]*v2.x[0]*v3.x[2]+v1.x[2]*v2.x[1]*v3.x[0]);
  return vol;
}
