#ifndef _MAPPING_H
#define _MAPPING_H

#include "wannorb.h"
#include "vector.h"

typedef struct __mapping{
  int nat;
  vector rvec;
} mapping;

int setup_mapping(mapping * map, vector * shift, wannorb * target, wannorb * source, int ntgt, int nsrc);
int setup_symm_mapping(mapping * map, vector * symm, vector * shift, wannorb * wann, int nwann);
#else
#endif
