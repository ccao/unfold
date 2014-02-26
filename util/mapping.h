#ifndef _MAPPING_H
#define _MAPPING_H

#include "vector.h"

typedef struct __mapping{
  int nat;
  vector rvec;
} mapping;

int setup_mapping(mapping * map, vector * target, vector * source, vector * shift, int * info_tgt, int * info_src, int ntgt, int nsrc);

#else
#endif
