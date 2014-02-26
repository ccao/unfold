#include <stdio.h>
#include "vector.h"
#include "mapping.h"

int setup_mapping(mapping * map, vector * target, vector * source, vector * shift, int * info_tgt, int * info_src, int ntgt, int nsrc) {
  int ii, jj;
  vector vt;

  for(ii=0; ii<ntgt; ii++) {
    if(shift!=NULL) {
      vector_add(&vt, target[ii], (*shift));
    }
    else {
      vt.x[0]=(target+ii)->x[0];
      vt.x[1]=(target+ii)->x[1];
      vt.x[2]=(target+ii)->x[2];
    }

    for(jj=0; jj<nsrc; jj++) {
      if(info_tgt && info_src) {  
        /* if multiple source can match target,
           then you need additional info from info_tgt & info_src
             to determine the real mapping */
        if(info_tgt[ii]!=info_src[jj])
          continue;
      }
      if( translate_match(&((map+ii)->rvec), vt, source[jj]) ) {
        (map+ii)->nat=jj;
        break;
      }
    }
    if(jj==nsrc) {
      return -(ii+1);
    }
  }

  return 0;
}
