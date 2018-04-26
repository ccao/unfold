#include <stdio.h>
#include "vector.h"
#include "mapping.h"

int setup_mapping(mapping * map, vector * target, vector * source, vector * shift, int * info_tgt, int * info_src, int ntgt, int nsrc) {
/*
 * This function finds mapping between two systems
 */
  int ii, jj;
  vector vt;

  for(ii=0; ii<ntgt; ii++) {
    if(shift!=NULL) {
      vector_add(&vt, target[ii], (*shift));
    }
    else {
      vt=target[ii];
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

int setup_symm_mapping(mapping * map, vector * symm, vector * shift, vector * site, int * info, int nsite) {
/*
 * This function performs a symmetry operation, and finds out the mapping of orbitals
 *   between the systems before and after the symmetry operation.
 *
 * output:
 *   map  :  The mapping between orbitals
 * input:
 *   symm :  The symmetry operation. symm[0:2] defines R, symm[3] defines T
 *   shift:  Additional translation, if desired
 *   site :  orbital sites
 *   info :  additional orbital labels
 *   nsite:  Number of orbital sites
 */

  int ii, jj;
  int result;

  vector * tgt;

  tgt=(vector *) malloc(sizeof(vector)*nsite);

  for (ii=0; ii<nsite; ii++) {
    /*
     * Performs the symmetry operation
     *   Each operation consists of a 
     *    rotation R followed by a
     *    translation T.
     *   x'=R*x+T
     */
    for (jj=0; jj<3; jj++ )
      tgt[ii].x[jj]=dot_product(symm[jj], site[ii])+symm[3].x[jj];
  }


  result=setup_mapping(map, tgt, site, shift, info, info, nsite, nsite);

  free(tgt);

  return result;
}
