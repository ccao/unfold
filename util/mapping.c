#include <stdio.h>
#include "vector.h"
#include "wannorb.h"
#include "mapping.h"

int setup_mapping(mapping * map, vector * shift, wannorb * target, wannorb * source, int ntgt, int nsrc) {
/*
 * This function finds mapping between two systems
 */
  int ii, jj;
  wannorb tmp;

  for(ii=0; ii<ntgt; ii++) {
    copy_wannorb(&tmp, target[ii]);

    if(shift!=NULL) {
      tmp.site=vector_add(tmp.site, (*shift));
    }

    for(jj=0; jj<nsrc; jj++) {
      if ( match_wannorb(&((map+ii)->rvec), tmp, source[jj]) ) {
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

int setup_symm_mapping(mapping * map, vector * symm, vector * shift, wannorb * wann, int nwann) {
/*
 * This function performs a symmetry operation, and finds out the mapping of orbitals
 *   between the systems before and after the symmetry operation.
 *
 * output:
 *   map  :  The mapping between orbitals
 * input:
 *   symm :  The symmetry operation. symm[0:2] defines R, symm[3] defines T
 *   shift:  Additional translation, if desired
 *   wann :  Wannier orbitals
 *   nwann:  Number of orbital sites
 */

  int ii, jj;
  int result;

  wannorb * tgt;

  tgt=(wannorb *) malloc(sizeof(wannorb)*nwann);

  for (ii=0; ii<nwann; ii++) {
    /*
     * Performs the symmetry operation
     *   Each operation consists of a 
     *    rotation R followed by a
     *    translation T.
     *   x'=R*x+T
     */
    init_wannorb(tgt+ii, NULL, wann[ii].l);
    symmop_wannorb(tgt+ii, wann[ii], (*shift), symm);
  }


  result=setup_mapping(map, shift, tgt, wann, nwann, nwann);

  for (ii=0; ii<nwann; ii++) {
    finalize_wannorb(tgt[ii]);
  }

  free(tgt);

  return result;
}
