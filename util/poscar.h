#include "constants.h"
#include "vector.h"

typedef struct __poscar{
  int nsp;
  int nat;
  int * nat_per_sp;
  vector cell[3];
  vector * tau;
} poscar;

void read_poscar_header(poscar * psc, FILE * fin);
void read_poscar(poscar * psc, char * fn);
void finalize_poscar(poscar psc);

