#include "wanndata.h"
#include "vector.h"
int read_orbdef(int ** norb_per_sp, char * fuc, char * fsc, FILE * fin);
void rotate_ham(wanndata * out, wanndata * in, vector * symm, vector * site, int * info);
