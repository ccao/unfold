#include <stdio.h>
#define MAXLEN 512
main () {
  int nkpt, nen;
  char line[MAXLEN];
  char * p;
  FILE * fin;
  int ik, ii;
  double dos;
  fin=fopen("unfold.log", "r");
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d %d", &nkpt, &nen);
  for(ik=0; ik<nkpt; ik++) {
    printf("\n");
    fgets(line, MAXLEN, fin);
    for(ii=0; ii<nen; ii++) {
      if(ii%20==0) {
        fgets(line, MAXLEN, fin);
        p=line;
      }
      sscanf(p, " %lf", &dos);
      printf("%5d%12.8f%12.8f\n", ik, (ii/1000.0), dos);
      p+=12;
    }
  }
  fclose(fin);
}
