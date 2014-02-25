#ifdef VECTOR_H
#else
#define VECTOR_H
typedef struct __vector {
  double x[3];
} vector;

int equal(vector v1, vector v2);
double distance(vector v1, vector v2, vector * Tmat);
void cross_product(vector * res, vector v1, vector v2);
double dot_product(vector v1, vector v2);
double volume_product(vector v1, vector v2, vector v3);

#endif
