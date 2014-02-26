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
void vector_add(vector * vr, vector v1, vector v2);
void vector_sub(vector * vr, vector v1, vector v2);
int translate_match(vector *rv, vector x1, vector x2);

#endif
