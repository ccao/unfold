#ifdef VECTOR_H
#else
#define VECTOR_H
typedef struct __vector {
  double x;
  double y;
  double z;
} vector;

void init_vector(vector * v, double * x, double * y, double * z);
int equal(vector v1, vector v2);
double distance(vector v1, vector v2, vector * Tmat);
vector cross_product(vector v1, vector v2);
double dot_product(vector v1, vector v2);
double volume_product(vector v1, vector v2, vector v3);
vector vector_scale(double a, vector v);
vector vector_multiply(vector v1, int * n);
vector vector_add(vector v1, vector v2);
vector vector_sub(vector v1, vector v2);
int translate_match(vector *rv, vector x1, vector x2);
int isenclosed(vector v1, vector v2);
int find_vector(vector v, vector * list, int nlist);
vector rotate_vector(vector in, vector * symm);

#endif
