//3d math libary
typedef double* V3;

static inline void v3_add(V3 c, V3 a, V3 b){
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

static inline void v3_sub(V3 c, V3 a, V3 b){
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

static inline void v3_scale(V3 c, V3 a, double b){
  c[0] = a[0] * b;
  c[1] = a[1] * b;
  c[2] = a[2] * b;
}

static inline double v3_dot(V3 a, V3 b){
  double c;
  c = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
  return c;
}

static inline void v3_cross(V3 c, V3 a, V3 b){
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void v3_assign(V3 c, double val1, double val2, double val3){
  c[0] = val1;
  c[1] = val2;
  c[2] = val3;
}

static inline double v3_distance(V3 a, V3 b){
  double c;
  c = (pow(a[0] - b[0], 2)) + (pow(a[1] - b[1], 2)) + (pow(a[2] - b[2], 2));
  c = sqrt(c);
  return c;
}

static inline void v3_normalize(V3 c, V3 a){
  double mag = (a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]);
  mag = sqrt(mag);
  c[0] = a[0]/mag;
  c[1] = a[1]/mag;
  c[2] = a[2]/mag;
}
