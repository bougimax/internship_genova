#include "quality_measure.h"

double energy_dirichlet(const pointType *p1, const pointType *p2,
                        const pointType *p3, const pointType *p4) {

  const vector3d v1(p1), v2(p2), v3(p3), v4(p4);
  vector3d j1(v1.c[0] + v3.c[0] - v2.c[0] - v4.c[0],
              v1.c[1] + v3.c[1] - v2.c[1] - v4.c[1],
              v1.c[2] + v3.c[2] - v2.c[2] - v4.c[2]);
  vector3d j2(v1.c[0] - v3.c[0] - v2.c[0] + v4.c[0],
              v1.c[1] - v3.c[1] - v2.c[1] + v4.c[1],
              v1.c[2] - v3.c[2] - v2.c[2] + v4.c[2]);
  vector3d j3((-v1.c[0]) + v3.c[0] - v2.c[0] + v4.c[0],
              (-v1.c[1]) + v3.c[1] - v2.c[1] + v4.c[1],
              (-v1.c[2]) + v3.c[2] - v2.c[2] + v4.c[2]);

  j1 *= 0.5;
  j2 *= 0.5;
  j3 *= 0.5;

  const double num = (j1 * j1) + (j2 * j2) + (j3 * j3);
  double det = j1.tripleProd(j3, j2);
  if (det <= 0) {
    return DBL_MAX;
  }

  return num / pow(det, (2.0 / 3.0));
}

double energy_dirichlet_big_float(const pointType *p1, const pointType *p2,
                        const pointType *p3, const pointType *p4) {

  const vec3d_bf v1(p1), v2(p2), v3(p3), v4(p4);
  vec3d_bf j1(v1.c[0] + v3.c[0] - v2.c[0] - v4.c[0],
              v1.c[1] + v3.c[1] - v2.c[1] - v4.c[1],
              v1.c[2] + v3.c[2] - v2.c[2] - v4.c[2]);
  vec3d_bf j2(v1.c[0] - v3.c[0] - v2.c[0] + v4.c[0],
              v1.c[1] - v3.c[1] - v2.c[1] + v4.c[1],
              v1.c[2] - v3.c[2] - v2.c[2] + v4.c[2]);
  vec3d_bf j3((-v1.c[0]) + v3.c[0] - v2.c[0] + v4.c[0],
              (-v1.c[1]) + v3.c[1] - v2.c[1] + v4.c[1],
              (-v1.c[2]) + v3.c[2] - v2.c[2] + v4.c[2]);

  j1 *= 0.5;
  j2 *= 0.5;
  j3 *= 0.5;

  const bigfloat num = (j1 * j1) + (j2 * j2) + (j3 * j3);
  bigfloat det = j1.tripleProd(j3, j2);
  if (det.get_d() <= 0) {
    return DBL_MAX;
  }

  return num.get_d() / pow(det.get_d(), (2.0 / 3.0));
}
