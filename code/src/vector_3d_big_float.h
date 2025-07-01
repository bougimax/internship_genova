#pragma once

#include "numeric_wrapper.h"
#include "vector_3d.h"
#include <cassert>

class vec3d_bf {
public:
  bigfloat c[3]; // 3 coordinates

  inline vec3d_bf() {}
  inline vec3d_bf(const bigfloat x, const bigfloat y, const bigfloat z) {
    c[0] = x;
    c[1] = y;
    c[2] = z;
  }
  inline vec3d_bf(const pointType *p) {
    if (p->isExplicit3D()) {
      const explicitPoint &e = p->toExplicit3D();
      c[0] = e.X();
      c[1] = e.Y();
      c[2] = e.Z();
    } else {
      bigfloat d;
      p->getBigfloatLambda(c[0], c[1], c[2], d);
    }
  }

  inline vec3d_bf operator+(const vec3d_bf &v) const {
    return vec3d_bf(c[0] + v.c[0], c[1] + v.c[1], c[2] + v.c[2]);
  }
  inline vec3d_bf operator-(const vec3d_bf &v) const {
    return vec3d_bf(c[0] - v.c[0], c[1] - v.c[1], c[2] - v.c[2]);
  }
  inline vec3d_bf operator*(const bigfloat d) const {
    return vec3d_bf(c[0] * d, c[1] * d, c[2] * d);
  }
  inline void operator*=(const bigfloat d) {
    c[0] = d * c[0];
    c[1] = d * c[1];
    c[2] = d * c[2];
  }

  inline bigfloat dot(const vec3d_bf &p) const {
    return (c[0] * p.c[0] + c[1] * p.c[1] + c[2] * p.c[2]);
  }
  inline vec3d_bf cross(const vec3d_bf &p) const {
    return vec3d_bf(c[1] * p.c[2] - c[2] * p.c[1],
                    c[2] * p.c[0] - c[0] * p.c[2],
                    c[0] * p.c[1] - c[1] * p.c[0]);
  }
  inline bigfloat tripleProd(const vec3d_bf &v2, const vec3d_bf &v3) const {
    return ((v2.c[0] * v3.c[1] * c[2]) - (v3.c[0] * v2.c[1] * c[2])) +
           ((v3.c[0] * c[1] * v2.c[2]) - (c[0] * v3.c[1] * v2.c[2])) +
           ((c[0] * v2.c[1] * v3.c[2]) - (v2.c[0] * c[1] * v3.c[2]));
  }

  inline bigfloat operator*(const vec3d_bf &d) const { return dot(d); }
  inline vec3d_bf operator&(const vec3d_bf &d) const { return cross(d); }

  // Squared length
  inline bigfloat sq_length() const { return dot(*this); }

  // Squared distance
  inline bigfloat dist_sq(const vec3d_bf &v) const {
    return ((*this) - v).sq_length();
  }
};
