#pragma once

/// <summary>
/// vector3d
/// This represents a floating-point representable 3D vector
/// along with a minimal set of necessary functions.
/// It is conservatively used as a fast replacement for slower exact methods.
/// </summary>

class vector3d {
public:
  double c[3]; // 3 coordinates

  inline vector3d() {}
  inline vector3d(const double x, const double y, const double z) {
    c[0] = x;
    c[1] = y;
    c[2] = z;
  }
  inline vector3d(const pointType *p) {
    p->getApproxXYZCoordinates(c[0], c[1], c[2]);
  }

  explicitPoint *toExplicitPoint() {
    return new explicitPoint(c[0], c[1], c[2]);
  }

  inline vector3d operator+(const vector3d &v) const {
    return vector3d(c[0] + v.c[0], c[1] + v.c[1], c[2] + v.c[2]);
  }
  inline vector3d operator-(const vector3d &v) const {
    return vector3d(c[0] - v.c[0], c[1] - v.c[1], c[2] - v.c[2]);
  }
  inline vector3d operator*(const double d) const {
    return vector3d(c[0] * d, c[1] * d, c[2] * d);
  }

  inline void operator+=(const vector3d &v) {
    c[0] += v.c[0];
    c[1] += v.c[1];
    c[2] += v.c[2];
  }
  inline void operator*=(const double d) {
    c[0] *= d;
    c[1] *= d;
    c[2] *= d;
  }

  inline double dot(const vector3d &p) const {
    return (c[0] * p.c[0] + c[1] * p.c[1] + c[2] * p.c[2]);
  }
  inline vector3d cross(const vector3d &p) const {
    return vector3d(c[1] * p.c[2] - c[2] * p.c[1],
                    c[2] * p.c[0] - c[0] * p.c[2],
                    c[0] * p.c[1] - c[1] * p.c[0]);
  }
  inline double tripleProd(const vector3d &v2, const vector3d &v3) const {
    return ((v2.c[0] * v3.c[1] * c[2]) - (v3.c[0] * v2.c[1] * c[2])) +
           ((v3.c[0] * c[1] * v2.c[2]) - (c[0] * v3.c[1] * v2.c[2])) +
           ((c[0] * v2.c[1] * v3.c[2]) - (v2.c[0] * c[1] * v3.c[2]));
  }

  inline double operator*(const vector3d &d) const { return dot(d); }
  inline vector3d operator&(const vector3d &d) const { return cross(d); }

  // Squared length
  inline double sq_length() const { return dot(*this); }

  // Squared distance
  inline double dist_sq(const vector3d &v) const {
    return ((*this) - v).sq_length();
  }

  // TRUE if r is in (or on border of) sphere having p-q as diameter
  static inline bool inSmallestSphere(const pointType *p, const pointType *q,
                                      const pointType *r) {
    return inSmallestSphere(vector3d(p), vector3d(q), vector3d(r));
  }

  static inline bool inSmallestSphere(const vector3d &pv, const vector3d &qv,
                                      const vector3d &rv) {
    return ((rv - pv).sq_length() + (rv - qv).sq_length()) <=
           (pv - qv).sq_length();
  }

  // TRUE if smallest sphere by p,q,r is larger than smallest sphere by p,q,s
  static inline bool hasLargerSphere(const pointType *p, const pointType *q,
                                     const pointType *r, const pointType *s) {
    return hasLargerSphere(vector3d(p), vector3d(q), vector3d(r), vector3d(s));
  }

  static inline bool hasLargerSphere(const vector3d &pv, const vector3d &qv,
                                     const vector3d &rv, const vector3d &sv) {
    const vector3d pms = pv - sv, qms = qv - sv, pmr = pv - rv, qmr = qv - rv;
    const double lens = pms.sq_length() * qms.sq_length();
    if (lens == 0)
      return true;
    const double lenr = pmr.sq_length() * qmr.sq_length();
    if (lenr == 0)
      return false;
    const double dots = pms.dot(qms);
    const double dotr = pmr.dot(qmr);

    return (dots * dots) * lenr < (dotr * dotr) * lens;
  }

  // TRUE if p is closer to q than to r
  static bool isCloserThan(const pointType *p, const pointType *q,
                           const pointType *r) {
    const vector3d pv(p), qv(q), rv(r);
    return pv.dist_sq(qv) < pv.dist_sq(rv);
  }

  // TRUE if distance p-q is at most twice the distance p-r
  static bool isAtMostTwiceDistanceThan(const pointType *p, const pointType *q,
                                        const pointType *r) {
    const vector3d pv(p), qv(q), rv(r);
    return pv.dist_sq(qv) * 4 < pv.dist_sq(rv);
  }
};

inline std::ostream &operator<<(std::ostream &os, const vector3d &p) {
  return os << (p.c[0]) << " " << (p.c[1]) << " " << (p.c[2]);
}
