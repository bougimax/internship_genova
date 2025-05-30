#ifndef _DELAUNAY_
#define _DELAUNAY_

#include "numeric_wrapper.h"
#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "updatable_priority_queue.h"
#include <assert.h>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#pragma intrinsic(fabs)

#define INFINITE_VERTEX UINT32_MAX
#define DT_UNKNOWN 0
#define DT_OUT 1
#define DT_IN 2

#define MARKBIT(m, twoPowBit) m |= ((uint32_t)twoPowBit)
#define UNMARKBIT(m, twoPowBit) m &= (~((uint32_t)twoPowBit))
#define ISMARKEDBIT(m, twoPowBit) m &((uint32_t)twoPowBit)

// Uncommenting the following macro definition makes the code use modified parts
// of hxt_SeqDel (Copyright (C) 2018 Célestin Marot). hxt_SeqDel is a 
// equential Delaunay triangulator hosted at
// https://git.immc.ucl.ac.be/hextreme/hxt_seqdel as of 2020. hxt_SeqDel is GPL
// licensed, meaning that if you uncomment the following line you accept the
// terms of the GPL license for the whole CDT code. If you need to use this code
// under the less restrictive LGPL license, please comment the following line.
// This will make the code slightly slower (from 1% to 3% depending on the
// cases). #define USE_MAROTS_METHOD

// Tetrahedral mesh data structure

class TetMesh {
public:
  typedef uint32_t vertex;
  typedef uint64_t corner;
  typedef uint64_t tetrahedra;
  typedef std::pair<vertex, vertex> edge;
  // General purpose fields
  std::vector<pointType *> vertices; // Vertices
  std::vector<tetrahedra> inc_tet;   // One tet incident upon each vertex
  std::vector<vertex>
      tet_node; // Tet corners -> [v_{0,0}, v_{0,1}, v_{0,2}, v_{0,3}; ... ;
                //                 v_{i,0}, v_{i,1}, v_{i,2}, v_{i,3}; ... ;
                //                 v_{n,0}, v_{n,1}, v_{n,2}, v_{n,3}]
                //                 Where v_{i,j} is the jth vertex of ith
                //                 tetrahedra
  std::vector<corner>
      tet_neigh; // Tet opposite corners ->
                 //               [c_{0,0}, c_{0,1}, c_{0,2}, c_{0,3}; ...;
                 //                 c_{i,0}, c_{i,1}, c_{i,2}, c_{i,3}; ... ;
                 //                 c_{n,0}, c_{n,1}, c_{n,2}, c_{n,3}]
                 //                 Where c_{i,j} is the corner index (relative
                 //                 to tet_node) of the opposite corner of
                 //                 vertex v_{i,j} in ith tetrahedra
  mutable std::vector<uint32_t> mark_tetrahedra;    // Marks on tets
  mutable std::vector<unsigned char> marked_vertex; // Marks on vertices

  std::map<vertex, vertex>
      temp_remap_vertex; // Temporary map for remapping when removing a vertex
  std::map<tetrahedra, tetrahedra>
      temp_remap_tetrahedra; // Temporary map for remapping when removing a
                             // tetrahedra
  std::vector<double>
      tetrahedras_energy; // Contains the energy for each tetrahedra t necessary
                          // to compute optimization pass, it's initialized with
                          // get_all_tets_energy()
  bool optimize_only_DT_IN = true;

  // Gift-wrapping fields
  std::vector<int> memo_o3d;
  std::vector<std::vector<int>>
      memo_o3d_v_origbndt; // i-th vector is {orient3d(original_cav_tri_1,v_i),
                           // ..., orient3d(original_cav_tri_n,v_i)}

  std::vector<uint64_t> Del_deleted;

  const bool has_outer_nvector3d
      /// This represents a floating-point representable 3D vector
      /// along with a minimal set of necessary functions.
      /// It is conservatively used as a fast replacement for slower exact
      /// methods.
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
      return hasLargerSphere(vector3d(p), vector3d(q), vector3d(r),
                             vector3d(s));
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
    static bool isAtMostTwiceDistanceThan(const pointType *p,
                                          const pointType *q,
                                          const pointType *r) {
      const vector3d pv(p), qv(q), rv(r);
      return pv.dist_sq(qv) * 4 < pv.dist_sq(rv);
    }
  };

  inline std::ostream &operator<<(std::ostream &os, const vector3d &p) {
    return os << (p.c[0]) << " " << (p.c[1]) << " " << (p.c[2]);
  }

#endif // _DELAUNAY_
