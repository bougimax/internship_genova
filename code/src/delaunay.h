#ifndef _DELAUNAY_
#define _DELAUNAY_

#include "numeric_wrapper.h"
#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "tetmesh.h"
#include "vector_3d.h"
#include "updatable_priority_queue.h"
#include <assert.h>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

using vertex = TetMesh::vertex;
using edge = TetMesh::edge;
using corner = TetMesh::corner;
using tetrahedra = TetMesh::tetrahedra;

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

class TetMeshTetrahedrizer {
public:
  TetMesh *mesh;

  // Gift-wrapping fields
  std::vector<int> memo_o3d;
  std::vector<std::vector<int>>
      memo_o3d_v_origbndt; // i-th vector is {orient3d(original_cav_tri_1,v_i),
                           // ..., orient3d(original_cav_tri_n,v_i)}

  std::vector<uint64_t> Del_deleted;

  // Constructor and destructor
  TetMeshTetrahedrizer(TetMesh &mesh_) { mesh = &mesh_; };
  ~TetMeshTetrahedrizer() {};
  // Init the mesh with a tet connecting four non coplanar points in vertices
  void init(uint32_t &unswap_k, uint32_t &unswap_l);
  // Create a Delaunay tetrahedrization by incremental insertion
  void tetrahedrize();

  // Predicates operating on vertex indexes
  int vOrient3D(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4) const {
    return -pointType::orient3D(*mesh->vertices[v1], *mesh->vertices[v2],
                                *mesh->vertices[v3], *mesh->vertices[v4]);
  }
  int vInSphere(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4,
                uint32_t v5) const {
    return -pointType::inSphere(*mesh->vertices[v1], *mesh->vertices[v2],
                                *mesh->vertices[v3], *mesh->vertices[v4],
                                *mesh->vertices[v5]);
  }
  // Marks internal tets ad DT_IN and external as DT_OUT and return the number
  // of internal tets. cornerMask must be TRUE for each corner whose opposite
  // face is a constraint.
  size_t markInnerTets(std::vector<bool> &cornerMask,
                       uint64_t single_start = UINT64_MAX);

  // marks a tet (identified by its first corner) as 'removed' and add it to the
  // queue for eventual deletion.
  void pushAndMarkDeletedTets(corner c) {
    Del_deleted.push_back(c);
    mesh->markToDelete(c);
  }
  // Clear deleted tets after insertions
  void removeDelTets();
  void removeManyDelTets();

  // Clear deleted vertices after removal
  void removeDelVertices();

  // Return TRUE if at least one tet becomes flat or inverted after having
  // snapped its vertices to their closest floating-point representable
  // positions. Init num_flipped and num_flattened with the overall number of
  // flips or flattings.
  bool hasBadSnappedOrientations(size_t &num_flipped,
                                 size_t &num_flattened) const;

  // Inserts an isolated vertex which is already in the vertices array.
  // ct is a hint for the algorithm to start searching the tet containing vi
  void insertExistingVertex(const uint32_t vi, uint64_t &ct);

  // Starting from 'c', move by adjacencies until a tet is found that
  // contains vertex v_id. Return that tet.
  uint64_t searchTetrahedron(corner c, const vertex v_id);

  // Use the order of the five cospherical points in 'indices' to
  // return a nonzero though coherent inSphere result.
  int symbolicPerturbation(uint32_t indices[5]) const;

  // This is as vInSphere(v[0], v[1], v[2], v[3], v_id) but is guaranteed to
  // return a nonzero value by relying on the symbolic perturbation above.
  int vertexInTetSphere(const uint32_t v[4], uint32_t v_id) const;

  // Same as above, but the four vertices are the vertices of 'tet'.
  int vertexInTetSphere(uint64_t tet, uint32_t v_id) const;

  // Collect all the vertices contained in the smallest sphere by ep0 and ep1
  // and return the one generating the largest circumcircle with ep0 and ep1.
  // Init tet with one tet having the encroaching point
  uint32_t findEncroachingPoint(const uint32_t ep0, const uint32_t ep1,
                                uint64_t &tet) const;

  // Start from c and turn around v1-v2 as long as adjacencies are well defined.
  // When an invalid adjacency is found, reinit it and exit.
  void seekAndSetMutualAdjacency(int p_o0, int p_o1, int p_o2,
                                 const uint32_t *v, uint64_t c, uint64_t o,
                                 const uint32_t *tet_node_data,
                                 uint64_t *tet_neigh_data);

  // Rebuild internal adjacencies for the cavity tet opposite to c
  void restoreLocalConnectivty(uint64_t c, const uint32_t *tet_node_data,
                               uint64_t *tet_neigh_data);

#ifdef USE_MAROTS_METHOD
  class DelTmp {
  public:
    uint32_t node[4];
    uint64_t bnd;

    DelTmp(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint64_t o)
        : node{a, b, c, d}, bnd(o) {}
  };

  std::vector<DelTmp> Del_tmp;
  uint64_t numDelTmp() const { return Del_tmp.size(); }
  void flushDelTmp() { Del_tmp.clear(); }
  uint64_t *delTmpVec() const { return (uint64_t *)Del_tmp.data(); }
  void bnd_push(uint32_t v_id, uint32_t node1, uint32_t node2, uint32_t node3,
                uint64_t bnd) {
    Del_tmp.push_back(DelTmp(v_id, node1, node2, node3, bnd));
  }

  void deleteInSphereTets(uint64_t tet, const uint32_t v_id);
  void tetrahedrizeHole(uint64_t *tet);
#endif

  // Set of functions implementing the face recovery by gift-wrapping
  void fill_memo_o3d_v_origbndt(const uint32_t v,
                                const std::vector<uint64_t> &original_bnd_tri);
  bool FAST_innerSegmentCrossesInnerTriangle(
      const uint32_t *s_ep, const uint64_t obndt_j,
      const std::vector<uint64_t> &original_bnd_tri);
  bool FAST_innerSegmentCrossesInnerTriangle(
      const pointType &cs0, const pointType &cs1, const pointType &cv0,
      const pointType &cv1, const pointType &cv2, int &o3d_tri_s0,
      int &o3d_tri_s1) const;
  bool
  aInnerTriASide_Crosses_InnerTriB(const pointType &vA0, const pointType &vA1,
                                   const pointType &vA2, const pointType &vB0,
                                   const pointType &vB1, const pointType &vB2);
  bool intersectionTEST_3(const pointType &u0, const pointType &u1,
                          const pointType &u2, const pointType &v0,
                          const pointType &v1, const pointType &v2,
                          const pointType &y, const int face_ori);
  bool isTetLocallyDelaunay(const uint32_t *tet_vrts,
                            const std::vector<uint32_t> &C_vrts,
                            const std::vector<uint64_t> &original_bnd_tri);
  bool isTetIntersecting(const uint32_t *tet_vrts,
                         const std::vector<uint64_t> &C_bnd_tri);
  void orient_bnd_tri(const uint64_t bnd_tri, uint32_t *v) const;
  bool is_the_connecting_vrt(const uint32_t *bnd_tri_v, const uint32_t w,
                             const std::vector<uint64_t> &C_bnd_tetfaces,
                             const std::vector<uint32_t> &C_vrts,
                             const std::vector<uint64_t> &original_C_bnd);
  void connect_bnd_tri(const uint64_t bnd_tri,
                       std::vector<uint64_t> &C_bnd_tetfaces,
                       std::vector<uint32_t> &C_vrts,
                       const std::vector<uint64_t> &original_C_bnd);

  void giftWrapping(const std::vector<uint32_t> &comm_vrts,
                    std::vector<uint32_t> &C1_vrts,
                    std::vector<uint32_t> &C2_vrts,
                    const std::vector<uint64_t> &C_bnd_tetface,
                    const uint64_t n_cav_tets, const uint64_t n_C1_bnd_tetface);
  bool isUpperCavityTet(const uint64_t t, std::vector<int> &v_orient) const;
  bool isLowerCavityTet(const uint64_t t, std::vector<int> &v_orient) const;
  void recoverFaceGiftWrap(std::vector<uint64_t> &i_tets,
                           std::vector<int> &v_orient);
};

#endif // _DELAUNAY_
