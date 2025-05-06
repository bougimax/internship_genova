#ifndef _DELAUNAY_
#define _DELAUNAY_

#include "numeric_wrapper.h"
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
// of hxt_SeqDel (Copyright (C) 2018 Célestin Marot). hxt_SeqDel is a sequential
// Delaunay triangulator hosted at
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
                 //                 vertex v_{i,j} in ith tetrahedra tetrahedra
  mutable std::vector<uint32_t> mark_tetrahedra;    // Marks on tets
  mutable std::vector<unsigned char> marked_vertex; // Marks on vertices

  // Gift-wrapping fields
  std::vector<int> memo_o3d;
  std::vector<std::vector<int>>
      memo_o3d_v_origbndt; // i-th vector is {orient3d(original_cav_tri_1,v_i),
                           // ..., orient3d(original_cav_tri_n,v_i)}

  std::vector<uint64_t> Del_deleted;

  const bool has_outer_vertices; // This is TRUE if mesh vertices must survive
                                 // after destruction

  // Constructor and destructor
  TetMesh() : has_outer_vertices(false) {};
  TetMesh(bool h) : has_outer_vertices(h) {};
  ~TetMesh() {
    if (!has_outer_vertices)
      flushVertices();
  };

  /////// Global functions ///////

  // Number of vertices (infinite vertex is not counted)
  uint32_t numVertices() const { return (uint32_t)vertices.size(); }

  // Number of tetrahedra including ghosts
  uint32_t numTets() const { return (uint32_t)(tet_node.size() >> 2); }

  // Number of non-ghost tetrahedra
  uint32_t countNonGhostTets() const {
    return numTets() - (uint32_t)std::count(tet_node.begin(), tet_node.end(),
                                            INFINITE_VERTEX);
  }

  // Fill the vertex vector with newly-created genericPoints
  void init_vertices(const double *coords, uint32_t num_v);

  // Destroy vertices
  void flushVertices() {
    for (pointType *p : vertices)
      delete p;
  }

  // Init the mesh with a tet connecting four non coplanar points in vertices
  void init(uint32_t &unswap_k, uint32_t &unswap_l);

  // Create a Delaunay tetrahedrization by incremental insertion
  void tetrahedrize();

  // Save the mesh to a .tet file
  // If inner_only is set, only tets tagged as DT_IN are saved
  bool saveTET(const char *filename, bool inner_only = false) const;

  // Save the mesh to a .mesh file (MEDIT format)
  // If inner_only is set, only tets tagged as DT_IN are saved
  bool saveMEDIT(const char *filename, bool inner_only = false) const;

  // As above, but uses a binary format to avoid rounding
  bool saveBinaryTET(const char *filename, bool inner_only = false) const;

  // Save the interface between DT_IN and DT_OUT as an OFF file
  bool saveBoundaryToOFF(const char *filename) const;

  // As above, but saves rational coordinates and distinguishes between inner
  // and outer tets
  bool saveRationalTET(const char *filename, bool inner_only = false);

  // Marks internal tets ad DT_IN and external as DT_OUT and return the number
  // of internal tets. cornerMask must be TRUE for each corner whose opposite
  // face is a constraint.
  size_t markInnerTets(std::vector<bool> &cornerMask,
                       uint64_t single_start = UINT64_MAX);

  // Clear deleted tets after insertions
  void removeDelTets();
  void removeManyDelTets();

  // Clear deleted vertices after removal
  void removeDelVertices();

  // Resize the whole structure to contain 'new_size' tets
  void resizeTets(uint64_t new_size);
  void reserveTets(uint64_t new_capacity);

  // Return TRUE if at least one tet becomes flat or inverted after having
  // snapped its vertices to their closest floating-point representable
  // positions. Init num_flipped and num_flattened with the overall number of
  // flips or flattings.
  bool hasBadSnappedOrientations(size_t &num_flipped,
                                 size_t &num_flattened) const;

  // Check whether the structure is coherent (use for debugging purposes)
  void checkMesh(bool checkDelaunay = true) const;

  /////// Local (element-based) functions ///////

  // Return ith vertex id of tetrahedra t
  vertex get_i_th_vertex_of_tetrahedra(tetrahedra t_index, uint32_t i) {
    assert(0 <= i && i < 4);
    return tet_node[(t_index << 2) + i];
  }
  // Return ith corner id of tetrahedra t
  corner get_i_th_corner_of_tetrahedra(tetrahedra t_index, uint32_t i) {
    assert(0 <= i && i < 4);
    return (t_index << 2) + i;
  }

  uint32_t get_index_of_corner_in_tet(corner c) const { return c & 3; }
  // Return tetrahedra index of the corner c
  tetrahedra get_tetrahedra_index_from_corner(corner c) {
    return (c & (~3)) >> 2;
  }

  // TRUE if tet is ghost
  bool isGhost(tetrahedra t_index) const {
    return tet_node[(t_index << 2) + 3] == INFINITE_VERTEX;
  }

  bool has_infinite_vertex(tetrahedra t_index) const {
    return (tet_node[(t_index << 2)] == INFINITE_VERTEX) ||
           (tet_node[(t_index << 2) + 1] == INFINITE_VERTEX) ||
           (tet_node[(t_index << 2) + 2] == INFINITE_VERTEX) ||
           (tet_node[(t_index << 2) + 3] == INFINITE_VERTEX);
  }

  // TRUE if t has vertex v
  bool tetHasVertex(tetrahedra t_index, vertex v) const;

  // Init 'ov' with the two vertices of tet which are not in 'v'
  void oppositeTetEdge(tetrahedra t_index, const vertex v[2],
                       vertex ov[2]) const;

  void oppositeTetEdgePair(tetrahedra t_index, const TetMesh::edge &edge,
                           TetMesh::edge &opposite_edge) const;

  // Let t and n be face-adjacent tets.
  // This function returns the corner in t which is opposite to n
  corner getCornerFromOppositeTet(tetrahedra t_index, tetrahedra n_index) const;

  // Return the i'th tet in neighbors 'n'
  inline corner getIthNeighbor(const corner *n, const uint64_t i) const {
    // n[i] recovers the ith opposite corner, then _ & (~3) recovers the tet
    // base of which this opposite corner is from
    return n[i] & (~3);
  }

  // Fill v with the three other vertices different from tet_node[c] of the
  // tetrahedra that contains corner c
  void getFaceVertices(corner c, vertex v[3]) const;

  // Fill 'nt' with the two tets that share the vertices v1,v2,v3
  bool getTetsFromFaceVertices(vertex v1, vertex v2, vertex v3,
                               tetrahedra *nt) const;

  // Return the corner of t which is opposite to its face with vertices v1,v2,v3
  corner tetOppositeCorner(tetrahedra t_index, vertex v1, vertex v2,
                           vertex v3) const;

  // Return the corner corresponding to vertex 'v' in the tet whose base corner
  // is tb
  corner tetCornerAtVertex(corner tb, vertex v) const {
    return ((tet_node[tb] == v) * (tb)) + ((tet_node[tb + 1] == v) * (tb + 1)) +
           ((tet_node[tb + 2] == v) * (tb + 2)) +
           ((tet_node[tb + 3] == v) * (tb + 3));

    // while (tet_node[tb] != v) tb++;
    // return tb;
  }

  corner get_corner_in_tet(tetrahedra t, vertex v) const {
    return tetCornerAtVertex(get_base_corner(t), v);
  }
  uint32_t get_index_of_vertex_in_tet(vertex v, tetrahedra t_index) const {
    return get_corner_in_tet(t_index, v) & 3;
  }

  corner get_base_corner(tetrahedra t) const { return t << 2; }

  // Set the adjacency between the two corners c1 and c2
  void setMutualNeighbors(const corner c1, const corner c2) {
    tet_neigh[c1] = c2;
    tet_neigh[c2] = c1;
  }

  // Direct pointer to nodes and neighs
  vertex *getTetNodes(corner c) { return tet_node.data() + c; }
  corner *getTetNeighs(corner c) { return tet_neigh.data() + c; }
  const vertex *getTetNodes(corner c) const { return tet_node.data() + c; }
  const corner *getTetNeighs(corner c) const { return tet_neigh.data() + c; }

  // tetNi is a sum modulo 3 - used to traverse the nodes of a tet
  static size_t tetN1(const size_t i) { return (i + 1) & 3; }
  static size_t tetN2(const size_t i) { return (i + 2) & 3; }
  static size_t tetN3(const size_t i) { return (i + 3) & 3; }

  // tetONi - as above, but results in a coherent orientation
  static size_t tetON1(const size_t i) { return tetN1(i); }
  static size_t tetON2(const size_t i) { return (i & 2) ^ 3; }
  static size_t tetON3(const size_t i) { return (i + 3) & 2; }

  // Push a new isolated vertex in the structure
  void pushVertex(pointType *p) {
    vertices.push_back(p);
    inc_tet.push_back(UINT64_MAX);
    marked_vertex.push_back(0);
  }

  // Inserts an isolated vertex which is already in the vertices array.
  // ct is a hint for the algorithm to start searching the tet containing vi
  void insertExistingVertex(const uint32_t vi, uint64_t &ct);

  // Starting from 'c', move by adjacencies until a tet is found that
  // contains vertex v_id. Return that tet.
  uint64_t searchTetrahedron(corner c, const vertex v_id);

  // Incident tetrahedra at a vertex
  void VT(vertex v, std::vector<tetrahedra> &vt) const;

  // Same as VT, but this one includes ghost tets as well
  void VTfull(vertex v, std::vector<tetrahedra> &vt) const;

  // Adjacent vertices
  void VV(vertex v, std::vector<vertex> &vv) const;

  // Incident tetrahedra at an edge
  void ET(vertex v1, vertex v2, std::vector<tetrahedra> &et) const;
  void ETfull(vertex v1, vertex v2, std::vector<tetrahedra> &et) const;

  // Incident tetrahedra at an edge represented as ordered sequence of corners
  void ETcorners(vertex v1, vertex v2, std::vector<corner> &et) const;

  // TRUE if v1 and v2 are connected by an edge
  bool hasEdge(vertex v1, vertex v2) const;

  // Swap the position of t1 and t2 in the structure and update all relations
  // accordingly
  void swapTets(const tetrahedra t1, const tetrahedra t2);

  // Mark/unmark/check one single bit in tet mask
  inline void mark_Tet_1(const uint64_t t) const {
    mark_tetrahedra[t] |= ((uint32_t)2);
  }
  inline void unmark_Tet_1(const uint64_t t) const {
    mark_tetrahedra[t] &= (~((uint32_t)2));
  }
  inline uint32_t is_marked_Tet_1(const uint64_t t) const {
    return mark_tetrahedra[t] & ((uint32_t)2);
  }
  inline void mark_Tet_2(const uint64_t t) const {
    mark_tetrahedra[t] |= ((uint32_t)4);
  }
  inline void unmark_Tet_2(const uint64_t t) const {
    mark_tetrahedra[t] &= (~((uint32_t)4));
  }
  inline uint32_t is_marked_Tet_2(const uint64_t t) const {
    return mark_tetrahedra[t] & ((uint32_t)4);
  }
  inline void mark_Tet_31(const uint64_t t) const {
    mark_tetrahedra[t] |= ((uint32_t)2147483648);
  }
  inline void unmark_Tet_31(const uint64_t t) const {
    mark_tetrahedra[t] &= (~((uint32_t)2147483648));
  }
  inline uint32_t is_marked_Tet_31(const uint64_t t) const {
    return mark_tetrahedra[t] & ((uint32_t)2147483648);
  }

  // Thes two functions mark/check one particular bit stating that a tet must be
  // deleted. Differently from above, here a tet is identified by its first
  // corner.
  void markToDelete(corner c) {
    mark_tetrahedra[c >> 2] |= ((uint32_t)1073741824);
  }
  bool isToDelete(corner c) const {
    return mark_tetrahedra[c >> 2] & ((uint32_t)1073741824);
  }

  // Marks a tet (identified by its first corner) as 'removed' and add it to the
  // queue for eventual deletion.
  void pushAndMarkDeletedTets(uint64_t c) {
    Del_deleted.push_back(c);
    markToDelete(c);
  }

  // Predicates operating on vertex indexes
  int vOrient3D(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4) const {
    return -pointType::orient3D(*vertices[v1], *vertices[v2], *vertices[v3],
                                *vertices[v4]);
  }
  int vInSphere(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4,
                uint32_t v5) const {
    return -pointType::inSphere(*vertices[v1], *vertices[v2], *vertices[v3],
                                *vertices[v4], *vertices[v5]);
  }

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

  bool optimizeNearDegenerateTets(bool verbose = false);

  /// Operations to optimize the mesh

  void optimizeMesh();

  // Execute first pass (refining) of optimization process as described in
  // sec 3.2 of tetwild MAX
  void first_pass(std::vector<double> &desired_lengths, double epsilon);

  double get_energy_from_splitting(tetrahedra tetrahedra, edge edge_to_split,
                                   pointType *potential_split_point);

  void first_pass_bis();

  // Execute second pass (coarsening) of optimization process as described in
  // sec 3.2 of tetwild MAX
  void second_pass();

  // Return TRUE if the tetrahedra t is fully inside the ball centered on v and
  // of radius length
  // MAX
  bool is_tet_in_sphere(vertex v, double length, tetrahedra t);

  // Split an edge ev0-ev1 into four subtets by inserting an isolated vertex v
  void splitEdge(vertex ev0, vertex ev1, vertex v);

  void splitEdgeBis(edge edge_to_split, vertex split_vertex);

  // 2-3 swap
  bool swapFace(uint64_t r, bool prevent_inversion,
                double min_energy = DBL_MAX);

  // Edge removal
  bool removeEdge(uint32_t v1, uint32_t v2, double min_energy = DBL_MAX);

  // Collapse an edge onto its first endpoint
  bool collapseOnV1(uint32_t v1, uint32_t v2, bool prevent_inversion,
                    double min_energy = DBL_MAX);

  // Fill 'bet' with boundary faces incident at v1-v2
  void boundaryETcorners(uint32_t v1, uint32_t v2,
                         std::vector<uint64_t> &bet) const;

  bool isOnBoundary(uint32_t v1, uint32_t v2) const {
    std::vector<uint64_t> bet;
    boundaryETcorners(v1, v2, bet);
    return !bet.empty();
  }

  bool isOnBoundary(uint32_t v) const {
    std::vector<uint64_t> bvt;
    boundaryVTcorners(v, bvt);
    return !bvt.empty();
  }

  // Fill 'bvt' with boundary faces incident at v
  void boundaryVTcorners(uint32_t v, std::vector<uint64_t> &bvt) const;

  // VV relation restricted to incident boundary triangles
  void boundaryVV(uint32_t v, std::vector<uint32_t> &bvv) const;

  // TRUE if v2 incident boundary triangles have no normals different
  // than those of boundary triangles incident at edge v1-v2.
  bool isDoubleFlatV2(uint32_t v1, uint32_t v2) const;

  double maxEnergyAtEdge(uint32_t v1, uint32_t v2) const;
  double maxEnergyAtFace(uint64_t f) const;
  double maxEnergyAtVertex(uint32_t v) const;

  bool isCollapsableOnV1(uint32_t v1, uint32_t v2) const {
    return (isDoubleFlatV2(v1, v2) || !isOnBoundary(v2));
  }

  void getMeshEdges(std::vector<std::pair<uint32_t, uint32_t>> &edges) const;

  size_t iterativelySwapMesh(double th_energy);

  double getTetEnergy(uint64_t t) const;

  double getTotalEnergy();
  double getMaxEnergy();
  double getMeanEnergy();

  void get_all_tets_energy(std::vector<double> &tets_energy);

  // Put in tets all the tetrahedras that are completely in the ball centered on
  // v of radius length
  // MAX
  void tets_in_ball(uint32_t v, double length, std::vector<uint64_t> &tets);

  void log_tetrahedra(tetrahedra t);
};

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

#endif // _DELAUNAY_
