#ifndef _TETMESH_
#define _TETMESH_
#include "numeric_wrapper.h"

#pragma intrinsic(fabs)

#define INFINITE_VERTEX UINT32_MAX
#define DT_UNKNOWN 0
#define DT_OUT 1
#define DT_IN 2

#define MARKBIT(m, twoPowBit) m |= ((uint32_t)twoPowBit)
#define UNMARKBIT(m, twoPowBit) m &= (~((uint32_t)twoPowBit))
#define ISMARKEDBIT(m, twoPowBit) m &((uint32_t)twoPowBit)

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

  bool has_outer_vertices;

  mutable std::vector<uint32_t> mark_tetrahedra;    // Marks on tets
  mutable std::vector<unsigned char> marked_vertex; // Marks on vertices

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
  // Thes two functions mark/check one particular bit stating that a tet must be
  // deleted. Differently from above, here a tet is identified by its first
  // corner.
  void markToDelete(corner c) {
    mark_tetrahedra[c >> 2] |= ((uint32_t)1073741824);
  }
  bool isToDelete(corner c) const {
    return mark_tetrahedra[c >> 2] & ((uint32_t)1073741824);
  }
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

  // Fill the vertex vector with newly-created genericPoints
  void init_vertices(const double *coords, uint32_t num_v);

  // Destroy vertices
  void flushVertices() {
    for (pointType *p : vertices)
      delete p;
  }
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

  // Resize the whole structure to contain 'new_size' tets
  void resizeTets(uint64_t new_size);
  void reserveTets(uint64_t new_capacity);

  // Check whether the structure is coherent (use for debugging purposes)
  void checkMesh(bool checkDelaunay = true);

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

  // Return the index of the corner in the tet he is in (it's literally c mod 4)
  uint64_t get_index_of_corner_in_tet(corner c) const { return c & 3; }

  // Return tetrahedra index of the corner c
  tetrahedra get_tetrahedra_index_from_corner(corner c) { return c >> 2; }
  tetrahedra get_tetrahedra_index_from_corner(corner c) const { return c >> 2; }

  // TRUE if tet is ghost
  bool isGhost(tetrahedra t_index) const {
    return tet_node[(t_index << 2) + 3] == INFINITE_VERTEX;
  }

  // TRUE if tet has an infinite vertex
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
    if (tet_node[tb] == v)
      return tb;
    if (tet_node[tb + 1] == v)
      return tb + 1;
    if (tet_node[tb + 2] == v)
      return tb + 2;
    if (tet_node[tb + 3] == v)
      return tb + 3;
    return UINT64_MAX;

    // while (tet_node[tb] != v) tb++;
    // return tb;
  }

  // Return the corner of v in the tetrahedra t
  corner get_corner_in_tet(tetrahedra t, vertex v) const {
    return tetCornerAtVertex(get_base_corner(t), v);
  }

  // Return the index of v in the tetrahedra t
  uint32_t get_index_of_vertex_in_tet(vertex v, tetrahedra t_index) const {
    return get_corner_in_tet(t_index, v) & 3;
  }

  // Get the corner basis of a tetrahedra
  corner get_base_corner(tetrahedra t) const { return t << 2; }

  // Get the corner basis from a corner
  corner get_base_corner_from_corner(corner c) const { return c & (~3); }

  // Set the adjacency between the two corners c1 and c2
  void setMutualNeighbors(const corner c1, const corner c2) {
    tet_neigh[c1] = c2;
    tet_neigh[c2] = c1;
  }

  // Direct pointer to nodes and neighs
  vertex *getTetNodes(corner c) { return tet_node.data() + c; }
  corner *getTetNeighs(corner c) { return tet_neigh.data() + c; }

  // Return a vector with points of t
  std::vector<pointType *> getTetPoints(tetrahedra t) {
    return {vertices[tet_node[t << 2]], vertices[tet_node[(t << 2) + 1]],
            vertices[tet_node[(t << 2) + 2]], vertices[tet_node[(t << 2) + 3]]};
  }
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

  // Incident tetrahedra at a vertex
  void VT(vertex v, std::vector<tetrahedra> &vt) const;

  // Same as VT, but this one includes ghost tets as well
  void VTfull(vertex v, std::vector<tetrahedra> &vt,
              bool verbose = false) const;

  // Adjacent vertices
  void VV(vertex v, std::vector<vertex> &vv);

  // Incident tetrahedra at an edge
  void ET(vertex v1, vertex v2, std::vector<tetrahedra> &et) const;
  void ETfull(vertex v1, vertex v2, std::vector<tetrahedra> &et,
              bool verbose = false) const;

  // Incident tetrahedra at an edge represented as ordered sequence of corners
  void ETcorners(vertex v1, vertex v2, std::vector<corner> &et,
                 bool verbose = false) const;

  void OneRing(TetMesh::edge e, std::vector<vertex> &one_ring,
               std::vector<tetrahedra> &incident_tetrahedras,
               bool verbose = false);

  // TRUE if v1 and v2 are connected by an edge
  bool hasEdge(vertex v1, vertex v2) const;

  // Swap the position of t1 and t2 in the structure and update all relations
  // accordingly
  void swapTets(const tetrahedra t1, const tetrahedra t2);

  // Fill 'bet' with boundary faces incident at v1-v2
  void boundaryETcorners(uint32_t v1, uint32_t v2, std::vector<uint64_t> &bet,
                         bool verbose = false) const;

  bool boundaryEdgePriority(TetMesh::edge e, TetMesh::edge f) {
    if (isOnBoundary(e.first, e.second)) {
      if (isOnBoundary(f.first, f.second)) {
        return (e.first < f.first ||
                (e.first == f.first && e.second < f.second));
      } else {
        return true;
      }
    } else {
      if (isOnBoundary(f.first, f.second)) {
        return false;
      } else {
        return (e.first < f.first ||
                (e.first == f.first && e.second < f.second));
      }
    }
  }

  bool isOnBoundary(uint32_t v1, uint32_t v2, bool verbose = false) const {
    std::vector<uint64_t> bet;
    boundaryETcorners(v1, v2, bet, verbose);
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

  // // TRUE if v2 incident boundary triangles have no normals different
  // // than those of boundary triangles incident at edge v1-v2.
  // bool isDoubleFlatV2(uint32_t v1, uint32_t v2) const;

  void getMeshEdges(std::vector<std::pair<uint32_t, uint32_t>> &edges) const;
  void get_edges_from_tetrahedras(
      std::vector<std::pair<uint32_t, uint32_t>> &all_edges,
      std::vector<tetrahedra> &tets) const;

  void log_tetrahedra(tetrahedra t);
};
#endif // _TETMESH_
