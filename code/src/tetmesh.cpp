#include "delaunay.h"
#include "numeric_wrapper.h"
#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/volume_mesh.h"
#include <algorithm>
#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <float.h>
#include <iomanip>
#include <iterator>
#include <memory>
#include <ostream>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <utility>

using namespace std;

void TetMesh::init_vertices(const double *coords, uint32_t num_v) {
  vertices.reserve(num_v);
  for (uint32_t i = 0; i < num_v; i++)
    vertices.push_back(
        new explicitPoint(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));
  inc_tet.resize(num_v, UINT64_MAX);
  marked_vertex.resize(num_v, 0);
}

void TetMesh::init(uint32_t &unswap_k, uint32_t &unswap_l) {
  const uint32_t n = numVertices();

  // Find non-coplanar vertices (we assume that no coincident vertices exist)
  int ori = 0;
  uint32_t i = 0, j = 1, k = 2, l = 3;

  for (; ori == 0 && k < n - 1; k++)
    for (l = k + 1; ori == 0 && l < n; l++)
      ori = vOrient3D(i, j, k, l);

  l--;
  k--;

  if (ori == 0)
    ip_error("TetMesh::init() - Input vertices do not define a volume.\n");

  unswap_k = k;
  unswap_l = l;
  std::swap(vertices[k], vertices[2]);
  k = 2;
  std::swap(vertices[l], vertices[3]);
  l = 3;

  if (ori < 0)
    std::swap(i, j); // Tets must have positive volume

  const uint32_t base_tet[] = {l, k, j, i,
                               l, j, k, INFINITE_VERTEX,
                               l, k, i, INFINITE_VERTEX,
                               l, i, j, INFINITE_VERTEX,
                               k, j, i, INFINITE_VERTEX};
  const uint64_t base_neigh[] = {19, 15, 11, 7, 18, 10, 13, 3, 17, 14,
                                 5,  2,  16, 6, 9,  1,  12, 8, 4,  0};

  resizeTets(5);
  std::memcpy(getTetNodes(0), base_tet, 20 * sizeof(uint32_t));
  std::memcpy(getTetNeighs(0), base_neigh, 20 * sizeof(uint64_t));

  // set the vertex-(one_of_the)incident-tetrahedron relation
  inc_tet[i] = inc_tet[j] = inc_tet[k] = inc_tet[l] = 0;
}

bool TetMesh::saveTET(const char *filename, bool inner_only) const {
  ofstream f(filename);

  if (!f) {
    std::cerr << "\nTetMesh::saveTET: Can't open file for writing.\n";
    return false;
  }

  f << numVertices() << " vertices\n";

  uint32_t ngnt = 0;
  for (uint32_t i = 0; i < numTets(); i++)
    if (mark_tetrahedra[i] == DT_IN)
      ngnt++;

  if (inner_only) {
    f << ngnt << " tets\n";
    for (uint32_t i = 0; i < numVertices(); i++)
      f << *vertices[i] << "\n";
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " "
          << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
  } else {
    f << ngnt << " inner tets\n";
    f << countNonGhostTets() - ngnt << " outer tets\n";
    for (uint32_t i = 0; i < numVertices(); i++)
      f << *vertices[i] << "\n";
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " "
          << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    for (uint32_t i = 0; i < numTets(); i++)
      if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " "
          << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
  }

  f.close();

  return true;
}

bool TetMesh::saveMEDIT(const char *filename, bool inner_only) const {
  ofstream f(filename);

  if (!f) {
    std::cerr << "\nTetMesh::saveMEDIT: Can't open file for writing.\n";
    return false;
  }

  f << "MeshVersionFormatted 2\nDimension\n3\n";

  f << "Vertices\n" << numVertices() << "\n";

  uint32_t ngnt = 0;
  for (uint32_t i = 0; i < numTets(); i++)
    if (mark_tetrahedra[i] == DT_IN)
      ngnt++;

  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  if (inner_only) {
    for (uint32_t i = 0; i < numVertices(); i++)
      f << *vertices[i] << " 1\n";
    f << "Tetrahedra\n" << ngnt << "\n";
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " "
          << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1
          << " 1\n";
  } else {
    for (uint32_t i = 0; i < numVertices(); i++)
      f << *vertices[i] << " 1\n";
    f << "Tetrahedra\n" << countNonGhostTets() << "\n";
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " "
          << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1
          << " 1\n";
    for (uint32_t i = 0; i < numTets(); i++)
      if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
        f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " "
          << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1
          << " 2\n";
  }

  f.close();

  return true;
}

bool TetMesh::saveBinaryTET(const char *filename, bool inner_only) const {
  ofstream f(filename, ios::binary);

  if (!f) {
    std::cerr << "\nTetMesh::saveBinaryTET: Can't open file for writing.\n";
    return false;
  }

  uint32_t num_v = numVertices(), num_t = 0;

  for (uint32_t i = 0; i < numTets(); i++)
    if (mark_tetrahedra[i] == DT_IN)
      num_t++;

  f << num_v << " vertices\n";

  if (inner_only) {
    f << num_t << " tets\n";
  } else {
    f << num_t << " inner tets\n";
    f << countNonGhostTets() - num_t << " outer tets\n";
  }

  double c[3];
  for (uint32_t i = 0; i < numVertices(); i++) {
    vertices[i]->getApproxXYZCoordinates(c[0], c[1], c[2], true);
    f.write((const char *)(&c), sizeof(double) * 3);
  }

  const uint32_t *tnd = tet_node.data();

  if (inner_only) {
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f.write((const char *)(tnd + i * 4), sizeof(uint32_t) * 4);
  } else {
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f.write((const char *)(tnd + i * 4), sizeof(uint32_t) * 4);
    for (uint32_t i = 0; i < numTets(); i++)
      if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
        f.write((const char *)(tnd + i * 4), sizeof(uint32_t) * 4);
  }

  f.close();

  return true;
}

bool TetMesh::saveBoundaryToOFF(const char *filename) const {
  ofstream f(filename);

  if (!f) {
    std::cerr << "\nTetMesh::saveBoundaryToOFF: Can't open file for writing.\n";
    return false;
  }

  f << "OFF\n" << numVertices() << " ";

  size_t num_tris = 0;
  for (uint64_t i = 0; i < tet_node.size(); i++)
    if (i > tet_neigh[i] &&
        mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2])
      num_tris++;

  f << num_tris << " 0\n";

  for (uint32_t i = 0; i < numVertices(); i++)
    f << *vertices[i] << "\n";

  uint32_t fv[3];
  for (uint64_t i = 0; i < tet_node.size(); i++)
    if (i > tet_neigh[i] &&
        mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2]) {
      getFaceVertices(i, fv);
      f << "3 " << fv[0] << " " << fv[1] << " " << fv[2] << "\n";
    }
  f.close();

  return true;
}

bool TetMesh::saveRationalTET(const char *filename, bool inner_only) {
#ifdef USE_INDIRECT_PREDS
  ofstream f(filename);

  if (!f) {
    std::cerr << "\nTetMesh::saveRationalTET: Can't open file for writing.\n";
    return false;
  }

  f << numVertices() << " vertices\n";

  uint32_t ngnt = 0;
  for (uint32_t i = 0; i < numTets(); i++)
    if (mark_tetrahedra[i] == DT_IN)
      ngnt++;

  if (inner_only) {
    f << ngnt << " tets\n";
    for (uint32_t i = 0; i < numVertices(); i++) {
      bigrational c[3];
      vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
      f << c[0] << " " << c[1] << " " << c[2] << "\n";
    }
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " "
          << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
  } else {
    f << ngnt << " inner tets\n";
    f << countNonGhostTets() - ngnt << " outer tets\n";
    for (uint32_t i = 0; i < numVertices(); i++) {
      bigrational c[3];
      vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
      f << c[0] << " " << c[1] << " " << c[2] << "\n";
    }
    for (uint32_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_IN)
        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " "
          << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    for (uint32_t i = 0; i < numTets(); i++)
      if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " "
          << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
  }

  f.close();
#endif

  return true;
}

bool TetMesh::tetHasVertex(tetrahedra t_index, vertex v) const {
  t_index <<= 2;
  return tet_node[t_index] == v || tet_node[t_index + 1] == v ||
         tet_node[t_index + 2] == v || tet_node[t_index + 3] == v;
}

void TetMesh::oppositeTetEdgePair(tetrahedra tet, const TetMesh::edge &edge,
                                  TetMesh::edge &opposite_edge) const {
  int i = 0, j = 0;
  tet <<= 2;
  while (i < 4) {
    const vertex w = tet_node[tet + i];
    if (w != edge.first && w != edge.second) {
      if (j == 0) {
        opposite_edge.first = w;
      } else {
        opposite_edge.second = w;
      }
      j++;
    }
    i++;
  }
  if (j != 2) {
    std::cout << "Found " << j << " opposite vertex for edge " << edge.first
              << ", " << edge.second << std::endl;
  }
  assert(j == 2);
}
void TetMesh::oppositeTetEdge(const tetrahedra tet, const vertex v[2],
                              vertex ov[2]) const {
  int i = 0, j = 0;
  while (i < 4) {
    const vertex w = tet_node[tet + i];
    if (w != v[0] && w != v[1])
      ov[j++] = w;
    i++;
  }
  assert(j == 2);
}

TetMesh::corner TetMesh::getCornerFromOppositeTet(tetrahedra t_index,
                                                  tetrahedra n_index) const {
  t_index <<= 2;
  for (int i = 0; i < 4; i++)
    if ((tet_neigh[t_index + i] >> 2) == n_index)
      return tet_neigh[t_index + i];
  assert(0);
  return UINT64_MAX;
}

void TetMesh::getFaceVertices(corner c, vertex v[3]) const {
  uint64_t tv = c & 3;
  const vertex *Node = tet_node.data() + (c - tv);
  v[0] = Node[(++tv) & 3];
  v[1] = Node[(++tv) & 3];
  v[2] = Node[(++tv) & 3];
}

bool TetMesh::getTetsFromFaceVertices(vertex v1, vertex v2, vertex v3,
                                      tetrahedra *nt) const {
  static std::vector<tetrahedra>
      vt; // Static to avoid reallocation at each call
  VTfull(v1, vt);
  int i = 0;
  for (tetrahedra t : vt)
    if (tetHasVertex(t, v2) && tetHasVertex(t, v3))
      nt[i++] = t;
  vt.clear();
  return (i == 2);
}

TetMesh::corner TetMesh::tetOppositeCorner(tetrahedra t, vertex v1, vertex v2,
                                           vertex v3) const {
  const corner tb = t << 2;
  const vertex *n = tet_node.data() + tb;
  for (int i = 0; i < 3; i++)
    if (n[i] != v1 && n[i] != v2 && n[i] != v3)
      return tet_neigh[tb + i];
  assert(n[3] != v1 && n[3] != v2 && n[3] != v3);
  return tet_neigh[tb + 3];
}

void TetMesh::resizeTets(uint64_t new_size) {
  mark_tetrahedra.resize(new_size, 0);
  new_size <<= 2;
  tet_node.resize(new_size);
  tet_neigh.resize(new_size);
}

void TetMesh::reserveTets(uint64_t new_capacity) {
  mark_tetrahedra.reserve(new_capacity);
  new_capacity <<= 2;
  tet_node.reserve(new_capacity);
  tet_neigh.reserve(new_capacity);
}

uint64_t TetMesh::searchTetrahedron(corner tet, const vertex v_id) {
  if (tet_node[tet + 3] == INFINITE_VERTEX)
    tet = getIthNeighbor(getTetNeighs(tet), 3);

  uint64_t i, f0 = 4;
  do {
    const uint32_t *Node = getTetNodes(tet);
    if (Node[3] == INFINITE_VERTEX)
      return tet;

    const uint64_t *Neigh = getTetNeighs(tet);
    for (i = 0; i < 4; i++)
      if (i != f0 && vOrient3D(Node[tetON1(i)], Node[tetON2(i)],
                               Node[tetON3(i)], v_id) < 0) {
        tet = getIthNeighbor(Neigh, i);
        f0 = Neigh[i] & 3;
        break;
      }
  } while (i != 4);

  return tet;
}

void TetMesh::VT(vertex v, std::vector<tetrahedra> &vt) const {
  // std::cout << "Trying to get incident tetrahedras to " << v << ", there is "
  //           << numVertices() << " vertices" << std::endl;

  static std::vector<corner>
      vt_queue; // Static to avoid reallocation at each call
  tetrahedra act_tet = inc_tet[v];
  corner act_corner;

  vt_queue.push_back(tetCornerAtVertex(get_base_corner(act_tet), v));
  mark_Tet_31(act_tet);

  for (size_t i = 0; i < vt_queue.size(); i++) {

    act_corner = vt_queue[i];

    const uint64_t sb = get_index_of_corner_in_tet(act_corner);
    const corner *tet_opposite_corner = tet_neigh.data() + act_corner - sb;

    for (int j = sb + 1; j < sb + 4; j++) {

      const corner opposite_corner = tet_opposite_corner[j & 3];
      const tetrahedra opposite_tet =
          get_tetrahedra_index_from_corner(opposite_corner);

      if (tet_node[opposite_corner] != INFINITE_VERTEX &&
          !is_marked_Tet_31(opposite_tet)) {

        vt_queue.push_back(
            tetCornerAtVertex(get_base_corner_from_corner(opposite_corner), v));
        mark_Tet_31(opposite_tet);
      }
    }
  }

  tetrahedra neigh_tet;

  for (corner c : vt_queue) {
    neigh_tet = get_tetrahedra_index_from_corner(c);
    unmark_Tet_31(neigh_tet);
    vt.push_back(neigh_tet);
  }
  vt_queue.clear();
}

void TetMesh::OneRing(TetMesh::edge e, std::vector<vertex> &one_ring,
                      std::vector<tetrahedra> &incident_tetrahedras) {
  ETfull(e.first, e.second, incident_tetrahedras);
  TetMesh::edge opposite_edge;
  one_ring.clear();

  for (tetrahedra t : incident_tetrahedras) {
    if (!has_infinite_vertex(t)) {
      oppositeTetEdgePair(t, e, opposite_edge);
      if (marked_vertex[opposite_edge.first] == 0) {
        marked_vertex[opposite_edge.first] |= 1;
        one_ring.push_back(opposite_edge.first);
      }
      if (marked_vertex[opposite_edge.second] == 0) {
        marked_vertex[opposite_edge.second] |= 1;
        one_ring.push_back(opposite_edge.second);
      }
    }
  }
  for (tetrahedra t : incident_tetrahedras) {
    if (!has_infinite_vertex(t)) {
      oppositeTetEdgePair(t, e, opposite_edge);
      marked_vertex[opposite_edge.first] = 0;
      marked_vertex[opposite_edge.second] = 0;
    }
  }
}

void TetMesh::VV(vertex v, std::vector<vertex> &vv) {

  static std::vector<corner> corner_queue,
      seen_corner; // Static to avoid reallocation at each call

  tetrahedra start_tet = inc_tet[v];
  assert(start_tet != UINT64_MAX);

  const corner tet_basis_corner = get_base_corner(start_tet);
  const corner start_corner = tetCornerAtVertex(tet_basis_corner, v);
  corner_queue.push_back(start_corner);
  mark_Tet_31(start_tet);

  corner act_corner;
  tetrahedra act_tetrahedra, next_potential_tet;
  uint32_t index_corner_in_tet;
  vertex neighbour;

  while (!corner_queue.empty()) {
    act_corner = corner_queue.back();
    seen_corner.push_back(act_corner);
    act_tetrahedra = get_tetrahedra_index_from_corner(act_corner);
    corner_queue.pop_back();

    index_corner_in_tet = get_index_of_corner_in_tet(act_corner);

    for (uint32_t i = index_corner_in_tet + 1; i < index_corner_in_tet + 4;
         i++) {
      neighbour = get_i_th_vertex_of_tetrahedra(act_tetrahedra, i & 3);
      if (neighbour != INFINITE_VERTEX && !(marked_vertex[neighbour] & 128)) {
        marked_vertex[neighbour] |= 128;
        vv.push_back(neighbour);
      }
      next_potential_tet = get_tetrahedra_index_from_corner(
          tet_neigh[get_i_th_corner_of_tetrahedra(act_tetrahedra, i & 3)]);
      if (!is_marked_Tet_31(next_potential_tet)) {
        corner_queue.push_back(
            tetCornerAtVertex(get_base_corner(next_potential_tet), v));
        mark_Tet_31(next_potential_tet);
      }
    }
  }

  for (corner c : seen_corner)
    unmark_Tet_31(get_tetrahedra_index_from_corner(c));
  for (vertex neigh : vv)
    marked_vertex[neigh] &= 127;

  seen_corner.clear();
  corner_queue.clear();
}

void TetMesh::ET(vertex v1, vertex v2, std::vector<tetrahedra> &et) const {
  std::vector<tetrahedra> incident_tetrahedras_v1;
  VT(v1, incident_tetrahedras_v1);
  for (tetrahedra t : incident_tetrahedras_v1) {
    if (tetHasVertex(t, v2))
      et.push_back(t);
  }
}

void TetMesh::ETfull(vertex v1, vertex v2, std::vector<tetrahedra> &et) const {
  std::vector<tetrahedra> incident_tetrahedras_v1;
  VTfull(v1, incident_tetrahedras_v1);
  for (tetrahedra t : incident_tetrahedras_v1) {
    if (tetHasVertex(t, v2))
      et.push_back(t);
  }
}

void TetMesh::ETcorners(vertex v1, vertex v2,
                        std::vector<tetrahedra> &et) const {
  uint64_t t;
  VTfull(v1, et);
  for (uint64_t s : et)
    if (tetHasVertex(s, v2)) {
      t = (s << 2);
      break;
    }

  while (tet_node[t] == v1 || tet_node[t] == v2)
    t++;

  et.clear();

  uint64_t c0 = t;
  do {
    et.push_back(t);                   // Add tet
    uint64_t oc = tet_neigh[t] & (~3); // Get next base
    uint32_t cv = tet_node[t];
    t &= (~3);
    while (tet_node[t] == v1 || tet_node[t] == v2 || tet_node[t] == cv)
      t++;
    t = tetCornerAtVertex(oc,
                          tet_node[t]); // Get corresp corner at opposite tet
  } while (t != c0);
}

void TetMesh::VTfull(vertex v, std::vector<tetrahedra> &vt) const {
  static std::vector<uint64_t>
      vt_queue; // Static to avoid reallocation at each call
  uint64_t s, t = inc_tet[v];
  vt_queue.push_back(t);
  mark_Tet_31(t);

  size_t num_corner = tet_neigh.size();

  while (!vt_queue.empty()) {
    t = vt_queue.back();
    vt_queue.pop_back();
    vt.push_back(t);
    t <<= 2;
    for (int i = 0; i < 4; i++) {
      if (t + i < num_corner) {
        s = tet_neigh[t + i] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) {
          vt_queue.push_back(s);
          mark_Tet_31(s);
        }
      }
    }
  }

  for (uint64_t t : vt)
    unmark_Tet_31(t);
}

bool TetMesh::hasEdge(vertex v1, vertex v2) const {
  static std::vector<uint64_t>
      vt_queue; // Static to avoid reallocation at each call
  uint64_t t = inc_tet[v1];
  const uint64_t tb = t << 2;
  if (tet_node[tb] == v2 || tet_node[tb + 1] == v2 || tet_node[tb + 2] == v2 ||
      tet_node[tb + 3] == v2)
    return true;

  vt_queue.push_back(tetCornerAtVertex(tb, v1));
  mark_Tet_31(t);

  for (size_t i = 0; i < vt_queue.size(); i++) {
    t = vt_queue[i];
    const uint64_t sb = t & 3;
    const uint64_t *tg = tet_neigh.data() + t - sb;
    for (int j = 1; j < 4; j++) {
      const uint64_t tb = tg[(sb + j) & 3];
      const uint64_t tbb = tb >> 2;
      const uint32_t w = tet_node[tb];
      if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
        vt_queue.push_back(tetCornerAtVertex(tbb << 2, v1));
        mark_Tet_31(tbb);
        if (w == v2) {
          for (uint64_t t : vt_queue)
            unmark_Tet_31(t >> 2);
          vt_queue.clear();
          return true;
        }
      }
    }
  }

  for (uint64_t t : vt_queue)
    unmark_Tet_31(t >> 2);
  vt_queue.clear();
  return false;
}

void TetMesh::swapTets(const tetrahedra t1, const tetrahedra t2) {
  if (t1 == t2)
    return;

  const corner t1_id = t1 << 2;
  const corner t2_id = t2 << 2;

  // update VT base relation
  for (int i = 0; i < 3; i++)
    if (inc_tet[tet_node[t1_id + i]] == t1)
      inc_tet[tet_node[t1_id + i]] = t2;
  if (tet_node[t1_id + 3] != INFINITE_VERTEX &&
      inc_tet[tet_node[t1_id + 3]] == t1)
    inc_tet[tet_node[t1_id + 3]] = t2;

  for (int i = 0; i < 3; i++)
    if (inc_tet[tet_node[t2_id + i]] == t2)
      inc_tet[tet_node[t2_id + i]] = t1;
  if (tet_node[t2_id + 3] != INFINITE_VERTEX &&
      inc_tet[tet_node[t2_id + 3]] == t2)
    inc_tet[tet_node[t2_id + 3]] = t1;

  // Update nodes and marks
  for (int i = 0; i < 4; i++)
    std::swap(tet_node[t1_id + i], tet_node[t2_id + i]);
  std::swap(mark_tetrahedra[t1], mark_tetrahedra[t2]);

  // update neigh-neigh relations
  const uint64_t ng1[] = {tet_neigh[t1_id + 0], tet_neigh[t1_id + 1],
                          tet_neigh[t1_id + 2], tet_neigh[t1_id + 3]};
  const uint64_t ng2[] = {tet_neigh[t2_id + 0], tet_neigh[t2_id + 1],
                          tet_neigh[t2_id + 2], tet_neigh[t2_id + 3]};

  for (int i = 0; i < 4; i++)
    if ((ng2[i] >> 2) != t1)
      tet_neigh[ng2[i]] = t1_id + i;
  for (int i = 0; i < 4; i++)
    if ((ng1[i] >> 2) != t2)
      tet_neigh[ng1[i]] = t2_id + i;

  for (int i = 0; i < 4; i++)
    if ((ng2[i] >> 2) != t1)
      tet_neigh[t1_id + i] = tet_neigh[t2_id + i];
    else
      tet_neigh[t1_id + i] = (tet_neigh[t2_id + i] & 3) + (t2 << 2);

  for (int i = 0; i < 4; i++)
    if ((ng1[i] >> 2) != t2)
      tet_neigh[t2_id + i] = ng1[i];
    else
      tet_neigh[t2_id + i] = (ng1[i] & 3) + (t1 << 2);
}

size_t TetMesh::markInnerTets(std::vector<bool> &cornerMask,
                              uint64_t single_start) {
  std::vector<uint64_t> C;

  // All ghosts are DT_OUT
  for (size_t i = 0; i < numTets(); i++)
    mark_tetrahedra[i] = (isGhost(i)) ? DT_OUT : DT_UNKNOWN;

  if (single_start != UINT64_MAX)
    C.push_back(single_start);
  else
    for (size_t i = 0; i < numTets(); i++)
      if (mark_tetrahedra[i] == DT_OUT)
        C.push_back(i);

  for (size_t i = 0; i < C.size(); i++) {
    uint64_t t = C[i];
    for (int j = 0; j < 4; j++) {
      const uint64_t n = tet_neigh[t * 4 + j];
      const uint64_t n2 = n >> 2;
      if (mark_tetrahedra[n2] == DT_UNKNOWN) {
        if (!cornerMask[n]) {
          mark_tetrahedra[n2] = mark_tetrahedra[t];
        } else {
          mark_tetrahedra[n2] =
              ((mark_tetrahedra[t] == DT_IN) ? (DT_OUT) : (DT_IN));
        }
        C.push_back(n2);
      }
    }
  }

  return std::count(mark_tetrahedra.begin(), mark_tetrahedra.end(), DT_IN);
}

bool TetMesh::hasBadSnappedOrientations(size_t &num_flipped,
                                        size_t &num_flattened) const {
  const uint32_t *tn = tet_node.data();
  const uint32_t *end = tn + tet_node.size();
  num_flipped = num_flattened = 0;
  explicitPoint v[4];
  while (tn < end) {
    if (tn[3] != INFINITE_VERTEX) {
      for (int i = 0; i < 4; i++) {
        const pointType *p = vertices[tn[i]];
        if (p->isExplicit3D())
          v[i] = p->toExplicit3D();
        else
          p->apapExplicit(v[i]);
      }
      const int o = pointType::orient3D(v[0], v[1], v[2], v[3]);
      if (o > 0)
        num_flipped++;
      else if (o == 0)
        num_flattened++;
    }
    tn += 4;
  }

  return (num_flipped || num_flattened);
}

void TetMesh::checkMesh(bool checkDelaunay) {
  size_t i;
  const uint32_t num_vertices = (uint32_t)vertices.size();
  // Check tet nodes
  for (i = 0; i < numTets(); i++)
    if (!isToDelete(i << 2)) {
      const uint32_t *tn = tet_node.data() + i * 4;
      if (tn[0] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[1] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[2] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[3] != INFINITE_VERTEX && tet_node[i * 4 + 3] >= num_vertices)
        ip_error("Wrong tet node!\n");
      if (tn[0] == tn[1] || tn[0] == tn[2] || tn[0] == tn[3] ||
          tn[1] == tn[2] || tn[1] == tn[3] || tn[2] == tn[3]) {

        std::cout << "Tet " << i << " has same vertex for different corners"
                  << std::endl;
        log_tetrahedra(i);
        ip_error("Wrong tet node indexes!\n");
      }
    }

  // Check neighbors
  for (i = 0; i < numTets() * 4; i++)
    if (!isToDelete(i)) {
      if (tet_neigh[i] >= tet_neigh.size()) {
        std::cout << "The opposite corner is out of bound for corner " << i
                  << std::endl;
        ip_error("Wrong neighbor!\n");
      }
      if (tet_neigh[tet_neigh[i]] != i) {
        std::cout << "opposite of opposite of " << i << " is not " << i
                  << " it is " << tet_neigh[tet_neigh[i]]
                  << " which is the opposite of " << tet_neigh[i] << std::endl;
        ip_error("Wrong neighbor!\n");
      }
    }

  // Check neighbor-node coherence
  for (i = 0; i < numTets() * 4; i++)
    if (!isToDelete(i)) {
      if (tetHasVertex(tet_neigh[i] >> 2, tet_node[i])) {
        std::cout << "The tetrahedra opposite to " << (i >> 2) << " which is "
                  << (tet_neigh[i] >> 2) << " both contains the vertex "
                  << tet_node[i] << " it is the opposite of " << i
                  << " which is " << tet_neigh[i] << std::endl;
        std::cout << "Here is tetrahedra " << (i >> 2) << std::endl;
        log_tetrahedra(i >> 2);
        std::cout << "Here is tetrahedra " << (tet_neigh[i] >> 2) << std::endl;
        log_tetrahedra(tet_neigh[i] >> 2);
        ip_error("Incoherent neighbor!\n");
      } else {
        uint32_t face[3];
        uint32_t opposite_face[3];
        getFaceVertices(i, face);
        getFaceVertices(tet_neigh[i], opposite_face);
        if (!tetHasVertex(tet_neigh[i] >> 2, face[0])) {
          std::cout << "Face opposed to corner " << i << " (vertex "
                    << tet_node[i] << ")" << " which is composed of " << face[0]
                    << ", " << face[1] << ", " << face[2]
                    << " is not the same opposite face to corner "
                    << tet_neigh[i] << " (vertex " << tet_node[tet_neigh[i]]
                    << ") which is composed of " << opposite_face[0] << ", "
                    << opposite_face[1] << ", " << opposite_face[2]
                    << std::endl;
          ip_error("Incoherent face at neighbors!\n");
        }
      }
    }

  // Check vt*
  for (i = 0; i < num_vertices; i++)
    if (inc_tet[i] != UINT64_MAX) {
      if (inc_tet[i] >= numTets()) {
        std::cout << "Got tetrahedra " << inc_tet[i] << " for vertex " << i
                  << " but there is only " << numTets() << " tetrahedras"
                  << std::endl;
        ip_error("Wrong vt* (out of range)!\n");
      }
      if (isGhost(inc_tet[i]))
        ip_error("Wrong vt* (ghost tet)!\n");
      const uint32_t *tn = tet_node.data() + inc_tet[i] * 4;
      if (tn[0] != i && tn[1] != i && tn[2] != i && tn[3] != i) {

        std::cout << "Vertex " << i << " says it belongs to tet " << inc_tet[i]
                  << " but does not" << std::endl;
        ip_error("Wrong vt*!\n");
      }
    }

  // Check marks
  // for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2))
  //    if (mark_tetrahedra[i])
  //        ip_error("Marked tet\n");

  // Check geometry
  for (i = 0; i < numTets(); i++)
    if (!isToDelete(i << 2)) {
      const uint32_t *tn = tet_node.data() + i * 4;
      if (tn[3] != INFINITE_VERTEX &&
          vOrient3D(tn[0], tn[1], tn[2], tn[3]) <= 0) {
        std::cout << "Found a tet that is degenerate :" << std::endl;
        log_tetrahedra(i);
        std::cout << "Here its energy :" << getTetEnergy(i) << std::endl;
        ip_error("Inverted/degn tet\n");
      }
    }

  // Check energy
  for (i = 0; i < numTets(); i++) {
    if (std::abs(getTetEnergy(i) - tetrahedras_energy[i]) > 0.000001) {
      std::cout << "Energy differs for tet " << i << " to " << getTetEnergy(i)
                << " for real energy to " << tetrahedras_energy[i]
                << " for stored one" << std::endl;
      ip_error("Incoherent energy\n");
    }
  }

  if (checkDelaunay) {
    for (size_t i = 0; i < numTets(); i++)
      if (!isToDelete(i << 2)) {
        const uint32_t *n = tet_node.data() + (i * 4);
        if (n[3] == INFINITE_VERTEX)
          continue;
        for (int j = 0; j < 4; j++) {
          uint32_t ov = tet_node[tet_neigh[i * 4 + j]];
          if (ov != INFINITE_VERTEX && vertexInTetSphere(n, ov) > 0)
            ip_error("Non delaunay\n");
        }
      }
  }

  // printf("checkMesh passed\n");
}

void TetMesh::log_tetrahedra(tetrahedra t) {
  std::cout << "Tetrahedra " << t << " is composed of vertices:" << std::endl;
  for (int i = 0; i < 4; i++) {
    std::cout << "corner : " << get_i_th_corner_of_tetrahedra(t, i)
              << ", vertex : " << get_i_th_vertex_of_tetrahedra(t, i)
              << ", coords : "
              << vector3d(vertices[get_i_th_vertex_of_tetrahedra(t, i)])
              << std::endl;
  }
}

void TetMesh::getMeshEdges(
    std::vector<std::pair<uint32_t, uint32_t>> &all_edges) const {
  for (size_t t = 0; t < numTets(); t++) {
    const uint32_t *tn = tet_node.data() + (t << 2);
    if (tn[3] == INFINITE_VERTEX)
      continue;
    for (int i = 0; i < 4; i++)
      for (int j = i + 1; j < 4; j++)
        if (tn[i] < tn[j])
          all_edges.push_back(std::pair<uint32_t, uint32_t>(tn[i], tn[j]));
        else
          all_edges.push_back(std::pair<uint32_t, uint32_t>(tn[j], tn[i]));
  }
  std::sort(all_edges.begin(), all_edges.end());
  all_edges.erase(std::unique(all_edges.begin(), all_edges.end()),
                  all_edges.end());
}
void TetMesh::get_edges_from_tetrahedras(
    std::vector<std::pair<uint32_t, uint32_t>> &all_edges,
    std::vector<tetrahedra> &tets) const {
  for (tetrahedra t : tets) {
    const uint32_t *tn = tet_node.data() + (t << 2);
    if (tn[3] == INFINITE_VERTEX)
      continue;
    for (int i = 0; i < 4; i++)
      for (int j = i + 1; j < 4; j++)
        if (tn[i] < tn[j])
          all_edges.push_back(std::pair<uint32_t, uint32_t>(tn[i], tn[j]));
        else
          all_edges.push_back(std::pair<uint32_t, uint32_t>(tn[j], tn[i]));
  }
  std::sort(all_edges.begin(), all_edges.end());
  all_edges.erase(std::unique(all_edges.begin(), all_edges.end()),
                  all_edges.end());
}

void TetMesh::boundaryETcorners(uint32_t v1, uint32_t v2,
                                std::vector<uint64_t> &et) const {
  ETcorners(v1, v2, et);
  for (size_t i = 0; i < et.size();)
    if (mark_tetrahedra[et[i] >> 2] == mark_tetrahedra[tet_neigh[et[i]] >> 2]) {
      std::swap(et[i], et[et.size() - 1]);
      et.pop_back();
    } else
      i++;
}

// Fill 'bvt' with boundary faces incident at v
void TetMesh::boundaryVTcorners(uint32_t v, std::vector<uint64_t> &bvt) const {
  std::vector<uint64_t> vt;
  VTfull(v, vt);
  for (uint64_t t : vt)
    for (int i = 0; i < 4; i++) {
      const uint64_t c = (t << 2) + i;
      const uint64_t n = tet_neigh[c];
      if (tet_node[c] != v && c < n &&
          mark_tetrahedra[t] != mark_tetrahedra[n >> 2])
        bvt.push_back(c);
    }
}

// VV relation restricted to incident boundary triangles
void TetMesh::boundaryVV(uint32_t v, std::vector<uint32_t> &bvv) const {
  std::vector<uint64_t> vt;
  VTfull(v, vt);
  for (uint64_t t : vt)
    for (int i = 0; i < 4; i++) {
      const uint64_t c = (t << 2) + i;
      const uint64_t n = tet_neigh[c];
      if (tet_node[c] != v && c < n &&
          mark_tetrahedra[t] != mark_tetrahedra[n >> 2]) {
        for (int j = 0; j < 3; j++) {
          const uint32_t v1 = tet_node[(t << 2) + ((i + j) & 3)];
          if (v1 != INFINITE_VERTEX && v1 != v && !(marked_vertex[v1] & 128)) {
            marked_vertex[v1] |= 128;
            bvv.push_back(v1);
          }
        }
      }
    }
  for (uint32_t w : bvv)
    marked_vertex[w] &= 127;
}

bool TetMesh::isDoubleFlatV2(uint32_t v1, uint32_t v2) const {
  std::vector<uint64_t> et;
  boundaryETcorners(v1, v2, et);

  std::vector<uint32_t> ov(et.size());
  uint32_t v[3];
  for (size_t i = 0; i < et.size(); i++) {
    getFaceVertices(et[i], v);
    for (int k = 0; k < 3; k++)
      if (v[k] != v1 && v[k] != v2)
        ov[i] = v[k];
  }

  // Now 'ov' contains opposite vertices of all boundary triangles incident at
  // v1-v2
  std::vector<uint32_t> vv;
  boundaryVV(v2, vv);

  for (uint32_t w : ov)
    marked_vertex[w] |= 128;

  // All the vertices in VV(v2)
  bool foundall = true;
  for (uint32_t o : vv)
    if (o != v1 && !(marked_vertex[o] & 128)) {
      bool found = false;
      for (uint32_t p : ov)
        if (vOrient3D(v1, v2, o, p) == 0) {
          found = true;
          break;
        }
      if (!found) {
        foundall = false;
        break;
      }
    }
  for (uint32_t w : ov)
    marked_vertex[w] &= 127;

  return foundall;
}
