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

void TetMesh::tetrahedrize() {
  uint32_t uk, ul;
  init(uk, ul); // First tet is made of vertices 0, 1, uk, ul

  // Need to unswap immediately to keep correct indexing and
  // ensure symbolic perturbation is coherent
  if (ul != 3) {
    std::swap(vertices[ul], vertices[3]);
    std::swap(inc_tet[ul], inc_tet[3]);
    for (uint32_t &tn : tet_node)
      if (tn == 3)
        tn = ul;
      else if (tn == ul)
        tn = 3;
  }

  if (uk != 2) {
    std::swap(vertices[uk], vertices[2]);
    std::swap(inc_tet[uk], inc_tet[2]);
    for (uint32_t &tn : tet_node)
      if (tn == 2)
        tn = uk;
      else if (tn == uk)
        tn = 2;
  }

  uint64_t ct = 0;
  for (uint32_t i = 2; i < numVertices(); i++)
    if (i != uk && i != ul)
      insertExistingVertex(i, ct);

  removeDelTets();
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

void TetMesh::removeManyDelTets() {
  uint64_t last = tet_node.size() - 4;
  while (isToDelete(last))
    last -= 4;
  for (uint64_t t : Del_deleted)
    if (t < last && isToDelete(t)) {
      for (int i = 0; i < 4; i++) {
        tet_node[t + i] = tet_node[last + i];
        const uint64_t n = tet_neigh[last + i];
        tet_neigh[t + i] = n;
        tet_neigh[n] = t + i;
        if (tet_node[last + i] != INFINITE_VERTEX &&
            inc_tet[tet_node[last + i]] == last >> 2)
          inc_tet[tet_node[last + i]] = t >> 2;
      }
      mark_tetrahedra[t >> 2] = mark_tetrahedra[last >> 2];
      last -= 4;
      while (isToDelete(last))
        last -= 4;
    }

  resizeTets((last + 4) >> 2);
  Del_deleted.clear();
}

#ifndef USE_MAROTS_METHOD
void TetMesh::removeDelTets() { removeManyDelTets(); }
#else
void TetMesh::removeDelTets() {
  uint64_t j;
  uint64_t tn = numTets();
  for (uint64_t i = 0; i < Del_deleted.size(); i++) {
    uint64_t to_delete = Del_deleted[i];
    uint64_t lastTet = (--tn) * 4;

    if (isToDelete(lastTet)) {
      for (j = i; j < Del_deleted.size(); j++)
        if (Del_deleted[j] == lastTet)
          break;

      Del_deleted[j] = Del_deleted[i];
    } else {
      for (j = 0; j < 4; j++) {
        tet_node[to_delete + j] = tet_node[lastTet + j];

        uint64_t neigh = tet_neigh[lastTet + j];
        tet_neigh[to_delete + j] = neigh;
        tet_neigh[neigh] = to_delete + j;

        if (tet_node[lastTet + j] != INFINITE_VERTEX &&
            inc_tet[tet_node[lastTet + j]] == lastTet >> 2)
          inc_tet[tet_node[lastTet + j]] = to_delete >> 2;
      }
      mark_tetrahedra[to_delete >> 2] = mark_tetrahedra[lastTet >> 2];
    }
  }
  resizeTets(tn);
  Del_deleted.clear();
}
#endif

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

int TetMesh::symbolicPerturbation(uint32_t indices[5]) const {
  int swaps = 0;
  int n = 5;
  int count;
  do {
    count = 0;
    n--;
    for (int i = 0; i < n; i++) {
      if (indices[i] > indices[i + 1]) {
        std::swap(indices[i], indices[i + 1]);
        count++;
      }
    }
    swaps += count;
  } while (count);

  n = vOrient3D(indices[1], indices[2], indices[3], indices[4]);
  if (n)
    return (swaps % 2) ? (-n) : n;

  n = vOrient3D(indices[0], indices[2], indices[3], indices[4]);
  return (swaps % 2) ? (n) : (-n);
}

int TetMesh::vertexInTetSphere(const uint32_t Node[4], uint32_t v_id) const {
  int det = vInSphere(Node[0], Node[1], Node[2], Node[3], v_id);
  if (det)
    return det;
  uint32_t nn[5] = {Node[0], Node[1], Node[2], Node[3], v_id};
  det = symbolicPerturbation(nn);
  if (det == 0.0)
    ip_error("Symbolic perturbation failed! Should not happen.\n");
  return det;
}

int TetMesh::vertexInTetSphere(uint64_t tet, uint32_t v_id) const {
  const uint32_t *Node = getTetNodes(tet);
  int det;

  if (Node[3] == INFINITE_VERTEX) {
    if ((det = vOrient3D(Node[0], Node[1], Node[2], v_id)) != 0)
      return det;
    const uint32_t nn[4] = {Node[0], Node[1], Node[2],
                            tet_node[tet_neigh[tet + 3]]};
    return -vertexInTetSphere(nn, v_id);
  } else
    return vertexInTetSphere(Node, v_id);
}

#ifdef USE_MAROTS_METHOD
void TetMesh::deleteInSphereTets(uint64_t tet, const uint32_t v_id) {
  pushAndMarkDeletedTets(tet);

  for (uint64_t t = Del_deleted.size() - 1; t < Del_deleted.size(); t++) {
    uint64_t tet = Del_deleted[t];
    uint64_t *Neigh = getTetNeighs(tet);
    uint32_t *Node = getTetNodes(tet);

    uint64_t neigh = getIthNeighbor(Neigh, 0);
    if (!isToDelete(neigh)) {
      if (vertexInTetSphere(neigh, v_id) < 0)
        bnd_push(v_id, Node[1], Node[2], Node[3], Neigh[0]);
      else
        pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 1);
    if (!isToDelete(neigh)) {
      if (vertexInTetSphere(neigh, v_id) < 0)
        bnd_push(v_id, Node[2], Node[0], Node[3], Neigh[1]);
      else
        pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 2);
    if (!isToDelete(neigh)) {
      if (vertexInTetSphere(neigh, v_id) < 0)
        bnd_push(v_id, Node[0], Node[1], Node[3], Neigh[2]);
      else
        pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 3);
    if (!isToDelete(neigh)) {
      if (vertexInTetSphere(neigh, v_id) < 0) {
        if (Node[1] < Node[2])
          bnd_push(v_id, Node[0], Node[2], Node[1], Neigh[3]);
        else
          bnd_push(v_id, Node[1], Node[0], Node[2], Neigh[3]);
      } else
        pushAndMarkDeletedTets(neigh);
    }
  }
}

void TetMesh::tetrahedrizeHole(uint64_t *tet) {
  uint64_t clength = Del_deleted.size(); // Num tets removed
  uint64_t blength = numDelTmp();        // Num tets to insert

  uint64_t tn = numTets();

  if (blength > clength) {
    for (uint64_t i = clength; i < blength; i++, tn++)
      Del_deleted.push_back(tn << 2);

    clength = blength;
    resizeTets(tn);
  }

  uint64_t start = clength - blength;

  for (uint64_t i = 0; i < blength; i++) {
    const uint64_t tet = Del_deleted[i + start];
    uint32_t *Node = getTetNodes(tet);

    Node[0] = Del_tmp[i].node[0];
    Node[1] = Del_tmp[i].node[1];
    Node[2] = Del_tmp[i].node[2];
    Node[3] = Del_tmp[i].node[3];

    uint64_t bnd = Del_tmp[i].bnd;
    tet_neigh[tet] = bnd;
    tet_neigh[bnd] = tet;
    Del_tmp[i].bnd = tet;

    mark_tetrahedra[tet >> 2] = 0;

    if (tet_node[tet + 3] != INFINITE_VERTEX)
      for (uint32_t j = 0; j < 4; j++)
        inc_tet[tet_node[tet + j]] = tet >> 2;
  }

  uint64_t tlength = 0;
  const uint64_t middle = blength * 3 / 2;

  uint64_t *Tmp = delTmpVec();
  const unsigned index[4] = {2, 3, 1, 2};

  for (uint64_t i = 0; i < blength; i++) {
    uint64_t tet = Del_deleted[start + i];
    const uint32_t *Node = getTetNodes(tet);

    for (uint64_t j = 0; j < 3; j++) {
      uint64_t key = ((uint64_t)Node[index[j]] << 32) + Node[index[j + 1]];
      tet++;

      uint64_t k;
      for (k = 0; k < tlength; k++)
        if (Tmp[k] == key)
          break;

      if (k == tlength) {
        Tmp[tlength] = (key >> 32) + (key << 32);
        Tmp[middle + tlength] = tet;
        tlength++;
      } else {
        uint64_t pairValue = Tmp[middle + k];
        tet_neigh[tet] = pairValue;
        tet_neigh[pairValue] = tet;
        tlength--;
        if (k < tlength) {
          Tmp[k] = Tmp[tlength];
          Tmp[middle + k] = Tmp[middle + tlength];
        }
      }
    }
  }

  flushDelTmp();
  *tet = Del_deleted[start];
  Del_deleted.resize(start);
}

void TetMesh::insertExistingVertex(const uint32_t vi, uint64_t &ct) {
  ct = searchTetrahedron(ct, vi);
  deleteInSphereTets(ct, vi);
  tetrahedrizeHole(&ct);
  uint64_t lt = ct;
  if (tet_node[lt + 3] == INFINITE_VERTEX)
    lt = tet_neigh[lt + 3];
  inc_tet[vi] = lt >> 2;
}

#else
// Start from c and turn around v1-v2 as long as adjacencies are well defined.
// When an invalid adjacency is found, reinit it and exit.
void TetMesh::seekAndSetMutualAdjacency(int p_o0, int p_o1, int p_o2,
                                        const uint32_t *v, uint64_t c,
                                        uint64_t o,
                                        const uint32_t *tet_node_data,
                                        uint64_t *tet_neigh_data) {
  const uint32_t ov = v[p_o0], v1 = v[p_o1], v2 = v[p_o2];
  o += p_o0;

  c &= (~3);
  while (tet_node_data[c] != ov)
    c++;

  for (;;) {
    uint64_t t = c;
    if ((c = tet_neigh_data[c]) == UINT64_MAX) {
      tet_neigh_data[t] = o;
      tet_neigh_data[o] = t;
      return;
    }
    const uint32_t w = tet_node_data[c];
    c &= (~3);
    while (tet_node_data[c] == v1 || tet_node_data[c] == v2 ||
           tet_node_data[c] == w)
      c++;
  }
}

// Rebuild internal adjacencies for the cavity tet opposite to c
void TetMesh::restoreLocalConnectivty(uint64_t c, const uint32_t *tet_node_data,
                                      uint64_t *tet_neigh_data) {
  const uint64_t o = tet_neigh_data[c];
  const uint32_t *v = tet_node_data + o;
  const uint64_t *n = tet_neigh_data + o;
  if (n[1] == UINT64_MAX)
    seekAndSetMutualAdjacency(1, 2, 3, v, c, o, tet_node_data, tet_neigh_data);
  if (n[2] == UINT64_MAX)
    seekAndSetMutualAdjacency(2, 1, 3, v, c, o, tet_node_data, tet_neigh_data);
  if (n[3] == UINT64_MAX)
    seekAndSetMutualAdjacency(3, 1, 2, v, c, o, tet_node_data, tet_neigh_data);
}

// Collect all tets whose circumsphere contains v_id and replace them
// with a star of new tets originating at v_id
void TetMesh::insertExistingVertex(const uint32_t v_id, uint64_t &tet) {
  static std::vector<uint64_t>
      cavityCorners; // Static to avoid reallocation on each call
  static const int fi[4][3] = {{2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
  uint32_t *tet_node_data = tet_node.data();
  uint64_t *tet_neigh_data = tet_neigh.data();

  // Move by adjacencies to find the tet containing v_id
  if (tet_node_data[tet + 3] == INFINITE_VERTEX)
    tet = tet_neigh_data[tet + 3] & (~3);

  uint64_t i, f0 = 4;
  do {
    const uint32_t *Node = tet_node_data + tet;
    if (Node[3] == INFINITE_VERTEX)
      break;

    for (i = 0; i < 4; i++)
      if (i != f0 && vOrient3D(Node[tetON1(i)], Node[tetON2(i)],
                               Node[tetON3(i)], v_id) < 0) {
        const uint64_t ni = tet_neigh_data[tet + i];
        tet = ni & (~3);
        f0 = ni & 3;
        break;
      }
  } while (i != 4);

  tet >>= 2;

  // Expand by adjacencies to collect all tets whose circumsphere contains v_id
  size_t first = Del_deleted.size();
  pushAndMarkDeletedTets(tet << 2);

  for (size_t i = first; i < Del_deleted.size(); i++) {
    const uint64_t *nb = tet_neigh_data + Del_deleted[i];
    const uint64_t *nl = nb + 4;

    for (; nb < nl; nb++) {
      const uint64_t n0 = *nb >> 2;
      uint32_t &mtn0 = mark_tetrahedra[n0];
      if (mtn0 == 0) {
        if (vertexInTetSphere(n0 << 2, v_id) < 0) {
          mtn0 = 2;
          cavityCorners.push_back(*nb);
        } else {
          pushAndMarkDeletedTets(n0 << 2);
        }
      } else if (mtn0 == 2)
        cavityCorners.push_back(*nb);
    }
  }

  // Resize the mesh to host the new tets
  uint64_t ntb, newpos = tet_node.size();
  if (cavityCorners.size() > Del_deleted.size()) {
    resizeTets(numTets() + (cavityCorners.size() - Del_deleted.size()));
    tet_node_data = tet_node.data();
    tet_neigh_data = tet_neigh.data();
  }

  // Create the new tets
  for (const uint64_t c : cavityCorners) {
    mark_tetrahedra[c >> 2] = 0;
    if (Del_deleted.empty()) {
      ntb = newpos;
      newpos += 4;
    } else {
      ntb = Del_deleted.back();
      Del_deleted.pop_back();
    }
    const uint64_t cb = c & 3;
    const uint32_t *cr = tet_node_data + (c - cb);
    uint32_t *cn = tet_node_data + ntb;
    *cn++ = v_id;
    *cn++ = cr[fi[cb][0]];
    *cn++ = cr[fi[cb][1]];
    *cn++ = cr[fi[cb][2]];

    tet_neigh_data[ntb] = c;
    tet_neigh_data[c] = ntb;
    tet_neigh_data[ntb + 1] = tet_neigh_data[ntb + 2] =
        tet_neigh_data[ntb + 3] = UINT64_MAX;

    ntb >>= 2;
    if ((*(--cn)) != INFINITE_VERTEX) {
      inc_tet[*cn] = ntb;
      inc_tet[*(--cn)] = ntb;
      inc_tet[*(--cn)] = ntb;
      inc_tet[v_id] = ntb;
    }
    mark_tetrahedra[ntb] = 0;
  }

  // Restore the connectivity within the cavity
  for (uint64_t c : cavityCorners)
    restoreLocalConnectivty(c, tet_node_data, tet_neigh_data);

  tet = tet_neigh_data[cavityCorners.back()];

  cavityCorners.clear();
}
#endif
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

// void TetMesh::VV(vertex v, std::vector<vertex> &vv) const {
//   static std::vector<uint64_t>
//       vt_queue; // Static to avoid reallocation at each call
//   uint64_t t = inc_tet[v];
//
//   assert(t != UINT64_MAX);
//   const uint64_t tb = t << 2;
//
//   const uint64_t s = tetCornerAtVertex(tb, v);
//   vt_queue.push_back(s);
//   mark_Tet_31(t);
//
//   const uint32_t *tn = tet_node.data() + tb;
//   const uint64_t sb = s & 3;
//   for (int j = 1; j < 4; j++) {
//     const uint32_t w = tn[(sb + j) & 3];
//     marked_vertex[w] |= 128;
//     vv.push_back(w);
//   }
//
//   for (size_t i = 0; i < vt_queue.size(); i++) {
//     t = vt_queue[i];
//     const uint64_t sb = t & 3;
//     const uint64_t *tg = tet_neigh.data() + t - sb;
//     for (int j = 1; j < 4; j++) {
//       const uint64_t tb = tg[(sb + j) & 3];
//       const uint64_t tbb = tb >> 2;
//       const uint32_t w = tet_node[tb];
//       if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
//         vt_queue.push_back(tetCornerAtVertex(tb & (~3), v));
//         mark_Tet_31(tbb);
//         if (!(marked_vertex[w] & 128)) {
//           marked_vertex[w] |= 128;
//           vv.push_back(w);
//         }
//       }
//     }
//   }
//
//   for (uint64_t t : vt_queue)
//     unmark_Tet_31(t >> 2);
//   vt_queue.clear();
//   for (uint32_t w : vv)
//     marked_vertex[w] &= 127;
// }

void TetMesh::ET(vertex v1, vertex v2, std::vector<tetrahedra> &et) const {
  // VT(v1, et);
  // for (size_t i = 0; i < et.size();)
  //   if (!tetHasVertex(et[i], v2)) {
  //     std::swap(et[i], et[et.size() - 1]);
  //     et.pop_back();
  //   } else
  //     i++;
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
  // for (size_t i = 0; i < et.size();)
  //   if (!tetHasVertex(et[i], v2)) {
  //     std::swap(et[i], et[et.size() - 1]);
  //     et.pop_back();
  //   } else
  //     i++;
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

uint32_t TetMesh::findEncroachingPoint(const uint32_t ep0, const uint32_t ep1,
                                       uint64_t &tet_e) const {
  static std::vector<uint64_t>
      enc_queue; // Static to avoid reallocation upon each call

  // Start collecting tetrahedra incident at the endpoints
  VT(ep0, enc_queue);

  for (uint64_t j : enc_queue)
    mark_Tet_1(j);

  const vector3d p0 = vertices[ep0];
  const vector3d p1 = vertices[ep1];
  const double eslen = (p0 - p1).sq_length();

  vector3d ep;
  uint32_t enc_pt_i = UINT32_MAX;

  marked_vertex[ep0] = marked_vertex[ep1] = 1;

  // Collect all encroaching points while expanding around insphere vertices
  for (uint32_t ti = 0; ti < enc_queue.size(); ti++) {
    const uint64_t tet = enc_queue[ti];
    const uint64_t tb = tet << 2;

    // Check each tet vertex for 'isphereness' and keep track of the one with
    // largest sphere
    const uint32_t *tn = tet_node.data() + tb;
    for (uint32_t i = 0; i < 4; i++) {
      const uint32_t ui = tn[i];
      if (!marked_vertex[ui]) {
        const vector3d &pui = vertices[ui];
        if (((pui - p0).sq_length() + (pui - p1).sq_length()) <= eslen) {
          marked_vertex[ui] = 1;
          if (enc_pt_i == UINT32_MAX ||
              vector3d::hasLargerSphere(p0, p1, pui, ep)) {
            ep = pui;
            enc_pt_i = ui;
            tet_e = tb;
          }
        } else
          marked_vertex[ui] = 2;
      }
    }

    const int nvmask[] = {
        (marked_vertex[tn[0]] == 1), (marked_vertex[tn[1]] == 1),
        (marked_vertex[tn[2]] == 1), (marked_vertex[tn[3]] == 1)};
    const int totmarkeda = nvmask[0] + nvmask[1] + nvmask[2] + nvmask[3];

    // Expand on adjacent tets if at least one common vertex is insphere
    const uint64_t *tg = tet_neigh.data() + tb;
    for (uint32_t i = 0; i < 4; i++) {
      const uint64_t nc = tg[i];
      const uint64_t n = nc >> 2;
      if (is_marked_Tet_1(n) == 2 || tet_node[nc] == INFINITE_VERTEX)
        continue;
      const int totmarked = totmarkeda - nvmask[i];
      if (totmarked) {
        mark_Tet_1(n);
        enc_queue.push_back(n);
      }
    }
  }

  // Clear all marks
  marked_vertex[ep0] = marked_vertex[ep1] = 0;
  for (uint64_t j : enc_queue) {
    unmark_Tet_1(j);
    j <<= 2;
    marked_vertex[tet_node[j++]] = 0;
    marked_vertex[tet_node[j++]] = 0;
    marked_vertex[tet_node[j++]] = 0;
    marked_vertex[tet_node[j]] = 0;
  }
  enc_queue.clear();

  return enc_pt_i;
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// M E S H   O P T I M I Z A T I O N   F U N C T I O N S
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Energy to be minimized for mesh optimization - no need to be exact
//
// This returns DBL_MAX for flipped or degenerate tets
// For a regular tetrahedron returns 3.
// For generic non-degenerate tetrahedra returns a value in the range [3,
// DBL_MAX]

double tetEnergyb(const pointType *p1, const pointType *p2, const pointType *p3,
                  const pointType *p4) {}

bool TetMesh::optimizeNearDegenerateTets(bool verbose) {

  std::vector<explicitPoint> ev(vertices.size());
  for (size_t i = 0; i < numVertices(); i++)
    vertices[i]->apapExplicit(ev[i]);

  bool iterate;
  size_t nflip, nflat;
  uint32_t max_iter = 10;

  do {
    if (verbose)
      printf("VERTICES: %u\n", numVertices());
    iterate = false;
    // First, remove zero-length edges
    for (uint64_t t = 0; t < numTets(); t++)
      if (!isToDelete(t << 2) && !isGhost(t)) {
        const uint32_t *tn = tet_node.data() + (t << 2);
        int j, k;
        for (j = 0; j < 4; j++) {
          for (k = j + 1; k < 4; k++)
            if (ev[tn[j]] == ev[tn[k]]) {
              if (collapseOnV1(tn[j], tn[k], true)) {
                j = 4;
                break;
              } else if (collapseOnV1(tn[k], tn[j], true)) {
                j = 4;
                break;
              }
            }
          if (k < 4) {
            iterate = true;
            break;
          }
        }
      }
    if (iterate) {
      removeManyDelTets();
      removeDelVertices();
    }

    // Second, swap tets to remove slivers
    iterativelySwapMesh(double(1UL << (2 * (max_iter - 1))));

    const uint32_t *tn = tet_node.data();
    const uint32_t *end = tn + tet_node.size();
    nflip = nflat = 0;
    while (tn < end) {
      if (tn[3] != INFINITE_VERTEX) {
        const int o =
            pointType::orient3D(ev[tn[0]], ev[tn[1]], ev[tn[2]], ev[tn[3]]);
        if (o > 0)
          nflip++;
        else if (o == 0)
          nflat++;
      }
      tn += 4;
    }

    iterate = (nflip || nflat);

    if (verbose)
      printf("ATTEMPT N.: %u - NUM DGN: %zu\n", max_iter, nflip + nflat);
  } while (--max_iter && iterate);

  removeManyDelTets();
  removeDelVertices();
  if (iterate)
    return false;

  // Do the actual snap rounding
  for (uint32_t v = 0; v < numVertices(); v++)
    if (!vertices[v]->isExplicit3D()) {
      explicitPoint *np = new explicitPoint(ev[v]);
      delete vertices[v];
      vertices[v] = np;
    }

  return true;
}

bool TetMesh::removeEdge(uint32_t v1, uint32_t v2, double pre_energy) {
  // THIS can be optimize as follows:
  // 1) Extract ET corners
  // 2) Find an et corner that satisfies the requirements
  // 3) If found, continue

  uint32_t newv = numVertices();
  pushVertex(NULL); // This is just a dummy vertex. No need for real coordinates
  const size_t num_tets_before = numTets();
  // split_edge(std::make_pair(v1, v2), newv);

  bool succeeds = false;
  static std::vector<uint32_t> vv;
  VV(newv, vv);

  for (uint32_t w : vv)
    if (w != v1 && w != v2 && collapseOnV1(w, newv, true, pre_energy)) {
      succeeds = true;
      break;
    }

  if (!succeeds) {
    collapseOnV1(v1, newv, false);
  }

  vertices.pop_back();
  marked_vertex.pop_back();
  inc_tet.pop_back();

  vv.clear();

  return succeeds;
}

// Checks:
// 0) Neither v1 nor v2 can be INFINITE_VERTEX
// 1) if v1 and v2 are on boundary then the edge must also be on boundary
// 2) tets that share v2 must keep a positive volume
bool TetMesh::collapseOnV1(uint32_t v1, uint32_t v2, bool prevent_inversion,
                           double th_energy) {
  if (v1 == INFINITE_VERTEX || v2 == INFINITE_VERTEX)
    return false;

  std::vector<uint64_t> vtf1, vtf2, v1nv, v2nv;
  bool v1_on_boundary = false, v2_on_boundary = false, e_on_boundary = false;

  // Check if v1 is on boundary
  VTfull(v1, vtf1);
  for (uint64_t t : vtf1)
    if (isGhost(t)) {
      v1_on_boundary = true;
      break;
    }

  // Check if v2 is on boundary
  VTfull(v2, vtf2);
  for (uint64_t t : vtf2)
    if (isGhost(t)) {
      v2_on_boundary = true;
      break;
    }

  for (uint64_t t : vtf2)      // For all tets incident at v2
    if (tetHasVertex(t, v1)) { // If one of the edge is the one
                               // we are going to collapse
      if (isGhost(t))
        e_on_boundary = true;
      const uint64_t tb = t << 2; // tet basis
      const uint64_t oc1 = tet_neigh[tetCornerAtVertex(tb, v2)];
      const uint64_t oc2 = tet_neigh[tetCornerAtVertex(tb, v1)];
      v1nv.push_back(oc1);
      v2nv.push_back(oc2);
      if (tet_node[oc1] == tet_node[oc2])
        return false;
    }

  if (v1_on_boundary && v2_on_boundary && !e_on_boundary)
    return false;

  if (prevent_inversion) {
    for (uint64_t t : vtf2)
      if (!tetHasVertex(t, v1) && !isGhost(t)) {
        const uint32_t *nn = tet_node.data() + (t << 2);
        uint32_t nn4[4] = {nn[0], nn[1], nn[2], nn[3]};
        nn4[tetCornerAtVertex((t << 2), v2) & 3] = v1;
        if (tetEnergyb(vertices[nn4[0]], vertices[nn4[1]], vertices[nn4[2]],
                       vertices[nn4[3]]) >= th_energy)
          return false;
        if (vOrient3D(nn4[0], nn4[1], nn4[2], nn4[3]) <= 0)
          return false;
      }
  }

  for (size_t i = 0; i < v1nv.size(); i++)
    setMutualNeighbors(v1nv[i], v2nv[i]);
  for (uint64_t t : vtf2)
    if (tetHasVertex(t, v1))
      pushAndMarkDeletedTets(t << 2);
    else
      tet_node[tetCornerAtVertex(t << 2, v2)] = v1;

  inc_tet[v1] = inc_tet[v2] = UINT64_MAX;

  for (uint64_t t : vtf1)
    if (!isGhost(t) && !isToDelete(t << 2)) {
      const uint64_t tb = t << 2;
      inc_tet[tet_node[tb]] = inc_tet[tet_node[tb + 1]] =
          inc_tet[tet_node[tb + 2]] = inc_tet[tet_node[tb + 3]] = t;
    }
  for (uint64_t t : vtf2)
    if (!isGhost(t) && !isToDelete(t << 2)) {
      const uint64_t tb = t << 2;
      inc_tet[tet_node[tb]] = inc_tet[tet_node[tb + 1]] =
          inc_tet[tet_node[tb + 2]] = inc_tet[tet_node[tb + 3]] = t;
    }

  if (inc_tet[v1] == UINT64_MAX || inc_tet[v2] == UINT64_MAX) {
    std::cout << "Did not managed to assigned tet to vertex v1 or v2"
              << std::endl;
  }

  return true;
}

size_t TetMesh::iterativelySwapMesh(double th_energy) {
  // std::vector<std::pair<uint32_t, uint32_t>>
  //     all_edges; // All the non-infinite mesh edges
  // getMeshEdges(all_edges);
  //
  // std::vector<edgeWithLength> ets;
  // for (auto &e : all_edges)
  //   if (!isOnBoundary(e.first, e.second)) {
  //     const vector3d v[2] = {vector3d(vertices[e.first]),
  //                            vector3d(vertices[e.second])};
  //     ets.push_back(edgeWithLength(e.first, e.second, (v[0].dist_sq(v[1]))));
  //   }
  // std::sort(ets.begin(), ets.end());
  //
  // size_t swapped_edges = 0;
  // for (size_t i = ets.size(); i > 0; i--) {
  //   const edgeWithLength &e = ets[i - 1];
  //   const double pre_energy = maxEnergyAtEdge(e.v1, e.v2);
  //   if (pre_energy >= th_energy && removeEdge(e.v1, e.v2, pre_energy))
  //     swapped_edges++;
  // }
  //
  // removeManyDelTets();
  //
  // size_t swapped_faces = 0;
  // for (size_t t = 0; t < tet_node.size(); t++) {
  //   const uint64_t tb = t >> 2, nb = tet_neigh[t] >> 2;
  //   if (!isGhost(tb) && !isGhost(nb) &&
  //       mark_tetrahedra[tb] == mark_tetrahedra[nb]) {
  //     const double pre_energy = maxEnergyAtFace(t);
  //     if (pre_energy >= th_energy && swap_face(t, true, pre_energy))
  //       swapped_faces++;
  //   }
  // }
  //
  // return swapped_edges + swapped_faces;
}

void TetMesh::removeDelVertices() {
  std::vector<uint64_t> vt;
  uint32_t last = numVertices() - 1;
  while (inc_tet[last] == UINT64_MAX)
    last--;

  for (uint32_t i = 0; i < last; i++)
    if (inc_tet[i] == UINT64_MAX) {
      // last is the first tail vertex to be maintained
      VTfull(last, vt);
      for (uint64_t t : vt) {
        t <<= 2;
        for (int j = 0; j < 4; j++)
          if (tet_node[t + j] == last)
            tet_node[t + j] = i;
      }
      vt.clear();

      std::swap(vertices[i], vertices[last]);
      std::swap(inc_tet[i], inc_tet[last]);
      std::swap(marked_vertex[i], marked_vertex[last]);
      while (inc_tet[last] == UINT64_MAX)
        last--;
    }
  last++;

  vertices.resize(last);
  inc_tet.resize(last);
  marked_vertex.resize(last);
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

double TetMesh::getTetEnergy(uint64_t t) const {
  const uint32_t *n1 = tet_node.data() + (t << 2);
  if (n1[0] == INFINITE_VERTEX || n1[1] == INFINITE_VERTEX ||
      n1[2] == INFINITE_VERTEX || n1[3] == INFINITE_VERTEX)
    return DBL_MAX;
  return tetEnergyb(vertices[n1[0]], vertices[n1[1]], vertices[n1[2]],
                    vertices[n1[3]]);
}

double TetMesh::getTotalEnergy() {
  const uint32_t num_tets = numTets();
  double total_energy = 0., current_tet_energy;
  for (uint32_t t = 0; t < num_tets; t++) {
    current_tet_energy = getTetEnergy(t);
    if (current_tet_energy != DBL_MAX)
      total_energy += current_tet_energy;
  }
  return total_energy;
}
double TetMesh::getMaxEnergy() {
  const uint32_t num_tets = numTets();
  double max_energy = 0., current_tet_energy;
  for (uint32_t t = 0; t < num_tets; t++) {
    if (mark_tetrahedra[t] == DT_IN) {
      current_tet_energy = tetrahedras_energy[t];
      if (current_tet_energy != DBL_MAX)
        max_energy = std::max(max_energy, current_tet_energy);
    }
  }
  return max_energy;
}
double TetMesh::getMeanEnergy() {
  const uint32_t num_tets = numTets();
  double mean_energy = 0., current_tet_energy;
  uint32_t num_real_tets = 0;
  for (uint32_t t = 0; t < num_tets; t++) {
    if (mark_tetrahedra[t] == DT_IN) {
      current_tet_energy = getTetEnergy(t);
      if (current_tet_energy != DBL_MAX) {
        num_real_tets++;
        mean_energy += current_tet_energy;
      }
    }
  }
  return mean_energy / num_real_tets;
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
