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

void TetMeshTetrahedrizer::init(uint32_t &unswap_k, uint32_t &unswap_l) {
  const uint32_t n = mesh->numVertices();

  // Find non-coplanar vertices (we assume that no coincident vertices exist)
  int ori = 0;
  uint32_t i = 0, j = 1, k = 2, l = 3;

  for (; ori == 0 && k < n - 1; k++)
    for (l = k + 1; ori == 0 && l < n; l++)
      ori = vOrient3D(i, j, k, l);

  l--;
  k--;

  if (ori == 0)
    ip_error("TetMeshTetrahedrizer::init() - Input vertices do not define a "
             "volume.\n");

  unswap_k = k;
  unswap_l = l;
  std::swap(mesh->vertices[k], mesh->vertices[2]);
  k = 2;
  std::swap(mesh->vertices[l], mesh->vertices[3]);
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

  mesh->resizeTets(5);
  std::memcpy(mesh->getTetNodes(0), base_tet, 20 * sizeof(uint32_t));
  std::memcpy(mesh->getTetNeighs(0), base_neigh, 20 * sizeof(uint64_t));

  // set the vertex-(one_of_the)incident-tetrahedron relation
  mesh->inc_tet[i] = mesh->inc_tet[j] = mesh->inc_tet[k] = mesh->inc_tet[l] = 0;
}

void TetMeshTetrahedrizer::tetrahedrize() {
  uint32_t uk, ul;
  init(uk, ul); // First tet is made of vertices 0, 1, uk, ul

  // Need to unswap immediately to keep correct indexing and
  // ensure symbolic perturbation is coherent
  if (ul != 3) {
    std::swap(mesh->vertices[ul], mesh->vertices[3]);
    std::swap(mesh->inc_tet[ul], mesh->inc_tet[3]);
    for (uint32_t &tn : mesh->tet_node)
      if (tn == 3)
        tn = ul;
      else if (tn == ul)
        tn = 3;
  }

  if (uk != 2) {
    std::swap(mesh->vertices[uk], mesh->vertices[2]);
    std::swap(mesh->inc_tet[uk], mesh->inc_tet[2]);
    for (uint32_t &tn : mesh->tet_node)
      if (tn == 2)
        tn = uk;
      else if (tn == uk)
        tn = 2;
  }

  uint64_t ct = 0;
  for (uint32_t i = 2; i < mesh->numVertices(); i++)
    if (i != uk && i != ul)
      insertExistingVertex(i, ct);

  removeDelTets();
}

size_t TetMeshTetrahedrizer::markInnerTets(std::vector<bool> &cornerMask,
                                           uint64_t single_start) {
  std::vector<uint64_t> C;

  // All ghosts are DT_OUT
  for (size_t i = 0; i < mesh->numTets(); i++)
    mesh->mark_tetrahedra[i] = (mesh->isGhost(i)) ? DT_OUT : DT_UNKNOWN;

  if (single_start != UINT64_MAX)
    C.push_back(single_start);
  else
    for (size_t i = 0; i < mesh->numTets(); i++)
      if (mesh->mark_tetrahedra[i] == DT_OUT)
        C.push_back(i);

  for (size_t i = 0; i < C.size(); i++) {
    uint64_t t = C[i];
    for (int j = 0; j < 4; j++) {
      const uint64_t n = mesh->tet_neigh[t * 4 + j];
      const uint64_t n2 = n >> 2;
      if (mesh->mark_tetrahedra[n2] == DT_UNKNOWN) {
        if (!cornerMask[n]) {
          mesh->mark_tetrahedra[n2] = mesh->mark_tetrahedra[t];
        } else {
          mesh->mark_tetrahedra[n2] =
              ((mesh->mark_tetrahedra[t] == DT_IN) ? (DT_OUT) : (DT_IN));
        }
        C.push_back(n2);
      }
    }
  }

  return std::count(mesh->mark_tetrahedra.begin(), mesh->mark_tetrahedra.end(),
                    DT_IN);
}

bool TetMeshTetrahedrizer::hasBadSnappedOrientations(
    size_t &num_flipped, size_t &num_flattened) const {
  const uint32_t *tn = mesh->tet_node.data();
  const uint32_t *end = tn + mesh->tet_node.size();
  num_flipped = num_flattened = 0;
  explicitPoint v[4];
  while (tn < end) {
    if (tn[3] != INFINITE_VERTEX) {
      for (int i = 0; i < 4; i++) {
        const pointType *p = mesh->vertices[tn[i]];
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
uint64_t TetMeshTetrahedrizer::searchTetrahedron(corner tet,
                                                 const vertex v_id) {
  if (mesh->tet_node[tet + 3] == INFINITE_VERTEX)
    tet = mesh->getIthNeighbor(mesh->getTetNeighs(tet), 3);

  uint64_t i, f0 = 4;
  do {
    const uint32_t *Node = mesh->getTetNodes(tet);
    if (Node[3] == INFINITE_VERTEX)
      return tet;

    const uint64_t *Neigh = mesh->getTetNeighs(tet);
    for (i = 0; i < 4; i++)
      if (i != f0 && vOrient3D(Node[mesh->tetON1(i)], Node[mesh->tetON2(i)],
                               Node[mesh->tetON3(i)], v_id) < 0) {
        tet = mesh->getIthNeighbor(Neigh, i);
        f0 = Neigh[i] & 3;
        break;
      }
  } while (i != 4);

  return tet;
}
void TetMeshTetrahedrizer::removeManyDelTets() {
  uint64_t last = mesh->tet_node.size() - 4;
  while (mesh->isToDelete(last))
    last -= 4;
  for (uint64_t t : Del_deleted)
    if (t < last && mesh->isToDelete(t)) {
      for (int i = 0; i < 4; i++) {
        mesh->tet_node[t + i] = mesh->tet_node[last + i];
        const uint64_t n = mesh->tet_neigh[last + i];
        mesh->tet_neigh[t + i] = n;
        mesh->tet_neigh[n] = t + i;
        if (mesh->tet_node[last + i] != INFINITE_VERTEX &&
            mesh->inc_tet[mesh->tet_node[last + i]] == last >> 2)
          mesh->inc_tet[mesh->tet_node[last + i]] = t >> 2;
      }
      mesh->mark_tetrahedra[t >> 2] = mesh->mark_tetrahedra[last >> 2];
      last -= 4;
      while (mesh->isToDelete(last))
        last -= 4;
    }

  mesh->resizeTets((last + 4) >> 2);
  Del_deleted.clear();
}

#ifndef USE_MAROTS_METHOD
void TetMeshTetrahedrizer::removeDelTets() { removeManyDelTets(); }
#else
void TetMeshTetrahedrizer::removeDelTets() {
  uint64_t j;
  uint64_t tn = mesh->numTets();
  for (uint64_t i = 0; i < Del_deleted.size(); i++) {
    uint64_t to_delete = Del_deleted[i];
    uint64_t lastTet = (--tn) * 4;

    if (mesh->isToDelete(lastTet)) {
      for (j = i; j < Del_deleted.size(); j++)
        if (Del_deleted[j] == lastTet)
          break;

      Del_deleted[j] = Del_deleted[i];
    } else {
      for (j = 0; j < 4; j++) {
        mesh->tet_node[to_delete + j] = mesh->tet_node[lastTet + j];

        uint64_t neigh = mesh->tet_neigh[lastTet + j];
        mesh->tet_neigh[to_delete + j] = neigh;
        mesh->tet_neigh[neigh] = to_delete + j;

        if (mesh->tet_node[lastTet + j] != INFINITE_VERTEX &&
            mesh->inc_tet[mesh->tet_node[lastTet + j]] == lastTet >> 2)
          mesh->inc_tet[mesh->tet_node[lastTet + j]] = to_delete >> 2;
      }
      mesh->mark_tetrahedra[to_delete >> 2] =
          mesh->mark_tetrahedra[lastTet >> 2];
    }
  }
  mesh->resizeTets(tn);
  Del_deleted.clear();
}
#endif

int TetMeshTetrahedrizer::symbolicPerturbation(uint32_t indices[5]) const {
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

int TetMeshTetrahedrizer::vertexInTetSphere(const uint32_t Node[4],
                                            uint32_t v_id) const {
  int det = vInSphere(Node[0], Node[1], Node[2], Node[3], v_id);
  if (det)
    return det;
  uint32_t nn[5] = {Node[0], Node[1], Node[2], Node[3], v_id};
  det = symbolicPerturbation(nn);
  if (det == 0.0)
    ip_error("Symbolic perturbation failed! Should not happen.\n");
  return det;
}

int TetMeshTetrahedrizer::vertexInTetSphere(uint64_t tet, uint32_t v_id) const {
  const uint32_t *Node = mesh->getTetNodes(tet);
  int det;

  if (Node[3] == INFINITE_VERTEX) {
    if ((det = vOrient3D(Node[0], Node[1], Node[2], v_id)) != 0)
      return det;
    const uint32_t nn[4] = {Node[0], Node[1], Node[2],
                            mesh->tet_node[mesh->tet_neigh[tet + 3]]};
    return -vertexInTetSphere(nn, v_id);
  } else
    return vertexInTetSphere(Node, v_id);
}

#ifdef USE_MAROTS_METHOD
void TetMeshTetrahedrizer::deleteInSphereTets(uint64_t tet,
                                              const uint32_t v_id) {
  pushAndMarkDeletedTets(tet);

  for (uint64_t t = Del_deleted.size() - 1; t < Del_deleted.size(); t++) {
    uint64_t tet = Del_deleted[t];
    uint64_t *Neigh = mesh->getTetNeighs(tet);
    uint32_t *Node = mesh->getTetNodes(tet);

    uint64_t neigh = mesh->getIthNeighbor(Neigh, 0);
    if (!mesh->isToDelete(neigh)) {
      if (vertexInTetSphere(neigh, v_id) < 0)
        bnd_push(v_id, Node[1], Node[2], Node[3], Neigh[0]);
      else
        pushAndMarkDeletedTets(neigh);
    }

    neigh = mesh->getIthNeighbor(Neigh, 1);
    if (!mesh->isToDelete(neigh)) {
      if (vertexInTetSphere(neigh, v_id) < 0)
        bnd_push(v_id, Node[2], Node[0], Node[3], Neigh[1]);
      else
        pushAndMarkDeletedTets(neigh);
    }

    neigh = mesh->getIthNeighbor(Neigh, 2);
    if (!mesh->isToDelete(neigh)) {
      if (vertexInTetSphere(neigh, v_id) < 0)
        bnd_push(v_id, Node[0], Node[1], Node[3], Neigh[2]);
      else
        pushAndMarkDeletedTets(neigh);
    }

    neigh = mesh->getIthNeighbor(Neigh, 3);
    if (!mesh->isToDelete(neigh)) {
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

void TetMeshTetrahedrizer::tetrahedrizeHole(uint64_t *tet) {
  uint64_t clength = Del_deleted.size(); // Num tets removed
  uint64_t blength = numDelTmp();        // Num tets to insert

  uint64_t tn = mesh->numTets();

  if (blength > clength) {
    for (uint64_t i = clength; i < blength; i++, tn++)
      Del_deleted.push_back(tn << 2);

    clength = blength;
    mesh->resizeTets(tn);
  }

  uint64_t start = clength - blength;

  for (uint64_t i = 0; i < blength; i++) {
    const uint64_t tet = Del_deleted[i + start];
    uint32_t *Node = mesh->getTetNodes(tet);

    Node[0] = Del_tmp[i].node[0];
    Node[1] = Del_tmp[i].node[1];
    Node[2] = Del_tmp[i].node[2];
    Node[3] = Del_tmp[i].node[3];

    uint64_t bnd = Del_tmp[i].bnd;
    mesh->tet_neigh[tet] = bnd;
    mesh->tet_neigh[bnd] = tet;
    Del_tmp[i].bnd = tet;

    mesh->mark_tetrahedra[tet >> 2] = 0;

    if (mesh->tet_node[tet + 3] != INFINITE_VERTEX)
      for (uint32_t j = 0; j < 4; j++)
        mesh->inc_tet[mesh->tet_node[tet + j]] = tet >> 2;
  }

  uint64_t tlength = 0;
  const uint64_t middle = blength * 3 / 2;

  uint64_t *Tmp = delTmpVec();
  const unsigned index[4] = {2, 3, 1, 2};

  for (uint64_t i = 0; i < blength; i++) {
    uint64_t tet = Del_deleted[start + i];
    const uint32_t *Node = mesh->getTetNodes(tet);

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
        mesh->tet_neigh[tet] = pairValue;
        mesh->tet_neigh[pairValue] = tet;
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

void TetMeshTetrahedrizer::insertExistingVertex(const uint32_t vi,
                                                uint64_t &ct) {
  ct = searchTetrahedron(ct, vi);
  deleteInSphereTets(ct, vi);
  tetrahedrizeHole(&ct);
  uint64_t lt = ct;
  if (mesh->tet_node[lt + 3] == INFINITE_VERTEX)
    lt = mesh->tet_neigh[lt + 3];
  mesh->inc_tet[vi] = lt >> 2;
}

#else
// Start from c and turn around v1-v2 as long as adjacencies are well defined.
// When an invalid adjacency is found, reinit it and exit.
void TetMeshTetrahedrizer::seekAndSetMutualAdjacency(
    int p_o0, int p_o1, int p_o2, const uint32_t *v, uint64_t c, uint64_t o,
    const uint32_t *tet_node_data, uint64_t *tet_neigh_data) {
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
void TetMeshTetrahedrizer::restoreLocalConnectivty(
    uint64_t c, const uint32_t *tet_node_data, uint64_t *tet_neigh_data) {
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
void TetMeshTetrahedrizer::insertExistingVertex(const uint32_t v_id,
                                                uint64_t &tet) {
  static std::vector<uint64_t>
      cavityCorners; // Static to avoid reallocation on each call
  static const int fi[4][3] = {{2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
  uint32_t *tet_node_data = mesh->tet_node.data();
  uint64_t *tet_neigh_data = mesh->tet_neigh.data();

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
      uint32_t &mtn0 = mesh->mark_tetrahedra[n0];
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
  uint64_t ntb, newpos = mesh->tet_node.size();
  if (cavityCorners.size() > Del_deleted.size()) {
    mesh->resizeTets(mesh->numTets() +
                     (cavityCorners.size() - Del_deleted.size()));
    tet_node_data = mesh->tet_node.data();
    tet_neigh_data = mesh->tet_neigh.data();
  }

  // Create the new tets
  for (const uint64_t c : cavityCorners) {
    mesh->mark_tetrahedra[c >> 2] = 0;
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
      mesh->inc_tet[*cn] = ntb;
      mesh->inc_tet[*(--cn)] = ntb;
      mesh->inc_tet[*(--cn)] = ntb;
      mesh->inc_tet[v_id] = ntb;
    }
    mesh->mark_tetrahedra[ntb] = 0;
  }

  // Restore the connectivity within the cavity
  for (uint64_t c : cavityCorners)
    restoreLocalConnectivty(c, tet_node_data, tet_neigh_data);

  tet = tet_neigh_data[cavityCorners.back()];

  cavityCorners.clear();
}
#endif

uint32_t TetMeshTetrahedrizer::findEncroachingPoint(const uint32_t ep0,
                                                    const uint32_t ep1,
                                                    uint64_t &tet_e) const {
  static std::vector<uint64_t>
      enc_queue; // Static to avoid reallocation upon each call

  // Start collecting tetrahedra incident at the endpoints
  mesh->VT(ep0, enc_queue);

  for (uint64_t j : enc_queue)
    mesh->mark_Tet_1(j);

  const vector3d p0 = mesh->vertices[ep0];
  const vector3d p1 = mesh->vertices[ep1];
  const double eslen = (p0 - p1).sq_length();

  vector3d ep;
  uint32_t enc_pt_i = UINT32_MAX;

  mesh->marked_vertex[ep0] = mesh->marked_vertex[ep1] = 1;

  // Collect all encroaching points while expanding around insphere vertices
  for (uint32_t ti = 0; ti < enc_queue.size(); ti++) {
    const uint64_t tet = enc_queue[ti];
    const uint64_t tb = tet << 2;

    // Check each tet vertex for 'isphereness' and keep track of the one with
    // largest sphere
    const uint32_t *tn = mesh->tet_node.data() + tb;
    for (uint32_t i = 0; i < 4; i++) {
      const uint32_t ui = tn[i];
      if (!mesh->marked_vertex[ui]) {
        const vector3d &pui = mesh->vertices[ui];
        if (((pui - p0).sq_length() + (pui - p1).sq_length()) <= eslen) {
          mesh->marked_vertex[ui] = 1;
          if (enc_pt_i == UINT32_MAX ||
              vector3d::hasLargerSphere(p0, p1, pui, ep)) {
            ep = pui;
            enc_pt_i = ui;
            tet_e = tb;
          }
        } else
          mesh->marked_vertex[ui] = 2;
      }
    }

    const int nvmask[] = {
        (mesh->marked_vertex[tn[0]] == 1), (mesh->marked_vertex[tn[1]] == 1),
        (mesh->marked_vertex[tn[2]] == 1), (mesh->marked_vertex[tn[3]] == 1)};
    const int totmarkeda = nvmask[0] + nvmask[1] + nvmask[2] + nvmask[3];

    // Expand on adjacent tets if at least one common vertex is insphere
    const uint64_t *tg = mesh->tet_neigh.data() + tb;
    for (uint32_t i = 0; i < 4; i++) {
      const uint64_t nc = tg[i];
      const uint64_t n = nc >> 2;
      if (mesh->is_marked_Tet_1(n) == 2 ||
          mesh->tet_node[nc] == INFINITE_VERTEX)
        continue;
      const int totmarked = totmarkeda - nvmask[i];
      if (totmarked) {
        mesh->mark_Tet_1(n);
        enc_queue.push_back(n);
      }
    }
  }

  // Clear all marks
  mesh->marked_vertex[ep0] = mesh->marked_vertex[ep1] = 0;
  for (uint64_t j : enc_queue) {
    mesh->unmark_Tet_1(j);
    j <<= 2;
    mesh->marked_vertex[mesh->tet_node[j++]] = 0;
    mesh->marked_vertex[mesh->tet_node[j++]] = 0;
    mesh->marked_vertex[mesh->tet_node[j++]] = 0;
    mesh->marked_vertex[mesh->tet_node[j]] = 0;
  }
  enc_queue.clear();

  return enc_pt_i;
}
