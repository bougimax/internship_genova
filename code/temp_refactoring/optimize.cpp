#include "optimize.h"
#include "delaunay.h"
using namespace TetMesh;

/// FIRST PASS
std::size_t edge_to_size_t(edge e) {
  return (static_cast<std::size_t>(e.first) << 32) | e.second;
}

edge size_t_to_edge(std::size_t value) {
  unsigned int first = static_cast<unsigned int>(value >> 32);
  unsigned int second = static_cast<unsigned int>(value & 0xFFFFFFFF);
  return std::make_pair(first, second);
}

void TetMeshOptimizer::first_pass() {

  std::cout << "Starting FIRST pass" << std::endl;

  std::vector<edge> edges;
  mesh_->getMeshEdges(edges);
  better_priority_queue::updatable_priority_queue<size_t, double>
      edges_to_split;

  for (edge e : edges) {
    std::unique_ptr<Split_edge_info> split_info = get_split_edge_info(e);
    if (split_info->is_good) {
      edges_to_split.push(edge_to_size_t(e), split_info->delta);
    }
  }

  std::cout << "Determined " << edges_to_split.size() << " edges to split"
            << std::endl;

  int num_splitted_edges = 0;

  while (!edges_to_split.empty()) {
    auto edge = edges_to_split.pop_value();
    if (edge.priority > 0) {
      split_edge_and_update(size_t_to_edge(edge.key), edges_to_split);
      num_splitted_edges++;
    } else
      break;
  }

  std::cout << "Actually splitted " << num_splitted_edges << " edges"
            << std::endl;

  std::cout << "Finished FIRST pass" << std::endl;
}

double TetMeshOptimizer::get_energy_from_splitting(tetrahedra tetrahedra,
                                                   edge edge,
                                                   pointType *split_point) {

  std::vector<pointType *> tetrahedra_points_v1 =
                               mesh_->getTetPoints(tetrahedra),
                           tetrahedra_points_v2 =
                               mesh_->getTetPoints(tetrahedra);
  tetrahedra_points_v1[mesh_->get_index_of_vertex_in_tet(
      edge.first, tetrahedra)] = split_point;
  tetrahedra_points_v2[mesh_->get_index_of_vertex_in_tet(
      edge.second, tetrahedra)] = split_point;

  double energy_tet_1 = tetEnergy(tetrahedra_points_v1);
  double energy_tet_2 = tetEnergy(tetrahedra_points_v2);

  // Prevent inversion
  if (!is_oriented_good(tetrahedra_points_v1) ||
      !is_oriented_good(tetrahedra_points_v2))
    return DBL_MAX;

  return std::max(energy_tet_1, energy_tet_2);
}

pointType *TetMeshOptimizer::get_split_point(edge e) {
  return ((vector3d(vertices[e.first]) + vector3d(vertices[e.second])) * 0.5)
      .toExplicitPoint();
}

std::unique_ptr<TetMeshOptimizer::Split_edge_info>
TetMeshOptimizer::get_split_edge_info(edge e) {

  std::unique_ptr<Split_edge_info> split_info =
      std::make_unique<Split_edge_info>();
  split_info->edge = e;

  std::vector<tetrahedra> incident_tetrahedras;
  mesh_->ETfull(e.first, e.second, incident_tetrahedras);

  double pre_transformation_energy = 0., post_transformation_energy = 0.,
         energy_of_split;

  // implicitPoint_LNC *potential_split_point =
  //     new implicitPoint_LNC(vertices[e.first]->toExplicit3D(),
  //                           vertices[e.second]->toExplicit3D(), 0.5);
  //  TODO : Implicit point not working great

  pointType *split_point = get_split_point(e);

  split_info->split_point = split_point;

  for (tetrahedra tet : incident_tetrahedras) {
    if (!mesh_->has_infinite_vertex(tet) &&
        (mesh_->mark_tetrahedra[tet] == DT_IN || !optimize_only_DT_IN)) {
      energy_of_split = get_energy_from_splitting(tet, e, split_point);
      if (energy_of_split == DBL_MAX) {
        split_info->is_good = false;
        return split_info;
      }
      pre_transformation_energy =
          std::max(pre_transformation_energy, tetrahedras_energy[tet]);
      post_transformation_energy =
          std::max(post_transformation_energy, energy_of_split);
    }
  }

  split_info->is_good = pre_transformation_energy > post_transformation_energy;

  if (split_info->is_good) {
    split_info->delta = pre_transformation_energy - post_transformation_energy;
    split_info->pre_energy = pre_transformation_energy;
    if (mesh_->isOnBoundary(e.first, e.second))
      split_info->prioritize = 1;
  }
  return split_info;
}

void TetMeshOptimizer::split_edge_and_update(
    edge e,
    better_priority_queue::updatable_priority_queue<size_t, double> &queue) {

  pointType *split_point = get_split_point(e);
  mesh_->pushVertex(split_point);

  std::vector<tetrahedra> impacted_tetrahedras;
  split_edge(e, vertices.size() - 1, impacted_tetrahedras);

  std::vector<edge> edges_to_update;
  get_edges_from_tetrahedras(edges_to_update, impacted_tetrahedras);

  for (edge edge_to_update : edges_to_update) {
    std::unique_ptr<Split_edge_info> split_info =
        get_split_edge_info(edge_to_update);
    if (split_info->is_good)
      queue.set(edge_to_size_t(edge_to_update), split_info->delta);
    else
      queue.update(edge_to_size_t(edge_to_update), -1);
  }
}

void TetMeshOptimizer::split_edge(
    edge edge_to_split, vertex split_vertex,
    std::vector<tetrahedra> &impacted_tetrahedras) {
  // Consider edge_to_split := v1 -- v2
  static std::vector<tetrahedra> incident_tetrahedras;
  incident_tetrahedras.clear();
  mesh_->ETfull(edge_to_split.first, edge_to_split.second,
                incident_tetrahedras);

  corner current_corner = mesh_->tet_node.size();
  corner start_corner = current_corner;
  vertex current_vertex;
  corner opposite_of_v1;

  for (tetrahedra tet : incident_tetrahedras) {
    impacted_tetrahedras.push_back(tet);
    mesh_->mark_tetrahedra.push_back(mesh_->mark_tetrahedra[tet]);
    impacted_tetrahedras.push_back(
        mesh_->get_tetrahedra_index_from_corner(current_corner));
    if (!mesh_->has_infinite_vertex(tet))
      inc_tet[split_vertex] = tet;

    opposite_of_v1 =
        mesh_->tet_neigh[get_corner_in_tet(tet, edge_to_split.first)];

    for (int i = 0; i < 4; i++) {
      current_vertex = mesh_->get_i_th_vertex_of_tetrahedra(tet, i);

      // Duplicate tetrahedra structure to add all the tetrahedras in which we
      // replace v1 by split_vertex
      if (current_vertex == edge_to_split.first) {
        // If we are on v1 we push split_vertex
        mesh_->tet_node.push_back(split_vertex);

        // We arrange corners

        mesh_->tet_neigh.push_back(opposite_of_v1);
        mesh_->tet_neigh[opposite_of_v1] = current_corner;

      } else {
        mesh_->tet_node.push_back(current_vertex);
        if (current_vertex == edge_to_split.second) {
          // Change inc_tet for v2 because v2 does not belong anymore to the
          // same tets
          if (!mesh_->has_infinite_vertex(tet))
            inc_tet[edge_to_split.second] =
                mesh_->get_tetrahedra_index_from_corner(current_corner);
          mesh_->tet_neigh.push_back(
              get_corner_in_tet(tet, edge_to_split.first));
          mesh_->tet_neigh[get_corner_in_tet(tet, edge_to_split.first)] =
              current_corner;
        } else {
          // Else reproduce tet structure at the end of mesh_->tet_node and
          // mesh_->tet_neigh
          corner opposite_former_corner =
              mesh_->tet_neigh[mesh_->get_i_th_corner_of_tetrahedra(tet, i)];
          tetrahedra tet_of_opp_former_corner =
              mesh_->get_tetrahedra_index_from_corner(opposite_former_corner);
          auto index_of_processing_tet =
              std::find(incident_tetrahedras.begin(),
                        incident_tetrahedras.end(), tet_of_opp_former_corner) -
              incident_tetrahedras.begin();

          uint32_t index_opposite_former_corner =
              get_index_of_corner_in_tet(opposite_former_corner);
          corner corner_to_add = start_corner + (index_of_processing_tet << 2) +
                                 index_opposite_former_corner;
          mesh_->tet_neigh.push_back(corner_to_add);
        }
      }
      current_corner++;
    }
    // Finally we replace v2 by split_vertex in the old tetrahedra, aside from
    // changing the opposite corner for v2 which we did in previous loop there
    // is nothing to change
    mesh_->tet_node[get_corner_in_tet(tet, edge_to_split.second)] =
        split_vertex;

    // Updating energy
    tetrahedras_energy.push_back(getTetEnergy(
        mesh_->get_tetrahedra_index_from_corner(current_corner - 4)));
    tetrahedras_energy[tet] = getTetEnergy(tet);
  }
}

/// SECOND PASS

void TetMeshOptimizer::second_pass() {

  std::cout << "\nStarting SECOND pass" << std::endl;

  std::vector<edge> edges;
  getMeshEdges(edges);
  std::set<vertex> deleted_vertices;
  better_priority_queue::updatable_priority_queue<size_t, double>
      edges_to_collapse;

  for (edge e : edges) {
    std::unique_ptr<Collapse_info> collapse_info = get_collapse_info(e);
    if (collapse_info->is_good) {
      edges_to_collapse.push(edge_to_size_t(collapse_info->edge),
                             collapse_info->delta);
    }
  }

  size_t num_edge_collapsed = 0;

  std::cout << "Determined " << edges_to_collapse.size() << " edges to collapse"
            << std::endl;

  while (!edges_to_collapse.empty()) {
    auto edge = edges_to_collapse.pop_value();
    if (edge.priority > 0) {
      collapse_on_v1_and_update(size_t_to_edge(edge.key), edges_to_collapse);
      num_edge_collapsed++;
    } else
      break;
  }

  std::cout << "Actually collapsed " << num_edge_collapsed << " edges"
            << std::endl;

  std::cout << "Finished SECOND pass" << std::endl;

  temp_remap_vertex.clear();
}

void TetMeshOptimizer::collapse_on_v1_and_update(
    edge edge_to_collapse,
    better_priority_queue::updatable_priority_queue<size_t, double> &queue) {
  std::vector<tetrahedra> impacted_tetrahedras;
  std::vector<edge> removed_edges;
  collapse_on_v1(edge_to_collapse, impacted_tetrahedras, removed_edges);

  // Remove edges that have been suppressed
  for (edge e : removed_edges) {
    queue.update(edge_to_size_t(e), -1);
  }

  std::vector<edge> edges_to_update;
  get_edges_from_tetrahedras(edges_to_update, impacted_tetrahedras);

  // Updates edges
  for (edge e : edges_to_update) {
    std::unique_ptr<Collapse_info> collapse_info = get_collapse_info(e);
    if (collapse_info->is_good) {
      queue.update(edge_to_size_t(std::make_pair((collapse_info->edge).second,
                                                 (collapse_info->edge).first)),
                   -1);
      queue.set(edge_to_size_t(collapse_info->edge), collapse_info->delta);
    } else {
      queue.update(edge_to_size_t(e), -1);
      queue.update(edge_to_size_t(std::make_pair(e.second, e.first)), -1);
    }
  }
}

void TetMeshOptimizer::mark_vertices(std::vector<vertex> &vertices_to_mark,
                                     unsigned char mark) {
  for (vertex v : vertices_to_mark)
    mesh_->marked_vertex[v] = mark;
}

bool TetMeshOptimizer::link_condition(edge e) {
  std::vector<vertex> one_ring_v1, one_ring_v2, one_ring_e,
      intersection_one_ring;
  std::vector<tetrahedra> incident_tetrahedras_e;
  edge opposite_edge;

  VV(e.first, one_ring_v1);
  VV(e.second, one_ring_v2);
  OneRing(e, one_ring_e, incident_tetrahedras_e);

  size_t size_intersection = 0;

  mark_vertices(one_ring_v2, 0);
  mark_vertices(one_ring_v1, 1);

  for (vertex neighbour_v2 : one_ring_v2) {
    if (mesh_->marked_vertex[neighbour_v2] == 1) {
      mesh_->marked_vertex[neighbour_v2] |= 2;
      size_intersection++;
    }
  }

  for (vertex neighbour_e : one_ring_e) {
    mesh_->marked_vertex[neighbour_e] |= 4;
    if ((mesh_->marked_vertex[neighbour_e] ^ 7) == 0) {
      size_intersection--;
    } else
      return false;
  }

  mark_vertices(one_ring_v2, 0);
  mark_vertices(one_ring_v1, 0);

  return size_intersection == 0;
}

double TetMeshOptimizer::get_energy_from_collapsing(
    edge e, std::vector<tetrahedra> &incident_tetrahedras_v2) {
  auto [v1, v2] = e;
  double post_transformation_energy = 0.;

  for (tetrahedra t2 : incident_tetrahedras_v2) {
    if (!mesh_->has_infinite_vertex(t2) &&
        (mesh_->mark_tetrahedra[t2] == DT_IN || !optimize_only_DT_IN) &&
        !mesh_->tetHasVertex(t2, v1)) {

      std::vector<pointType *> t2_points = mesh_->getTetPoints(t2);
      t2_points[mesh_->get_index_of_vertex_in_tet(v2, t2)] = vertices[v1];

      if (!is_oriented_good(t2_points))
        return DBL_MAX;

      post_transformation_energy =
          std::max(post_transformation_energy, tetEnergy(t2_points));
    }
  }

  return post_transformation_energy;
}

std::unique_ptr<TetMeshOptimizer::Collapse_info>
TetMeshOptimizer::get_collapse_info(edge &e) {

  std::unique_ptr<Collapse_info> collapse_info =
      std::make_unique<Collapse_info>();
  // Consider e := v1 -- v2
  auto [v1, v2] = e;

  if (!link_condition(e) || mesh_->isOnBoundary(v1, v2)) {
    collapse_info->is_good = false;
    return collapse_info;
  }

  std::vector<tetrahedra> incident_tetrahedras_v1, incident_tetrahedras_v2;
  VTfull(v1, incident_tetrahedras_v1);
  VTfull(v2, incident_tetrahedras_v2);

  double pre_transformation_energy = 0., post_transformation_energy_v1,
         post_transformation_energy_v2, post_transformation_energy;

  // Compute pre-tranformation energy
  for (tetrahedra t : incident_tetrahedras_v1) {
    if (!mesh_->has_infinite_vertex(t) &&
        (mesh_->mark_tetrahedra[t] == DT_IN || !optimize_only_DT_IN) &&
        mesh_->tetHasVertex(t, v2))
      pre_transformation_energy =
          std::max(pre_transformation_energy, tetrahedras_energy[t]);
  }

  post_transformation_energy_v1 =
      get_energy_from_collapsing(e, incident_tetrahedras_v2);
  post_transformation_energy_v2 = get_energy_from_collapsing(
      std::make_pair(v2, v1), incident_tetrahedras_v1);

  post_transformation_energy =
      std::min(post_transformation_energy_v1, post_transformation_energy_v2);

  collapse_info->is_good =
      post_transformation_energy < pre_transformation_energy;

  if (collapse_info->is_good) {
    collapse_info->delta =
        pre_transformation_energy - post_transformation_energy;
    collapse_info->pre_energy = pre_transformation_energy;
    if (post_transformation_energy_v1 < post_transformation_energy_v2 &&
        !mesh_->isOnBoundary(v2)) {
      collapse_info->edge = e;
    } else if (post_transformation_energy_v2 < post_transformation_energy_v1 &&
               !mesh_->isOnBoundary(v1)) {
      collapse_info->edge = std::make_pair(v2, v1);
    } else {
      collapse_info->is_good = false;
    }
  }
  return collapse_info;
}

bool TetMeshOptimizer::collapse_on_v1(
    edge e, std::vector<tetrahedra> &impacted_tetrahedras,
    std::vector<edge> &removed_edges) {
  auto [v1, v2] = e;

  if (v1 == INFINITE_VERTEX || v2 == INFINITE_VERTEX ||
      mesh_->isOnBoundary(v2) || !hasEdge(v1, v2)) {
    return false;
  }

  std::vector<tetrahedra> incident_tetrahedras, incident_tetrahedras_v2;
  std::vector<vertex> neighbour_v2;
  mesh_->ETfull(v1, v2, incident_tetrahedras);
  VTfull(v2, incident_tetrahedras_v2);
  VV(v2, neighbour_v2);

  for (vertex neighbour : neighbour_v2) {
    removed_edges.push_back(std::make_pair(v2, neighbour));
    removed_edges.push_back(std::make_pair(neighbour, v2));
  }

  edge opposite_edge;
  corner opposite_corner_1, opposite_corner_2, potential_corner_to_change;
  tetrahedra t1, t2, next_to_delete;

  std::vector<tetrahedra> tet_to_delete;

  if (incident_tetrahedras.empty()) {
    return false;
  }

  for (tetrahedra t : incident_tetrahedras) {
    mesh_->oppositeTetEdgePair(t, std::make_pair(v1, v2), opposite_edge);
    opposite_corner_1 = mesh_->tet_neigh[get_corner_in_tet(t, v2)];
    opposite_corner_2 = mesh_->tet_neigh[get_corner_in_tet(t, v1)];
    t1 = mesh_->get_tetrahedra_index_from_corner(opposite_corner_1);
    t2 = mesh_->get_tetrahedra_index_from_corner(opposite_corner_2);

    setMutualNeighbors(opposite_corner_1, opposite_corner_2);
    if (!isGhost(t1)) {
      inc_tet[v1] = t1;
    }

    inc_tet[v2] = UINT64_MAX;
    if (opposite_edge.first != INFINITE_VERTEX && !isGhost(t1)) {
      inc_tet[opposite_edge.first] = t1;
      inc_tet[opposite_edge.second] = t1;
    }
    if (opposite_edge.second != INFINITE_VERTEX && !isGhost(t2)) {
      inc_tet[opposite_edge.first] = t2;
      inc_tet[opposite_edge.second] = t2;
    }
    tet_to_delete.push_back(t);
  }

  for (tetrahedra t : incident_tetrahedras_v2) {
    potential_corner_to_change = get_corner_in_tet(t, v2);
    if (potential_corner_to_change != UINT64_MAX) {
      mesh_->tet_node[potential_corner_to_change] = v1;
      // Update energy for t2
      tetrahedras_energy[t] = getTetEnergy(t);
    }
  }

  // Delete in descendant order (it could be possible to swap t with a tet at
  // the end that is also to delete)

  std::sort(tet_to_delete.begin(), tet_to_delete.end());

  while (!tet_to_delete.empty()) {
    next_to_delete = tet_to_delete.back();
    tet_to_delete.pop_back();
    remove_tetrahedra(next_to_delete);
  }

  // Add the edges that have been impacted by removing v2 which implies to remap
  // the last vertex to v2, in this case all the edges incident to last_vertex
  // are now incident to v2
  std::vector<vertex> impacted_vertex_remap;
  vertex last_vertex = mesh_->numVertices() - 1;
  if (v2 != last_vertex) {
    VV(last_vertex, impacted_vertex_remap);
    VTfull(last_vertex, impacted_tetrahedras);
    for (vertex neighbour : impacted_vertex_remap) {
      removed_edges.push_back(std::make_pair(last_vertex, neighbour));
      removed_edges.push_back(std::make_pair(neighbour, last_vertex));
    }
  }

  remove_vertex(v2);

  remap_vertex(v1);

  VTfull(v1, impacted_tetrahedras);
  temp_remap_vertex.clear();

  return true;
}

/// THIRD PASS

void TetMeshOptimizer::third_pass() {

  // To iterate over all faces once we can represent uniquely a face with the
  // two opposite corner they induce this way we iterate over all corners then
  // we are considering the corner only if it's smaller than its
  // mesh_->tet_neigh then we add those pair of corner in a queue for swapping
  // if it's worth it, then when we are popping from the queue we just verify if
  // both corner are still opposed (i.e. they haven't been touched by an other
  // swap)

  std::cout << "\nStarting THIRD pass" << std::endl;

  third_pass_edge();

  third_pass_face();

  std::cout << "Finished THIRD pass" << std::endl;
}

void TetMeshOptimizer::third_pass_face() {
  better_priority_queue::updatable_priority_queue<corner, double> faces_to_swap;

  for (corner face = 0; face < mesh_->tet_node.size(); face++) {
    std::unique_ptr<Swap_face_info> swap_info = get_swap_face_info(face);
    if (swap_info->is_good)
      faces_to_swap.push(face, swap_info->delta);
  }

  std::cout << "Determined " << faces_to_swap.size() << " faces to swap"
            << std::endl;

  uint32_t num_face_swapped = 0;

  while (!faces_to_swap.empty()) {
    auto act_info = faces_to_swap.pop_value();
    if (act_info.priority > 0) {
      swap_face_and_update(act_info.key, faces_to_swap);
      num_face_swapped++;
    } else
      break;
  }

  std::cout << "Actually swapped " << num_face_swapped << " faces\n"
            << std::endl;
}

void TetMeshOptimizer::third_pass_edge() {
  std::vector<edge> edges;
  getMeshEdges(edges);
  better_priority_queue::updatable_priority_queue<size_t,
                                                  std::pair<double, vertex>>
      edges_to_swap;

  for (edge e : edges) {
    std::unique_ptr<Swap_edge_info> swap_info = get_swap_edge_info(e);
    if (swap_info->is_good)
      edges_to_swap.push(
          edge_to_size_t(e),
          std::make_pair(swap_info->delta, swap_info->collapse_vertex));
  }

  std::cout << "Determined " << edges_to_swap.size() << " edges to swap"
            << std::endl;

  int num_edges_swapped = 0;

  while (!edges_to_swap.empty()) {
    auto act_info = edges_to_swap.pop_value();
    if (act_info.priority.first > 0) {
      swap_edge_and_update(size_t_to_edge(act_info.key),
                           act_info.priority.second, edges_to_swap);
      num_edges_swapped++;
    } else
      break;
  }
  std::cout << "Actually swapped " << num_edges_swapped << " edges"
            << std::endl;
}

void TetMeshOptimizer::swap_face_and_update(
    corner face_to_swap,
    better_priority_queue::updatable_priority_queue<corner, double> &queue) {
  std::vector<corner> impacted_faces;
  swap_face(face_to_swap, impacted_faces);

  for (corner face : impacted_faces) {
    std::unique_ptr<Swap_face_info> swap_info = get_swap_face_info(face);
    if (swap_info->is_good) {
      queue.set(face, swap_info->delta);
      queue.update(mesh_->tet_neigh[face], -1);
    } else {
      queue.update(face, -1);
    }
  }
}
void TetMeshOptimizer::swap_edge_and_update(
    edge edge_to_swap, vertex collapse_vertex,
    better_priority_queue::updatable_priority_queue<
        size_t, std::pair<double, vertex>> &queue) {
  std::vector<edge> impacted_edges;
  swap_edge(edge_to_swap, collapse_vertex, impacted_edges);

  for (edge e : impacted_edges) {
    std::unique_ptr<Swap_edge_info> swap_info = get_swap_edge_info(e);
    if (swap_info->is_good) {
      queue.set(edge_to_size_t(e),
                std::make_pair(swap_info->delta, swap_info->collapse_vertex));
      queue.update(edge_to_size_t(std::make_pair(e.second, e.first)),
                   std::make_pair(-1, 0));
    } else {
      queue.update(edge_to_size_t(e), std::make_pair(-1, 0));
      queue.update(edge_to_size_t(std::make_pair(e.second, e.first)),
                   std::make_pair(-1, 0));
    }
  }
}

void TetMeshOptimizer::swap_edge(edge edge_to_swap, vertex collapse_vertex,
                                 std::vector<edge> &impacted_edges) {
  std::vector<vertex> one_ring;
  std::vector<tetrahedra> incident_tetrahedras;
  std::vector<vertex> neigh_collapse;
  edge opposite_edge;

  OneRing(edge_to_swap, one_ring, incident_tetrahedras);

  // Add impacted edges that are connected to one of the two endpoint of the
  // edge to swap, then also add newly created edges
  for (vertex neigh : one_ring) {
    impacted_edges.push_back(std::make_pair(edge_to_swap.first, neigh));
    impacted_edges.push_back(std::make_pair(edge_to_swap.second, neigh));
    if (neigh != collapse_vertex) {
      impacted_edges.push_back(std::make_pair(neigh, collapse_vertex));
    }
  }

  // Add impacted edges that are around the edge we are going to swap
  for (tetrahedra t : incident_tetrahedras) {
    mesh_->oppositeTetEdgePair(t, edge_to_swap, opposite_edge);
    impacted_edges.push_back(opposite_edge);
  }

  // Do the swapping
  pointType *split_point = get_split_point(edge_to_swap);
  mesh_->pushVertex(split_point);

  std::vector<tetrahedra> impacted_tetrahedras;
  split_edge(edge_to_swap, mesh_->numVertices() - 1, impacted_tetrahedras);

  std::vector<edge> removed_edges;
  collapse_on_v1(std::make_pair(collapse_vertex, mesh_->numVertices() - 1),
                 impacted_tetrahedras, removed_edges);
}

std::unique_ptr<TetMeshOptimizer::Swap_face_info>
TetMeshOptimizer::get_swap_face_info(corner face) {
  std::unique_ptr<Swap_face_info> swap_info =
      std::make_unique<Swap_face_info>();
  double pre_transformation_energy = 0, post_transformation_energy = DBL_MAX;

  if (face < mesh_->tet_neigh[face] &&
      !mesh_->has_infinite_vertex(
          mesh_->get_tetrahedra_index_from_corner(face)) &&
      mesh_->tet_node[mesh_->tet_neigh[face]] != INFINITE_VERTEX &&
      ((mesh_->mark_tetrahedra[mesh_->get_tetrahedra_index_from_corner(face)] ==
            DT_IN &&
        mesh_->mark_tetrahedra[mesh_->get_tetrahedra_index_from_corner(
            mesh_->tet_neigh[face])] == DT_IN) ||
       (mesh_->mark_tetrahedra[mesh_->get_tetrahedra_index_from_corner(face)] ==
            DT_OUT &&
        mesh_->mark_tetrahedra[mesh_->get_tetrahedra_index_from_corner(
            mesh_->tet_neigh[face])] == DT_OUT &&
        !optimize_only_DT_IN))) {
    pre_transformation_energy = std::max(
        tetrahedras_energy[mesh_->get_tetrahedra_index_from_corner(face)],
        tetrahedras_energy[mesh_->get_tetrahedra_index_from_corner(
            mesh_->tet_neigh[face])]);
    post_transformation_energy = get_energy_from_swapping_face(face);
  }

  swap_info->is_good = post_transformation_energy < pre_transformation_energy;

  if (swap_info->is_good) {
    swap_info->delta = pre_transformation_energy - post_transformation_energy;
    swap_info->pre_energy = pre_transformation_energy;
    swap_info->face = std::make_pair(face, mesh_->tet_neigh[face]);
  }

  return swap_info;
}

double TetMeshOptimizer::get_energy_from_swapping_face(corner face) {

  const uint64_t index_in_tet = face & 3;
  const corner tet_basis = face - index_in_tet;

  const uint64_t r0 = tet_basis + tetON1(index_in_tet),
                 r1 = tet_basis + tetON3(index_in_tet),
                 r2 = tet_basis + tetON2(index_in_tet);
  const uint32_t c0 = mesh_->tet_node[r0], c1 = mesh_->tet_node[r1],
                 c2 = mesh_->tet_node[r2], c3 = mesh_->tet_node[face];

  const uint64_t g00 = mesh_->tet_neigh[r0], g01 = mesh_->tet_neigh[r1],
                 g02 = mesh_->tet_neigh[r2];

  const uint64_t orx = mesh_->tet_neigh[face];
  const uint64_t opp = orx & (~3);
  const uint64_t or0 = tetCornerAtVertex(opp, c0);
  const uint64_t or1 = tetCornerAtVertex(opp, c1);
  const uint64_t or2 = tetCornerAtVertex(opp, c2);

  const uint64_t g10 = mesh_->tet_neigh[or0], g11 = mesh_->tet_neigh[or1],
                 g12 = mesh_->tet_neigh[or2];

  const uint32_t oc = mesh_->tet_node[orx];

  return std::max(
      tetEnergy(vertices[c3], vertices[oc], vertices[c1], vertices[c2]),
      std::max(
          tetEnergy(vertices[c3], vertices[c0], vertices[oc], vertices[c2]),
          tetEnergy(vertices[c3], vertices[c0], vertices[c1], vertices[oc])));
}

std::unique_ptr<TetMeshOptimizer::Swap_edge_info>
TetMeshOptimizer::get_swap_edge_info(edge e) {

  std::unique_ptr<Swap_edge_info> swap_info =
      std::make_unique<Swap_edge_info>();
  swap_info->edge = e;

  if (mesh_->isOnBoundary(e.first, e.second)) {
    swap_info->is_good = false;
    return swap_info;
  }

  std::vector<tetrahedra> incident_tetrahedras;
  std::vector<vertex> one_ring;

  OneRing(e, one_ring, incident_tetrahedras);

  double energy_of_swapping = DBL_MAX, current_energy_of_swapping,
         pre_transformation_energy = 0;
  edge opposite_edge;
  vertex min_collapse_point;
  uint32_t i_v1, i_v2, num_tet_seen;

  for (vertex collapse_point : one_ring) {
    num_tet_seen = 0;
    current_energy_of_swapping = 0;

    for (tetrahedra t : incident_tetrahedras) {
      if (mesh_->has_infinite_vertex(t)) {
        current_energy_of_swapping = DBL_MAX;
        break;
      }
      num_tet_seen++;
      pre_transformation_energy =
          std::max(tetrahedras_energy[t], pre_transformation_energy);
      mesh_->oppositeTetEdgePair(t, e, opposite_edge);
      if (opposite_edge.first != collapse_point &&
          opposite_edge.second != collapse_point) {
        std::vector<pointType *> t_v1_points = mesh_->getTetPoints(t);
        std::vector<pointType *> t_v2_points = t_v1_points;
        i_v1 = mesh_->get_index_of_vertex_in_tet(e.first, t);
        t_v1_points[i_v1] = vertices[collapse_point];
        i_v2 = mesh_->get_index_of_vertex_in_tet(e.second, t);
        t_v2_points[i_v2] = vertices[collapse_point];
        if (!is_oriented_good(t_v1_points) || !is_oriented_good(t_v2_points)) {
          current_energy_of_swapping = DBL_MAX;
          break;
        } else {
          current_energy_of_swapping = std::max(
              current_energy_of_swapping,
              std::max(tetEnergy(t_v1_points), tetEnergy(t_v2_points)));
        }
      }
    }

    if (num_tet_seen <= 2) {
      current_energy_of_swapping = DBL_MAX;
    }

    if (current_energy_of_swapping < energy_of_swapping) {
      energy_of_swapping = current_energy_of_swapping;
      min_collapse_point = collapse_point;
    }
  }

  swap_info->is_good = pre_transformation_energy > energy_of_swapping;

  if (swap_info->is_good) {
    swap_info->delta = pre_transformation_energy - energy_of_swapping;
    swap_info->pre_energy = pre_transformation_energy;
    swap_info->collapse_vertex = min_collapse_point;
    if (mesh_->isOnBoundary(e.first) && mesh_->isOnBoundary(e.second))
      swap_info->prioritize = 1;
  }

  return swap_info;
}

// 2-3 swap
// See blender file to understand indices (it's non sense otherwise)
// When it says a face it's a corner that represent the opposite face in its
// tetrahedra
bool TetMeshOptimizer::swap_face(uint64_t face,
                                 std::vector<corner> &impacted_faces,
                                 bool prevent_inversion, double th_energy) {
  const uint64_t b2 = mesh_->tet_node.size();
  const size_t newsize = mesh_->tet_node.size() + 4;

  const uint64_t index_in_tet = face & 3;
  const corner tet_basis = face - index_in_tet;

  const uint64_t r0 = tet_basis + tetON1(index_in_tet),
                 r1 = tet_basis + tetON3(index_in_tet),
                 r2 = tet_basis + tetON2(index_in_tet);
  const uint32_t c0 = mesh_->tet_node[r0], c1 = mesh_->tet_node[r1],
                 c2 = mesh_->tet_node[r2], c3 = mesh_->tet_node[face];

  const uint64_t g00 = mesh_->tet_neigh[r0], g01 = mesh_->tet_neigh[r1],
                 g02 = mesh_->tet_neigh[r2];

  const uint64_t orx = mesh_->tet_neigh[face];
  const uint64_t opp = orx & (~3);
  const uint64_t or0 = tetCornerAtVertex(opp, c0);
  const uint64_t or1 = tetCornerAtVertex(opp, c1);
  const uint64_t or2 = tetCornerAtVertex(opp, c2);

  const uint64_t g10 = mesh_->tet_neigh[or0], g11 = mesh_->tet_neigh[or1],
                 g12 = mesh_->tet_neigh[or2];

  const uint32_t oc = mesh_->tet_node[orx];

  if (prevent_inversion) {
    // Verify that the swap does not invert any tet
    if (vOrient3D(c3, oc, c1, c2) <= 0 || vOrient3D(c3, c0, oc, c2) <= 0 ||
        vOrient3D(c3, c0, c1, oc) <= 0)
      return false;
  }

  mesh_->tet_node.resize(newsize);
  mesh_->tet_neigh.resize(newsize);
  mesh_->mark_tetrahedra.resize(newsize >> 2,
                                mesh_->mark_tetrahedra[tet_basis >> 2]);

  uint32_t *tn = getTetNodes(tet_basis);
  *tn++ = c3;
  *tn++ = oc;
  *tn++ = c1;
  *tn++ = c2;
  tetrahedras_energy[mesh_->get_tetrahedra_index_from_corner(tet_basis)] =
      getTetEnergy(mesh_->get_tetrahedra_index_from_corner(tet_basis));
  tn = getTetNodes(opp);
  *tn++ = c3;
  *tn++ = c0;
  *tn++ = oc;
  *tn++ = c2;
  tetrahedras_energy[mesh_->get_tetrahedra_index_from_corner(opp)] =
      getTetEnergy(mesh_->get_tetrahedra_index_from_corner(opp));
  tn = getTetNodes(b2);
  *tn++ = c3;
  *tn++ = c0;
  *tn++ = c1;
  *tn++ = oc;
  tetrahedras_energy.push_back(
      getTetEnergy(mesh_->get_tetrahedra_index_from_corner(b2)));

  uint64_t *tg = getTetNeighs(tet_basis);
  *tg++ = g10;
  *tg++ = g00;
  *tg++ = opp + 1;
  *tg++ = b2 + 1;
  tg = getTetNeighs(opp);
  *tg++ = g11;
  *tg++ = tet_basis + 2;
  *tg++ = g01;
  *tg++ = b2 + 2;
  tg = getTetNeighs(b2);
  *tg++ = g12;
  *tg++ = tet_basis + 3;
  *tg++ = opp + 3;
  *tg++ = g02;

  mesh_->tet_neigh[g00] = tet_basis + 1;
  mesh_->tet_neigh[g01] = opp + 2;
  mesh_->tet_neigh[g02] = b2 + 3;
  mesh_->tet_neigh[g10] = tet_basis;
  mesh_->tet_neigh[g11] = opp;
  mesh_->tet_neigh[g12] = b2;

  impacted_faces.push_back(g00);
  impacted_faces.push_back(g01);
  impacted_faces.push_back(g02);
  impacted_faces.push_back(g10);
  impacted_faces.push_back(g11);
  impacted_faces.push_back(g12);
  impacted_faces.push_back(tet_basis + 1);
  impacted_faces.push_back(opp + 2);
  impacted_faces.push_back(b2 + 3);
  impacted_faces.push_back(tet_basis);
  impacted_faces.push_back(opp);
  impacted_faces.push_back(b2);

  inc_tet[c0] = opp >> 2;
  inc_tet[c1] = tet_basis >> 2;

  return true;
}

/// FOURTH PASS

void TetMeshOptimizer::fourth_pass() {
  std::cout << "\nStarting FOURTH pass" << std::endl;

  uint32_t n_vertices = mesh_->numVertices();

  better_priority_queue::updatable_priority_queue<vertex, double>
      points_to_move;

  for (vertex v = 0; v < n_vertices; v++) {
    std::unique_ptr<Move_info> move_info = get_move_info(v);
    if (move_info->is_good)
      points_to_move.push(v, move_info->delta);
  }

  std::cout << "Determined " << points_to_move.size() << " points to move"
            << std::endl;

  size_t num_vertex_moved = 0;

  while (!points_to_move.empty()) {
    auto act_info = points_to_move.pop_value();
    if (act_info.priority > 0) {
      move_and_update(act_info.key, points_to_move);
      num_vertex_moved++;
    } else
      break;
  }

  std::cout << "Actually moved " << num_vertex_moved << " points" << std::endl;

  std::cout << "Finished FOURTH pass" << std::endl;
}

void TetMeshOptimizer::move_and_update(
    vertex v,
    better_priority_queue::updatable_priority_queue<vertex, double> &queue) {
  std::vector<vertex> impacted_vertices;
  std::vector<tetrahedra> incident_tetrahedras;

  VV(v, impacted_vertices);
  VT(v, incident_tetrahedras);

  move_vertex(v, get_barycenter(v, incident_tetrahedras));

  for (vertex vertex_to_update : impacted_vertices) {
    std::unique_ptr<Move_info> move_info = get_move_info(vertex_to_update);
    if (move_info->is_good)
      queue.set(vertex_to_update, move_info->delta);
    else
      queue.update(vertex_to_update, -1);
  }
}

pointType *TetMeshOptimizer::get_barycenter(
    vertex v, std::vector<tetrahedra> &incident_tetrahedras) {
  vector3d barycenter = vector3d(0., 0., 0.);
  uint32_t num_neighbour = 0;
  for (tetrahedra tet : incident_tetrahedras) {
    for (int i = 0; i < 4; i++) {
      if (mesh_->get_i_th_vertex_of_tetrahedra(tet, i) != INFINITE_VERTEX)
        mesh_->marked_vertex[mesh_->get_i_th_vertex_of_tetrahedra(tet, i)] = 0;
    }
  }

  for (tetrahedra tet : incident_tetrahedras) {
    if (mesh_->mark_tetrahedra[tet] == DT_IN) {
      for (int i = 0; i < 4; i++) {
        vertex act_vertex = mesh_->get_i_th_vertex_of_tetrahedra(tet, i);
        if (act_vertex != v && mesh_->marked_vertex[act_vertex] == 0) {
          mesh_->marked_vertex[act_vertex] = 1;
          barycenter += vector3d(vertices[act_vertex]);
          num_neighbour++;
        }
      }
    }
  }

  if (num_neighbour < 2)
    throw std::runtime_error("Can't compute barycenter for a vertex that have "
                             "less that 2 neighbour");
  barycenter *= (1. / (double)num_neighbour);

  for (tetrahedra tet : incident_tetrahedras) {
    for (int i = 0; i < 4; i++) {
      if (mesh_->get_i_th_vertex_of_tetrahedra(tet, i) != INFINITE_VERTEX)
        mesh_->marked_vertex[mesh_->get_i_th_vertex_of_tetrahedra(tet, i)] = 0;
    }
  }

  return barycenter.toExplicitPoint();
}

std::unique_ptr<TetMeshOptimizer::Move_info>
TetMeshOptimizer::get_move_info(vertex v) {
  std::unique_ptr<Move_info> move_info = std::make_unique<Move_info>();
  move_info->vertex = v;

  if (mesh_->isOnBoundary(v)) {
    move_info->is_good = false;
    return move_info;
  }

  std::vector<tetrahedra> incident_tetrahedras;
  VT(v, incident_tetrahedras);

  pointType *barycenter;
  try {
    barycenter = get_barycenter(v, incident_tetrahedras);
    move_info->barycenter = barycenter;
  } catch (const std::exception &e) {
    move_info->is_good = false;
    return move_info;
  }

  double pre_transformation_energy = 0., post_transformation_energy = 0.;

  for (tetrahedra tet : incident_tetrahedras) {
    if (mesh_->mark_tetrahedra[tet] == DT_OUT ||
        mesh_->has_infinite_vertex(tet))
      continue;
    pre_transformation_energy =
        std::max(pre_transformation_energy, tetrahedras_energy[tet]);
  }

  for (tetrahedra tet : incident_tetrahedras) {
    if (mesh_->mark_tetrahedra[tet] == DT_OUT ||
        mesh_->has_infinite_vertex(tet))
      continue;
    std::vector<pointType *> tet_points = mesh_->getTetPoints(tet);
    tet_points[mesh_->get_index_of_vertex_in_tet(v, tet)] = barycenter;

    if (!is_oriented_good(tet_points)) {
      move_info->is_good = false;
      return move_info;
    }
    post_transformation_energy =
        std::max(post_transformation_energy, tetEnergy(tet_points));
  }

  move_info->is_good = pre_transformation_energy > post_transformation_energy;

  if (move_info->is_good) {
    move_info->delta = pre_transformation_energy - post_transformation_energy;
    move_info->pre_energy = pre_transformation_energy;
  }

  return move_info;
}

void TetMeshOptimizer::move_vertex(vertex v, pointType *coord_to_move) {
  std::vector<tetrahedra> incident_tetrahedras;
  VTfull(v, incident_tetrahedras);

  vertices[v] = coord_to_move;

  for (tetrahedra t : incident_tetrahedras)
    tetrahedras_energy[t] = getTetEnergy(t);
}

/// FIFTH PASS

void TetMeshOptimizer::fifth_pass() {
  std::cout << "\nStarting FIFTH pass" << std::endl;

  std::vector<std::unique_ptr<Split_tetrahedra_info>> tetrahedras_to_split;

  for (tetrahedra t = 0; t < mesh_->numTets(); t++) {
    if (tetrahedras_energy[t] > 100) {
      std::unique_ptr<Split_tetrahedra_info> split_info =
          get_split_tetrahedra_info(t);
      if (split_info->is_good) {
        tetrahedras_to_split.push_back(std::move(split_info));
      }
    }
  }
  std::cout << "Determined " << tetrahedras_to_split.size()
            << " tetrahedras to split" << std::endl;

  size_t num_tetrahedra_splitted = 0;

  for (auto &i : tetrahedras_to_split) {
    if (try_to_split_tetrahedra(std::move(i))) {

      std::cout << "Already splitted " << ++num_tetrahedra_splitted << " tets"
                << std::endl;
    }
  }
  std::cout << "Actually splitted " << num_tetrahedra_splitted << " tetrahedras"
            << std::endl;
  std::cout << "Finished FIFTH pass" << std::endl;
}

std::unique_ptr<TetMeshOptimizer::Split_tetrahedra_info>
TetMeshOptimizer::get_split_tetrahedra_info(tetrahedra t) {
  std::unique_ptr<Split_tetrahedra_info> split_info =
      std::make_unique<Split_tetrahedra_info>();
  if (mesh_->has_infinite_vertex(t)) {
    split_info->is_good = false;
    return split_info;
  }
  edge edge_1, edge_2;
  int best_edge = 0;
  double pre_transformation_energy = 0., post_transformation_energy = DBL_MAX,
         current_post_transformation_energy;
  for (int i = 1; i < 4; i++) {
    current_post_transformation_energy = 0.;
    edge_1 = std::make_pair(mesh_->get_i_th_vertex_of_tetrahedra(t, 0),
                            mesh_->get_i_th_vertex_of_tetrahedra(t, i));
    mesh_->oppositeTetEdgePair(t, edge_1, edge_2);
    std::unique_ptr<Split_edge_info> split_info_1 = get_split_edge_info(edge_1);
    std::unique_ptr<Split_edge_info> split_info_2 = get_split_edge_info(edge_2);
    pre_transformation_energy =
        std::max(split_info_1->pre_energy, split_info_2->pre_energy);
    for (const auto &[tet, energy] : split_info_1->energy_per_tet) {
      if (tet != t)
        current_post_transformation_energy =
            std::max(current_post_transformation_energy, energy);
    }
    current_post_transformation_energy =
        std::max(current_post_transformation_energy,
                 tetEnergyPositive(vertices[(split_info_1->edge).first],
                                   vertices[(split_info_2->edge).first],
                                   split_info_1->split_point,
                                   split_info_2->split_point));
    current_post_transformation_energy =
        std::max(current_post_transformation_energy,
                 tetEnergyPositive(vertices[(split_info_1->edge).first],
                                   vertices[(split_info_2->edge).second],
                                   split_info_1->split_point,
                                   split_info_2->split_point));
    current_post_transformation_energy =
        std::max(current_post_transformation_energy,
                 tetEnergyPositive(vertices[(split_info_1->edge).second],
                                   vertices[(split_info_2->edge).first],
                                   split_info_1->split_point,
                                   split_info_2->split_point));
    current_post_transformation_energy =
        std::max(current_post_transformation_energy,
                 tetEnergyPositive(vertices[(split_info_1->edge).second],
                                   vertices[(split_info_2->edge).second],
                                   split_info_1->split_point,
                                   split_info_2->split_point));
    if (current_post_transformation_energy < post_transformation_energy) {
      best_edge = i;
      post_transformation_energy = current_post_transformation_energy;
      split_info->split_point_1 = split_info_1->split_point;
      split_info->split_point_2 = split_info_2->split_point;
    }
  }
  split_info->is_good = post_transformation_energy < pre_transformation_energy;
  if (split_info->is_good) {
    split_info->pre_energy = pre_transformation_energy;
    split_info->delta = pre_transformation_energy - post_transformation_energy;
    split_info->edge_1 =
        std::make_pair(mesh_->get_i_th_vertex_of_tetrahedra(t, 0),
                       mesh_->get_i_th_vertex_of_tetrahedra(t, best_edge));
    mesh_->oppositeTetEdgePair(t, split_info->edge_1, split_info->edge_2);
    split_info->tetrahedra = t;
  }

  return split_info;
}

bool TetMeshOptimizer::try_to_split_tetrahedra(
    std::unique_ptr<Split_tetrahedra_info> info) {
  if (get_split_tetrahedra_info(info->tetrahedra)->is_good) {
    mesh_->pushVertex(info->split_point_1);
    // split_edge(info->edge_1, mesh_->numVertices() - 1);
    // mesh_->pushVertex(info->split_point_2);
    // split_edge(info->edge_2, mesh_->numVertices() - 1);
    return true;
  } else
    return false;
}

void TetMeshOptimizer::getMeshEdges(
    std::vector<std::pair<uint32_t, uint32_t>> &all_edges) const {
  for (size_t t = 0; t < mesh_->numTets(); t++) {
    const uint32_t *tn = mesh_->tet_node.data() + (t << 2);
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
void TetMeshOptimizer::get_edges_from_tetrahedras(
    std::vector<std::pair<uint32_t, uint32_t>> &all_edges,
    std::vector<tetrahedra> &tets) const {
  for (tetrahedra t : tets) {
    const uint32_t *tn = mesh_->tet_node.data() + (t << 2);
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

double TetMeshOptimizer::getTetEnergy(uint64_t t) const {
  const uint32_t *n1 = mesh_->tet_node.data() + (t << 2);
  if (n1[0] == INFINITE_VERTEX || n1[1] == INFINITE_VERTEX ||
      n1[2] == INFINITE_VERTEX || n1[3] == INFINITE_VERTEX)
    return DBL_MAX;
  return tetEnergy(vertices[n1[0]], vertices[n1[1]], vertices[n1[2]],
                   vertices[n1[3]]);
}

double TetMeshOptimizer::getTotalEnergy() {
  const uint32_t num_tets = mesh_->numTets();
  double total_energy = 0., current_tet_energy;
  for (uint32_t t = 0; t < num_tets; t++) {
    current_tet_energy = getTetEnergy(t);
    if (current_tet_energy != DBL_MAX)
      total_energy += current_tet_energy;
  }
  return total_energy;
}
double TetMeshOptimizer::getMaxEnergy() {
  const uint32_t num_tets = mesh_->numTets();
  double max_energy = 0., current_tet_energy;
  for (uint32_t t = 0; t < num_tets; t++) {
    if (mesh_->mark_tetrahedra[t] == DT_IN) {
      current_tet_energy = tetrahedras_energy[t];
      if (current_tet_energy != DBL_MAX)
        max_energy = std::max(max_energy, current_tet_energy);
    }
  }
  return max_energy;
}
double TetMeshOptimizer::getMeanEnergy() {
  const uint32_t num_tets = mesh_->numTets();
  double mean_energy = 0., current_tet_energy;
  uint32_t num_real_tets = 0;
  for (uint32_t t = 0; t < num_tets; t++) {
    if (mesh_->mark_tetrahedra[t] == DT_IN) {
      current_tet_energy = getTetEnergy(t);
      if (current_tet_energy != DBL_MAX) {
        num_real_tets++;
        mean_energy += current_tet_energy;
      }
    }
  }
  return mean_energy / num_real_tets;
}
