#include "optimize.h"
#include "tetmesh.h"
#include <cfloat>
#include <cstdint>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

/// HELPERS

template <class T> void remove_duplicates(std::vector<T> &vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

std::size_t edge_to_size_t(edge e) {
  if (e.first < e.second)
    return (static_cast<std::size_t>(e.first) << 32) | e.second;
  else
    return (static_cast<std::size_t>(e.second) << 32) | e.first;
}

edge make_edge(vertex v1, vertex v2) {
  if (v1 < v2)
    return std::make_pair(v1, v2);
  else
    return std::make_pair(v2, v1);
}

edge size_t_to_edge(std::size_t value) {
  unsigned int first = static_cast<unsigned int>(value >> 32);
  unsigned int second = static_cast<unsigned int>(value & 0xFFFFFFFF);
  return std::make_pair(first, second);
}
void TetMeshOptimizer::get_all_tets_energy() {
  const uint32_t num_tets = mesh_->numTets();
  tetrahedras_energy.clear();
  tetrahedras_energy.reserve(num_tets);
  for (tetrahedra t = 0; t < num_tets; t++) {
    tetrahedras_energy.push_back(get_quality_measure(t));
  }
}

bool is_oriented_good(const pointType *p1, const pointType *p2,
                      const pointType *p3, const pointType *p4) {
  return -pointType::orient3D(*p1, *p2, *p3, *p4) >= 0;
}

bool is_oriented_good(std::vector<pointType *> &points) {
  pointType *p1 = points[0];
  pointType *p2 = points[1];
  pointType *p3 = points[2];
  pointType *p4 = points[3];
  return -pointType::orient3D(*p1, *p2, *p3, *p4) >= 0;
}
bool is_oriented_good(const std::vector<pointType *> &points) {
  pointType *p1 = points[0];
  pointType *p2 = points[1];
  pointType *p3 = points[2];
  pointType *p4 = points[3];
  return -pointType::orient3D(*p1, *p2, *p3, *p4) >= 0;
}

double TetMeshOptimizer::get_quality_measure(const pointType *p1,
                                             const pointType *p2,
                                             const pointType *p3,
                                             const pointType *p4) {
  if (!is_oriented_good(p1, p2, p3, p4))
    return DBL_MAX;
  else
    return quality_measure_(p1, p2, p3, p4);
}

double TetMeshOptimizer::get_quality_measure(const tetrahedra t) {
  if (mesh_->has_infinite_vertex(t))
    return DBL_MAX;
  std::vector<pointType *> points = mesh_->getTetPoints(t);
  if (!is_oriented_good(points))
    return DBL_MAX;
  else
    return quality_measure_(
        mesh_->vertices[mesh_->get_i_th_vertex_of_tetrahedra(t, 0)],
        mesh_->vertices[mesh_->get_i_th_vertex_of_tetrahedra(t, 1)],
        mesh_->vertices[mesh_->get_i_th_vertex_of_tetrahedra(t, 2)],
        mesh_->vertices[mesh_->get_i_th_vertex_of_tetrahedra(t, 3)]);
}

double TetMeshOptimizer::get_quality_measure(
    const std::vector<pointType *> &tet_points) {
  if (!is_oriented_good(tet_points))
    return DBL_MAX;
  else
    return quality_measure_(tet_points[0], tet_points[1], tet_points[2],
                            tet_points[3]);
}

void TetMeshOptimizer::log_message(std::string message) {
  auto now = std::chrono::system_clock::now();
  std::time_t now_c = std::chrono::system_clock::to_time_t(now);
  std::tm tm = *std::localtime(&now_c);

  std::ostringstream oss;
  oss << std::put_time(&tm, "%F %T"); // équivalent ISO: YYYY-MM-DD HH:MM:SS
  std::string formatted_time = oss.str();

  *time_log << formatted_time << ": " << message << std::endl;
}
void TetMeshOptimizer::log_energy() {
  *mean_energy_log << getMeanEnergy() << std::endl;
  *max_energy_log << getMaxEnergy() << std::endl;
}

bool TetMeshOptimizer::optimize() {
  if (verbose) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Starting optimizing pass\n" << std::endl;
  }
  if (log)
    log_message("Start optimization pass");

  mesh_->checkMesh(false);

  std::unique_ptr<UpdatableQueue<Split_edge_info, double, size_t>> split_queue;
  if (random)
    split_queue =
        std::make_unique<UpdatableRandomQueue<Split_edge_info, size_t>>();
  else
    split_queue = std::make_unique<
        UpdatablePriorityQueue<Split_edge_info, double, size_t>>();

  if (log)
    log_message("Start first pass");

  bool splitted = first_pass(optim_field_edge_split, *split_queue);

  if (log) {
    log_energy();
    log_message("End first pass");
  }

  mesh_->checkMesh(false);

  std::unique_ptr<UpdatableQueue<Collapse_info, double, size_t>> collapse_queue;
  if (random)
    collapse_queue =
        std::make_unique<UpdatableRandomQueue<Collapse_info, size_t>>();
  else
    collapse_queue = std::make_unique<
        UpdatablePriorityQueue<Collapse_info, double, size_t>>();

  if (log)
    log_message("Start second pass");

  bool collapsed = second_pass(optim_field_collapse, *collapse_queue);

  for (const auto &[nom, temps] : tempsCumules) {
    std::cout << nom << " a pris " << temps << " µs au total.\n";
  }
  tempsCumules.clear();

  if (log) {
    log_energy();
    log_message("End second pass");
  }

  mesh_->checkMesh(false);

  bool swapped = third_pass(optim_field_edge_swap, optim_field_face_swap);

  for (const auto &[nom, temps] : tempsCumules) {
    std::cout << nom << " a pris " << temps << " µs au total.\n";
  }
  tempsCumules.clear();

  mesh_->checkMesh(false);

  std::unique_ptr<UpdatableQueue<Move_info, double, vertex>> move_queue;
  if (random)
    move_queue = std::make_unique<UpdatableRandomQueue<Move_info, vertex>>();
  else
    move_queue =
        std::make_unique<UpdatablePriorityQueue<Move_info, double, vertex>>();

  if (log)
    log_message("Start fourth pass");

  bool moved = fourth_pass(optim_field_vertex_move, *move_queue);

  if (log) {
    log_energy();
    log_message("End fourth pass");
  }

  mesh_->checkMesh(false);

  std::unique_ptr<UpdatableQueue<Split_tetrahedra_info, double, tetrahedra>>
      split_tet_queue;
  if (random)
    split_tet_queue = std::make_unique<
        UpdatableRandomQueue<Split_tetrahedra_info, tetrahedra>>();
  else
    split_tet_queue = std::make_unique<
        UpdatablePriorityQueue<Split_tetrahedra_info, double, tetrahedra>>();

  if (log)
    log_message("Start fifth pass");

  bool tet_splitted = fifth_pass(optim_field_tet_split, *split_tet_queue);

  for (const auto &[nom, temps] : tempsCumules) {
    std::cout << nom << " a pris " << temps << " µs au total.\n";
  }
  tempsCumules.clear();

  if (log)
    log_message("End fifth pass");

  mesh_->checkMesh(false);

  if (verbose) {
    std::cout << "\nFinished optimizing pass" << std::endl;
    std::cout << "-----------------------------------------------\n"
              << std::endl;
  }

  if (log)
    log_message("End optimization pass");

  return splitted || collapsed || swapped || moved || tet_splitted;
}

/// FIRST PASS
bool TetMeshOptimizer::first_pass(
    double Split_edge_info::*field_to_optimize,
    UpdatableQueue<Split_edge_info, double, size_t> &queue) {

  if (verbose)
    std::cout << "Starting FIRST pass" << std::endl;

  std::vector<edge> edges;
  mesh_->getMeshEdges(edges);

  for (edge e : edges) {
    std::unique_ptr<Split_edge_info> split_info = get_split_edge_info(e);
    if (split_info->is_good) {
      queue.push(std::move(split_info), (*split_info).*field_to_optimize);
    }
  }

  if (verbose)
    std::cout << "Determined " << queue.size() << " edges to split"
              << std::endl;

  int num_splitted_edges = 0;

  while (!queue.empty()) {
    auto split_info = queue.pop();
    if (split_info.has_value()) {
      split_edge_and_update(std::move(split_info.value()), queue,
                            field_to_optimize);
      num_splitted_edges++;
    } else
      break;
  }

  if (verbose)
    std::cout << "Actually splitted " << num_splitted_edges << " edges"
              << std::endl;

  if (verbose)
    std::cout << "Finished FIRST pass" << std::endl;
  return num_splitted_edges > 0;
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

  double energy_tet_1 = get_quality_measure(tetrahedra_points_v1);
  double energy_tet_2 = get_quality_measure(tetrahedra_points_v2);

  // Prevent inversion
  if (!is_oriented_good(tetrahedra_points_v1) ||
      !is_oriented_good(tetrahedra_points_v2))
    return DBL_MAX;

  return std::max(energy_tet_1, energy_tet_2);
}

pointType *TetMeshOptimizer::get_split_point(edge e, double t) {
  return (((vector3d(mesh_->vertices[e.first])) * t) +
          ((vector3d(mesh_->vertices[e.second])) * (1 - t)))
      .toExplicitPoint();
}

std::unique_ptr<TetMeshOptimizer::Split_edge_info>
TetMeshOptimizer::get_split_edge_info(edge e) {

  std::unique_ptr<Split_edge_info> split_info =
      std::make_unique<Split_edge_info>();
  split_info->edge = e;
  split_info->id = edge_to_size_t(e);
  split_info->is_good = false;

  std::vector<tetrahedra> incident_tetrahedras;
  mesh_->ETfull(e.first, e.second, incident_tetrahedras);

  double pre_transformation_energy = 0., post_transformation_energy = DBL_MAX,
         energy_of_split;

  // implicitPoint_LNC *potential_split_point =
  //     new implicitPoint_LNC(mesh_->vertices[e.first]->toExplicit3D(),
  //                           mesh_->vertices[e.second]->toExplicit3D(), 0.5);
  //  TODO : Implicit point not working great

  int LOD = 2;

  for (int i = 1; i < LOD; i++) {
    pointType *split_point = get_split_point(e, (double)i / LOD);
    double current_pre_transformation_energy = 0.,
           current_post_transformation_energy = 0.;

    for (tetrahedra tet : incident_tetrahedras) {
      if (!mesh_->has_infinite_vertex(tet) &&
          (mesh_->mark_tetrahedra[tet] == DT_IN || !optimize_only_DT_IN)) {
        energy_of_split = get_energy_from_splitting(tet, e, split_point);
        if (energy_of_split == DBL_MAX) {
          split_info->is_good = false;
          return split_info;
        }
        current_pre_transformation_energy = std::max(
            current_pre_transformation_energy, tetrahedras_energy[tet]);
        current_post_transformation_energy =
            std::max(current_post_transformation_energy, energy_of_split);
      }
    }

    if (current_pre_transformation_energy != DBL_MAX &&
        current_post_transformation_energy <
            current_pre_transformation_energy) {
      split_info->is_good = true;
      if (current_post_transformation_energy < post_transformation_energy) {
        split_info->split_point = split_point;
        split_info->delta = current_pre_transformation_energy -
                            current_post_transformation_energy;
        split_info->pre_energy = current_pre_transformation_energy;
        if (mesh_->isOnBoundary(e.first, e.second))
          split_info->prioritize = 1;
      }
    }
  }

  return split_info;
}

void TetMeshOptimizer::split_edge_and_update(
    std::unique_ptr<Split_edge_info> split_info,
    UpdatableQueue<Split_edge_info, double, size_t> &queue,
    double Split_edge_info::*field_to_optimize) {

  auto split_result = std::make_unique<Split_edge_result>();
  split_edge(std::move(split_info), split_result.get());

  for (edge edge_to_update : split_result->impacted_edges) {
    std::unique_ptr<Split_edge_info> split_info =
        get_split_edge_info(edge_to_update);
    if (split_info->is_good)
      queue.set(std::move(split_info), (*split_info).*field_to_optimize);
    else
      queue.remove(split_info->id);
  }
}

void TetMeshOptimizer::split_edge(std::unique_ptr<Split_edge_info> split_info,
                                  Split_edge_result *split_result) {

  split_result->success = false;

  edge edge_to_split = split_info->edge;
  vertex split_vertex = split_info->split_vertex;
  pointType *split_point;

  if (split_vertex == INFINITE_VERTEX) {
    if (!split_info->split_point) {
      split_point = get_split_point(edge_to_split);
    } else {
      split_point = split_info->split_point;
    }
    mesh_->pushVertex(split_point);
    split_vertex = mesh_->numVertices() - 1;
  }

  std::unordered_set<vertex> seen_vertices;
  split_result->impacted_edges.push_back(
      make_edge(edge_to_split.first, split_vertex));
  split_result->impacted_edges.push_back(
      make_edge(edge_to_split.second, split_vertex));

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
    split_result->impacted_tetrahedras.push_back(tet);
    mesh_->mark_tetrahedra.push_back(mesh_->mark_tetrahedra[tet]);
    split_result->impacted_tetrahedras.push_back(
        mesh_->get_tetrahedra_index_from_corner(current_corner));
    if (!mesh_->has_infinite_vertex(tet))
      mesh_->inc_tet[split_vertex] = tet;

    opposite_of_v1 =
        mesh_->tet_neigh[mesh_->get_corner_in_tet(tet, edge_to_split.first)];

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
          // Change mesh_->inc_tet for v2 because v2 does not belong anymore to
          // the same tets
          if (!mesh_->has_infinite_vertex(tet))
            mesh_->inc_tet[edge_to_split.second] =
                mesh_->get_tetrahedra_index_from_corner(current_corner);
          mesh_->tet_neigh.push_back(
              mesh_->get_corner_in_tet(tet, edge_to_split.first));
          mesh_->tet_neigh[mesh_->get_corner_in_tet(tet, edge_to_split.first)] =
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
              mesh_->get_index_of_corner_in_tet(opposite_former_corner);
          corner corner_to_add = start_corner + (index_of_processing_tet << 2) +
                                 index_opposite_former_corner;
          mesh_->tet_neigh.push_back(corner_to_add);

          if (!seen_vertices.contains(current_vertex) &&
              current_vertex != INFINITE_VERTEX) {
            for (vertex v :
                 {edge_to_split.first, edge_to_split.second, split_vertex})
              split_result->impacted_edges.push_back(
                  make_edge(current_vertex, v));
            seen_vertices.insert(current_vertex);
          }
        }
      }
      current_corner++;
    }
    // Finally we replace v2 by split_vertex in the old tetrahedra, aside from
    // changing the opposite corner for v2 which we did in previous loop there
    // is nothing to change
    mesh_->tet_node[mesh_->get_corner_in_tet(tet, edge_to_split.second)] =
        split_vertex;

    // Updating energy
    tetrahedras_energy.push_back(get_quality_measure(
        mesh_->get_tetrahedra_index_from_corner(current_corner - 4)));
    tetrahedras_energy[tet] = get_quality_measure(tet);
  }

  if (mesh_->inc_tet[split_vertex] == UINT64_MAX)
    throw std::runtime_error("Not managed to init inc_tet for split_vertex");

  remove_duplicates(split_result->impacted_edges);
  remove_duplicates(split_result->impacted_tetrahedras);

  split_result->success = true;
}

/// SECOND PASS

bool TetMeshOptimizer::second_pass(
    double Collapse_info::*field_to_optimize,
    UpdatableQueue<Collapse_info, double, size_t> &queue) {

  if (verbose)
    std::cout << "\nStarting SECOND pass" << std::endl;

  std::vector<edge> edges;
  mesh_->getMeshEdges(edges);

  for (edge e : edges) {
    std::unique_ptr<Collapse_info> collapse_info = get_collapse_info(e);
    if (collapse_info->is_good) {
      queue.push(std::move(collapse_info), (*collapse_info).*field_to_optimize);
    }
  }

  size_t num_edge_collapsed = 0;

  if (verbose)
    std::cout << "Determined " << queue.size() << " edges to collapse"
              << std::endl;

  while (!queue.empty()) {
    auto collapse_info = queue.pop();
    if (collapse_info.has_value()) {
      collapse_and_update(std::move(collapse_info.value()), queue,
                          field_to_optimize);
      num_edge_collapsed++;
    } else
      break;
  }

  if (verbose)
    std::cout << "Actually collapsed " << num_edge_collapsed << " edges"
              << std::endl;

  if (verbose)
    std::cout << "Finished SECOND pass" << std::endl;

  temp_remap_vertex.clear();
  return num_edge_collapsed > 0;
}

void TetMeshOptimizer::collapse_and_update(
    std::unique_ptr<Collapse_info> collapse_info,
    UpdatableQueue<Collapse_info, double, size_t> &queue,
    double Collapse_info::*field_to_optimize) {

  // ChronoAuto chrono(__func__);

  auto collapse_info_check = get_collapse_info(collapse_info->edge, false);

  if (!collapse_info_check->is_good) {
    throw std::runtime_error("Check for collapsing is not good");
  }

  auto collapse_result = std::make_unique<Collapse_result>();
  collapse(std::move(collapse_info), collapse_result.get(), false);

  // Remove edges that have been suppressed
  for (edge e : collapse_result->removed_edges) {
    queue.remove(edge_to_size_t(e));
  }

  for (edge e : collapse_result->impacted_edges) {
    std::unique_ptr<Collapse_info> collapse_info = get_collapse_info(e);
    if (collapse_info->is_good) {
      queue.set(std::move(collapse_info), (*collapse_info).*field_to_optimize);
    } else {
      queue.remove(edge_to_size_t(e));
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

  mesh_->VV(e.first, one_ring_v1);
  mesh_->VV(e.second, one_ring_v2);
  mesh_->OneRing(e, one_ring_e, incident_tetrahedras_e);

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
    edge e, std::vector<tetrahedra> &incident_tetrahedras_v2, bool debug) {

  // ChronoAuto chrono(__func__);
  if (incident_tetrahedras_v2.size() == 0) {
    return DBL_MAX;
  }

  auto [v1, v2] = e;
  double post_transformation_energy = 0.;

  int num_tet_taken_into_account = 0;

  for (tetrahedra t2 : incident_tetrahedras_v2) {
    if (!mesh_->has_infinite_vertex(t2) &&
        (mesh_->mark_tetrahedra[t2] == DT_IN || !optimize_only_DT_IN) &&
        !mesh_->tetHasVertex(t2, v1)) {

      std::vector<pointType *> t2_points = mesh_->getTetPoints(t2);
      t2_points[mesh_->get_index_of_vertex_in_tet(v2, t2)] =
          mesh_->vertices[v1];

      if (!is_oriented_good(t2_points)) {
        return DBL_MAX;
      }

      double energy_prevision = get_quality_measure(t2_points);

      post_transformation_energy =
          std::max(post_transformation_energy, get_quality_measure(t2_points));

      num_tet_taken_into_account++;
    }
  }

  if (num_tet_taken_into_account == 0)
    post_transformation_energy = DBL_MAX;

  return post_transformation_energy;
}

std::unique_ptr<TetMeshOptimizer::Collapse_info>
TetMeshOptimizer::get_collapse_info(edge &e, bool debug) {
  // ChronoAuto chrono(__func__);

  auto collapse_info = std::make_unique<Collapse_info>();

  collapse_info->id = edge_to_size_t(e);
  collapse_info->edge = e;
  // Consider e := v1 -- v2
  auto [v1, v2] = e;

  if (!link_condition(e) || mesh_->isOnBoundary(v1, v2)) {
    collapse_info->is_good = false;
    return collapse_info;
  }

  std::vector<tetrahedra> incident_tetrahedras_v1, incident_tetrahedras_v2;
  mesh_->VTfull(v1, incident_tetrahedras_v1);
  mesh_->VTfull(v2, incident_tetrahedras_v2);

  double pre_transformation_energy_v1 = 0., pre_transformation_energy_v2 = 0.,
         post_transformation_energy_v1, post_transformation_energy_v2,
         post_transformation_energy;

  // Compute pre-tranformation energy
  for (tetrahedra t : incident_tetrahedras_v2) {
    if (!mesh_->has_infinite_vertex(t) &&
        (mesh_->mark_tetrahedra[t] == DT_IN || !optimize_only_DT_IN))
      pre_transformation_energy_v1 =
          std::max(pre_transformation_energy_v1, tetrahedras_energy[t]);
  }

  for (tetrahedra t : incident_tetrahedras_v1) {
    if (!mesh_->has_infinite_vertex(t) &&
        (mesh_->mark_tetrahedra[t] == DT_IN || !optimize_only_DT_IN))
      pre_transformation_energy_v2 =
          std::max(pre_transformation_energy_v2, tetrahedras_energy[t]);
  }

  post_transformation_energy_v1 =
      get_energy_from_collapsing(e, incident_tetrahedras_v2, debug);
  post_transformation_energy_v2 = get_energy_from_collapsing(
      std::make_pair(v2, v1), incident_tetrahedras_v1, debug);

  collapse_info->is_good =
      (post_transformation_energy_v1 < pre_transformation_energy_v1 &&
       pre_transformation_energy_v1 != DBL_MAX &&
       post_transformation_energy_v1 != DBL_MAX) ||
      (post_transformation_energy_v2 < pre_transformation_energy_v2 &&
       pre_transformation_energy_v2 != DBL_MAX &&
       post_transformation_energy_v2 != DBL_MAX);

  if (collapse_info->is_good) {
    double delta_v1 =
        pre_transformation_energy_v1 - post_transformation_energy_v1;
    double delta_v2 =
        pre_transformation_energy_v2 - post_transformation_energy_v2;
    if (delta_v1 > delta_v2) {
      collapse_info->delta = delta_v1;
      collapse_info->pre_energy = pre_transformation_energy_v1;
    } else {
      collapse_info->delta = delta_v2;
      collapse_info->pre_energy = pre_transformation_energy_v2;
    }

    if (delta_v1 > delta_v2 && !mesh_->isOnBoundary(v2)) {
      collapse_info->collapse_vertex = v1;
    } else if (delta_v1 <= delta_v2 && !mesh_->isOnBoundary(v1)) {
      collapse_info->collapse_vertex = v2;
    } else {
      collapse_info->is_good = false;
    }
  }
  return collapse_info;
}

void TetMeshOptimizer::collapse(std::unique_ptr<Collapse_info> collapse_info,
                                Collapse_result *collapse_result, bool debug) {

  // ChronoAuto chrono(__func__);

  auto [v1, v2] = collapse_info->edge;

  if (v1 != collapse_info->collapse_vertex) {
    v2 = v1;
    v1 = collapse_info->collapse_vertex;
  }

  collapse_result->success = false;

  if (v1 == INFINITE_VERTEX || v2 == INFINITE_VERTEX ||
      mesh_->isOnBoundary(v2) || !mesh_->hasEdge(v1, v2) ||
      !link_condition(collapse_info->edge)) {
    return;
  }

  std::unordered_set<vertex> seen_vertices;

  std::vector<tetrahedra> incident_tetrahedras, incident_tetrahedras_v2,
      deformed_tetrahedras;
  std::vector<vertex> neighbour_v2;
  mesh_->ETfull(v1, v2, incident_tetrahedras);
  mesh_->VTfull(v2, incident_tetrahedras_v2);
  mesh_->VV(v2, neighbour_v2);

  for (vertex neighbour : neighbour_v2) {
    collapse_result->removed_edges.push_back(make_edge(v2, neighbour));
  }

  edge opposite_edge;
  corner opposite_corner_1, opposite_corner_2, potential_corner_to_change;
  tetrahedra t1, t2, next_to_delete;

  std::vector<tetrahedra> tet_to_delete;

  if (incident_tetrahedras.empty()) {
    return;
  }

  mesh_->inc_tet[v1] = UINT64_MAX;
  mesh_->inc_tet[v2] = UINT64_MAX;

  for (tetrahedra t : incident_tetrahedras) {

    mesh_->oppositeTetEdgePair(t, make_edge(v1, v2), opposite_edge);

    opposite_corner_1 = mesh_->tet_neigh[mesh_->get_corner_in_tet(t, v2)];
    opposite_corner_2 = mesh_->tet_neigh[mesh_->get_corner_in_tet(t, v1)];

    t1 = mesh_->get_tetrahedra_index_from_corner(opposite_corner_1);
    t2 = mesh_->get_tetrahedra_index_from_corner(opposite_corner_2);

    mesh_->setMutualNeighbors(opposite_corner_1, opposite_corner_2);

    if (!mesh_->isGhost(t1)) {
      mesh_->inc_tet[v1] = t1;
    }

    if (opposite_edge.first != INFINITE_VERTEX && !mesh_->isGhost(t1)) {
      mesh_->inc_tet[opposite_edge.first] = t1;
      mesh_->inc_tet[opposite_edge.second] = t1;
    }
    if (opposite_edge.second != INFINITE_VERTEX && !mesh_->isGhost(t2)) {
      mesh_->inc_tet[opposite_edge.first] = t2;
      mesh_->inc_tet[opposite_edge.second] = t2;
    }
    tet_to_delete.push_back(t);
  }

  for (tetrahedra t : incident_tetrahedras_v2) {
    collapse_result->impacted_tetrahedras.push_back(t);
    potential_corner_to_change = mesh_->get_corner_in_tet(t, v2);
    if (potential_corner_to_change != UINT64_MAX &&
        !mesh_->tetHasVertex(t, v1)) {
      if (!mesh_->isGhost(t)) {
        mesh_->inc_tet[v1] = t;
      }
      mesh_->tet_node[potential_corner_to_change] = v1;
      // Update energy for t2
      tetrahedras_energy[t] = get_quality_measure(t);
      deformed_tetrahedras.push_back(t);
    }
  }

  if (mesh_->inc_tet[v1] == UINT64_MAX) {
    throw std::runtime_error("Inc_tet for " + std::to_string(v1) +
                             " cannot be properly set after collapsing edge " +
                             std::to_string(v1) + ", " + std::to_string(v2));
  }

  // Delete in descendant order (it could be possible to swap t with a tet at
  // the end that is also to delete)

  std::sort(tet_to_delete.begin(), tet_to_delete.end());

  while (!tet_to_delete.empty()) {
    next_to_delete = tet_to_delete.back();
    tet_to_delete.pop_back();
    collapse_result->removed_tetrahedras.push_back(next_to_delete);
    remove_tetrahedra(next_to_delete);
  }

  // Add the edges that have been impacted by removing v2 which implies to remap
  // the last vertex to v2, in this case all the edges incident to last_vertex
  // are now incident to v2
  std::vector<vertex> impacted_vertex_remap;
  vertex last_vertex = mesh_->numVertices() - 1;
  if (v2 != last_vertex) {
    mesh_->VV(last_vertex, impacted_vertex_remap);
    mesh_->VTfull(last_vertex, collapse_result->impacted_tetrahedras);
    for (vertex neighbour : impacted_vertex_remap) {
      collapse_result->removed_edges.push_back(
          make_edge(last_vertex, neighbour));
    }
  }

  remove_vertex(v2);

  remap_vertex(v1);

  std::vector<vertex> neighbours;

  for (tetrahedra t : deformed_tetrahedras) {
    remap_tetrahedra(t);
    collapse_result->impacted_tetrahedras.push_back(t);
    for (int i = 0; i < 4; i++) {
      vertex current_vertex = mesh_->get_i_th_vertex_of_tetrahedra(t, i);
      if (!seen_vertices.contains(current_vertex) &&
          current_vertex != INFINITE_VERTEX) {
        seen_vertices.insert(current_vertex);
        neighbours.clear();
        mesh_->VV(current_vertex, neighbours);
        for (vertex neighbour : neighbours) {
          collapse_result->impacted_edges.push_back(
              make_edge(current_vertex, neighbour));
          collapse_result->set_impacted_edges.insert(
              make_edge(current_vertex, neighbour));
        }
      }
    }
  }

  for (vertex current_vertex : neighbour_v2) {
    if (!seen_vertices.contains(current_vertex) &&
        current_vertex != INFINITE_VERTEX &&
        current_vertex < mesh_->numVertices()) {
      seen_vertices.insert(current_vertex);
      neighbours.clear();
      mesh_->VV(current_vertex, neighbours);
      for (vertex neighbour : neighbours) {
        collapse_result->impacted_edges.push_back(
            make_edge(current_vertex, neighbour));
      }
    }
  }

  mesh_->VTfull(v1, collapse_result->impacted_tetrahedras);
  temp_remap_vertex.clear();

  remove_duplicates(collapse_result->impacted_edges);
  remove_duplicates(collapse_result->impacted_tetrahedras);
  remove_duplicates(collapse_result->removed_edges);
  remove_duplicates(collapse_result->removed_tetrahedras);

  collapse_result->success = true;
}

/// THIRD PASS

bool TetMeshOptimizer::third_pass(
    double Swap_edge_info::*field_to_optimize_edge,
    double Swap_face_info::*field_to_optimize_face) {

  // To iterate over all faces once we can represent uniquely a face with the
  // two opposite corner they induce this way we iterate over all corners then
  // we are considering the corner only if it's smaller than its
  // mesh_->tet_neigh then we add those pair of corner in a queue for swapping
  // if it's worth it, then when we are popping from the queue we just verify if
  // both corner are still opposed (i.e. they haven't been touched by an other
  // swap)

  if (verbose)
    std::cout << "\nStarting THIRD pass" << std::endl;

  std::unique_ptr<UpdatableQueue<Swap_edge_info, double, size_t>> edge_queue;
  if (random)
    edge_queue =
        std::make_unique<UpdatableRandomQueue<Swap_edge_info, size_t>>();
  else
    edge_queue = std::make_unique<
        UpdatablePriorityQueue<Swap_edge_info, double, size_t>>();

  if (log)
    log_message("Start third pass edge");

  bool swapped_edge = third_pass_edge(field_to_optimize_edge, *edge_queue);

  if (log) {
    log_energy();
    log_message("End third pass edge");
  }

  std::unique_ptr<UpdatableQueue<Swap_face_info, double, corner>> face_queue;
  if (random)
    face_queue =
        std::make_unique<UpdatableRandomQueue<Swap_face_info, corner>>();
  else
    face_queue = std::make_unique<
        UpdatablePriorityQueue<Swap_face_info, double, corner>>();

  if (log)
    log_message("Start third pass face");

  bool swapped_face = third_pass_face(field_to_optimize_face, *face_queue);

  if (log) {
    log_energy();
    log_message("End third pass face");
  }

  if (verbose)
    std::cout << "Finished THIRD pass" << std::endl;

  return swapped_edge || swapped_face;
}

bool TetMeshOptimizer::third_pass_face(
    double Swap_face_info::*field_to_optimize,
    UpdatableQueue<Swap_face_info, double, corner> &queue) {

  for (corner face = 0; face < mesh_->tet_node.size(); face++) {
    std::unique_ptr<Swap_face_info> swap_info = get_swap_face_info(face);
    if (swap_info->is_good)
      queue.push(std::move(swap_info), (*swap_info).*field_to_optimize);
  }

  if (verbose)
    std::cout << "Determined " << queue.size() << " faces to swap" << std::endl;

  uint32_t num_face_swapped = 0;

  while (!queue.empty()) {
    auto swap_info = queue.pop();
    if (swap_info.has_value()) {
      swap_face_and_update(std::move(swap_info.value()), queue,
                           field_to_optimize);
      num_face_swapped++;
    } else
      break;
  }

  if (verbose)
    std::cout << "Actually swapped " << num_face_swapped << " faces"
              << std::endl;

  return num_face_swapped > 0;
}

bool TetMeshOptimizer::third_pass_edge(
    double Swap_edge_info::*field_to_optimize,
    UpdatableQueue<Swap_edge_info, double, size_t> &queue) {
  std::vector<edge> edges;
  mesh_->getMeshEdges(edges);

  for (edge e : edges) {
    std::unique_ptr<Swap_edge_info> swap_info = get_swap_edge_info(e);
    if (swap_info->is_good)
      queue.push(std::move(swap_info), (*swap_info).*field_to_optimize);
  }

  if (verbose)
    std::cout << "Determined " << queue.size() << " edges to swap" << std::endl;

  int num_edges_swapped = 0;

  while (!queue.empty()) {
    auto swap_info = queue.pop();
    if (swap_info.has_value()) {
      swap_edge_and_update(std::move(swap_info.value()), queue,
                           field_to_optimize);
      num_edges_swapped++;
    } else
      break;
  }
  if (verbose)
    std::cout << "Actually swapped " << num_edges_swapped << " edges"
              << std::endl;
  return num_edges_swapped > 0;
}

void TetMeshOptimizer::swap_face_and_update(
    std::unique_ptr<Swap_face_info> swap_info,
    UpdatableQueue<Swap_face_info, double, corner> &queue,
    double Swap_face_info::*field_to_optimize) {

  auto swap_result = std::make_unique<Swap_face_result>();
  swap_face(std::move(swap_info), swap_result.get());

  for (corner face : swap_result->impacted_faces) {
    std::unique_ptr<Swap_face_info> swap_info = get_swap_face_info(face);
    queue.remove(mesh_->tet_neigh[face]);
    if (swap_info->is_good) {
      queue.set(std::move(swap_info), (*swap_info).*field_to_optimize);
    } else {
      queue.remove(face);
    }
  }
}
void TetMeshOptimizer::swap_edge_and_update(
    std::unique_ptr<Swap_edge_info> swap_info,
    UpdatableQueue<Swap_edge_info, double, size_t> &queue,
    double Swap_edge_info::*field_to_optimize) {
  // ChronoAuto chrono(__func__);
  std::vector<edge> impacted_edges;
  std::vector<tetrahedra> impacted_tets, removed_tets;
  auto swap_result = std::make_unique<Swap_edge_result>();
  swap_edge(std::move(swap_info), swap_result.get());

  for (edge e : swap_result->impacted_edges) {
    std::unique_ptr<Swap_edge_info> swap_info = get_swap_edge_info(e);
    if (swap_info->is_good) {
      queue.set(std::move(swap_info), (*swap_info).*field_to_optimize);
    } else {
      queue.remove(edge_to_size_t(e));
    }
  }
}

void TetMeshOptimizer::swap_edge(std::unique_ptr<Swap_edge_info> swap_info,
                                 Swap_edge_result *swap_result) {
  swap_result->success = false;

  edge edge_to_swap = swap_info->edge;
  vertex collapse_vertex = swap_info->collapse_vertex;

  std::vector<vertex> one_ring;
  std::vector<tetrahedra> incident_tetrahedras;
  std::vector<vertex> neigh_collapse;
  edge opposite_edge;

  mesh_->OneRing(edge_to_swap, one_ring, incident_tetrahedras);

  // Add impacted edges that are connected to one of the two endpoint of the
  // edge to swap, then also add newly created edges
  for (vertex neigh : one_ring) {
    swap_result->impacted_edges.push_back(make_edge(edge_to_swap.first, neigh));
    swap_result->impacted_edges.push_back(
        make_edge(edge_to_swap.second, neigh));
    if (neigh != collapse_vertex) {
      swap_result->impacted_edges.push_back(make_edge(neigh, collapse_vertex));
    }
  }

  // Add impacted edges that are around the edge we are going to swap
  for (tetrahedra t : incident_tetrahedras) {
    mesh_->oppositeTetEdgePair(t, edge_to_swap, opposite_edge);
    swap_result->impacted_edges.push_back(opposite_edge);
  }

  // Do the swapping
  std::unique_ptr<Split_edge_result> split_result =
      std::make_unique<Split_edge_result>();
  split_edge(std::move(std::make_unique<Split_edge_info>(edge_to_swap)),
             split_result.get());

  std::unique_ptr<Collapse_result> collapse_result =
      std::make_unique<Collapse_result>();

  edge edge_to_collapse = make_edge(collapse_vertex, mesh_->numVertices() - 1);

  if (!link_condition(edge_to_collapse)) {
    edge edge_to_close =
        make_edge(edge_to_swap.second, mesh_->numVertices() - 1);

    std::unique_ptr<Collapse_info> collapse_info =
        std::make_unique<Collapse_info>(edge_to_size_t(edge_to_close),
                                        edge_to_close, edge_to_swap.second);

    collapse(std::move(collapse_info), collapse_result.get());
    swap_result->impacted_edges.clear();
    return;
  }

  std::unique_ptr<Collapse_info> collapse_info =
      std::make_unique<Collapse_info>(edge_to_size_t(edge_to_collapse),
                                      edge_to_collapse, collapse_vertex);

  collapse(std::move(collapse_info), collapse_result.get());

  swap_result->impacted_edges.reserve(
      swap_result->impacted_edges.size() +
      distance(collapse_result->impacted_edges.begin(),
               collapse_result->impacted_edges.end()));
  swap_result->impacted_edges.insert(swap_result->impacted_edges.end(),
                                     collapse_result->impacted_edges.begin(),
                                     collapse_result->impacted_edges.end());

  remove_duplicates(swap_result->impacted_edges);
  remove_duplicates(swap_result->removed_edges);

  swap_result->success = true;
}

std::unique_ptr<TetMeshOptimizer::Swap_face_info>
TetMeshOptimizer::get_swap_face_info(corner face) {
  std::unique_ptr<Swap_face_info> swap_info =
      std::make_unique<Swap_face_info>();
  swap_info->id = face;
  swap_info->face = face;
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
  }

  return swap_info;
}

double TetMeshOptimizer::get_energy_from_swapping_face(corner face) {

  const uint64_t index_in_tet = face & 3;
  const corner tet_basis = face - index_in_tet;

  const uint64_t r0 = tet_basis + mesh_->tetON1(index_in_tet),
                 r1 = tet_basis + mesh_->tetON3(index_in_tet),
                 r2 = tet_basis + mesh_->tetON2(index_in_tet);
  const uint32_t c0 = mesh_->tet_node[r0], c1 = mesh_->tet_node[r1],
                 c2 = mesh_->tet_node[r2], c3 = mesh_->tet_node[face];

  const uint64_t g00 = mesh_->tet_neigh[r0], g01 = mesh_->tet_neigh[r1],
                 g02 = mesh_->tet_neigh[r2];

  const uint64_t orx = mesh_->tet_neigh[face];
  const uint64_t opp = orx & (~3);
  const uint64_t or0 = mesh_->tetCornerAtVertex(opp, c0);
  const uint64_t or1 = mesh_->tetCornerAtVertex(opp, c1);
  const uint64_t or2 = mesh_->tetCornerAtVertex(opp, c2);

  const uint64_t g10 = mesh_->tet_neigh[or0], g11 = mesh_->tet_neigh[or1],
                 g12 = mesh_->tet_neigh[or2];

  const uint32_t oc = mesh_->tet_node[orx];

  return std::max(
      get_quality_measure(mesh_->vertices[c3], mesh_->vertices[oc],
                          mesh_->vertices[c1], mesh_->vertices[c2]),
      std::max(get_quality_measure(mesh_->vertices[c3], mesh_->vertices[c0],
                                   mesh_->vertices[oc], mesh_->vertices[c2]),
               get_quality_measure(mesh_->vertices[c3], mesh_->vertices[c0],
                                   mesh_->vertices[c1], mesh_->vertices[oc])));
}

std::unique_ptr<TetMeshOptimizer::Swap_edge_info>
TetMeshOptimizer::get_swap_edge_info(edge e, bool verbose) {
  // ChronoAuto chrono(__func__);

  std::unique_ptr<Swap_edge_info> swap_info =
      std::make_unique<Swap_edge_info>();
  swap_info->edge = e;
  swap_info->id = edge_to_size_t(e);

  if (mesh_->isOnBoundary(e.first, e.second, verbose)) {
    swap_info->is_good = false;
    return swap_info;
  }

  std::vector<tetrahedra> incident_tetrahedras;
  std::vector<vertex> one_ring;

  mesh_->OneRing(e, one_ring, incident_tetrahedras, verbose);

  double energy_of_swapping = DBL_MAX, current_energy_of_swapping,
         pre_transformation_energy = 0;
  edge opposite_edge;
  vertex min_collapse_point;
  uint32_t i_v1, i_v2, num_tet_seen, num_neigh_on_ring;
  std::vector<vertex> neighbour_collapse_point;

  for (vertex collapse_point : one_ring) {
    num_tet_seen = 0;
    num_neigh_on_ring = 0;
    current_energy_of_swapping = 0;
    neighbour_collapse_point.clear();
    mesh_->VV(collapse_point, neighbour_collapse_point);

    for (vertex neigh_collapse_point : neighbour_collapse_point) {
      for (vertex ring_vertex : one_ring) {
        if (neigh_collapse_point == ring_vertex) {
          num_neigh_on_ring++;
          break;
        }
      }
    }
    if (num_neigh_on_ring != 2) {
      current_energy_of_swapping = DBL_MAX;
    } else {

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
          t_v1_points[i_v1] = mesh_->vertices[collapse_point];
          i_v2 = mesh_->get_index_of_vertex_in_tet(e.second, t);
          t_v2_points[i_v2] = mesh_->vertices[collapse_point];
          if (!is_oriented_good(t_v1_points) ||
              !is_oriented_good(t_v2_points)) {
            current_energy_of_swapping = DBL_MAX;
            break;
          } else {
            current_energy_of_swapping =
                std::max(current_energy_of_swapping,
                         std::max(get_quality_measure(t_v1_points),
                                  get_quality_measure(t_v2_points)));
          }
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
  }

  return swap_info;
}

// 2-3 swap
// See blender file to understand indices (it's non sense otherwise)
// When it says a face it's a corner that represent the opposite face in its
// tetrahedra
void TetMeshOptimizer::swap_face(std::unique_ptr<Swap_face_info> swap_info,
                                 Swap_face_result *swap_result) {
  swap_result->success = false;

  corner face = swap_info->face;

  const uint64_t b2 = mesh_->tet_node.size();
  const size_t newsize = mesh_->tet_node.size() + 4;

  const uint64_t index_in_tet = face & 3;
  const corner tet_basis = face - index_in_tet;

  const uint64_t r0 = tet_basis + mesh_->tetON1(index_in_tet),
                 r1 = tet_basis + mesh_->tetON3(index_in_tet),
                 r2 = tet_basis + mesh_->tetON2(index_in_tet);
  const uint32_t c0 = mesh_->tet_node[r0], c1 = mesh_->tet_node[r1],
                 c2 = mesh_->tet_node[r2], c3 = mesh_->tet_node[face];

  const uint64_t g00 = mesh_->tet_neigh[r0], g01 = mesh_->tet_neigh[r1],
                 g02 = mesh_->tet_neigh[r2];

  const uint64_t orx = mesh_->tet_neigh[face];
  const uint64_t opp = orx & (~3);
  const uint64_t or0 = mesh_->tetCornerAtVertex(opp, c0);
  const uint64_t or1 = mesh_->tetCornerAtVertex(opp, c1);
  const uint64_t or2 = mesh_->tetCornerAtVertex(opp, c2);

  const uint64_t g10 = mesh_->tet_neigh[or0], g11 = mesh_->tet_neigh[or1],
                 g12 = mesh_->tet_neigh[or2];

  const uint32_t oc = mesh_->tet_node[orx];

  if (!is_oriented_good(mesh_->vertices[c3], mesh_->vertices[c3],
                        mesh_->vertices[c3], mesh_->vertices[c3]) ||
      !is_oriented_good(mesh_->vertices[c3], mesh_->vertices[c3],
                        mesh_->vertices[c3], mesh_->vertices[c3]) ||
      !is_oriented_good(mesh_->vertices[c3], mesh_->vertices[c3],
                        mesh_->vertices[c3], mesh_->vertices[c3]))
    return;

  mesh_->tet_node.resize(newsize);
  mesh_->tet_neigh.resize(newsize);
  mesh_->mark_tetrahedra.resize(newsize >> 2,
                                mesh_->mark_tetrahedra[tet_basis >> 2]);

  uint32_t *tn = mesh_->getTetNodes(tet_basis);
  *tn++ = c3;
  *tn++ = oc;
  *tn++ = c1;
  *tn++ = c2;
  tetrahedras_energy[mesh_->get_tetrahedra_index_from_corner(tet_basis)] =
      get_quality_measure(mesh_->get_tetrahedra_index_from_corner(tet_basis));
  tn = mesh_->getTetNodes(opp);
  *tn++ = c3;
  *tn++ = c0;
  *tn++ = oc;
  *tn++ = c2;
  tetrahedras_energy[mesh_->get_tetrahedra_index_from_corner(opp)] =
      get_quality_measure(mesh_->get_tetrahedra_index_from_corner(opp));
  tn = mesh_->getTetNodes(b2);
  *tn++ = c3;
  *tn++ = c0;
  *tn++ = c1;
  *tn++ = oc;
  tetrahedras_energy.push_back(
      get_quality_measure(mesh_->get_tetrahedra_index_from_corner(b2)));

  uint64_t *tg = mesh_->getTetNeighs(tet_basis);
  *tg++ = g10;
  *tg++ = g00;
  *tg++ = opp + 1;
  *tg++ = b2 + 1;
  tg = mesh_->getTetNeighs(opp);
  *tg++ = g11;
  *tg++ = tet_basis + 2;
  *tg++ = g01;
  *tg++ = b2 + 2;
  tg = mesh_->getTetNeighs(b2);
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

  for (int i = 0; i < 4; i++) {
    corner c = tet_basis + i;
    swap_result->impacted_faces.push_back(c);
    swap_result->impacted_faces.push_back(mesh_->tet_neigh[c]);
  }
  for (int i = 0; i < 4; i++) {
    corner c = opp + i;
    swap_result->impacted_faces.push_back(c);
    swap_result->impacted_faces.push_back(mesh_->tet_neigh[c]);
  }
  for (int i = 0; i < 4; i++) {
    corner c = b2 + i;
    swap_result->impacted_faces.push_back(c);
    swap_result->impacted_faces.push_back(mesh_->tet_neigh[c]);
  }

  mesh_->inc_tet[c0] = opp >> 2;
  mesh_->inc_tet[c1] = tet_basis >> 2;

  remove_duplicates(swap_result->impacted_faces);

  swap_result->success = true;
}

/// FOURTH PASS

bool TetMeshOptimizer::fourth_pass(
    double Move_info::*field_to_optimize,
    UpdatableQueue<Move_info, double, vertex> &queue) {
  if (verbose)
    std::cout << "\nStarting FOURTH pass" << std::endl;

  uint32_t n_vertices = mesh_->numVertices();

  for (vertex v = 0; v < n_vertices; v++) {
    std::unique_ptr<Move_info> move_info = get_move_info(v);
    if (move_info->is_good)
      queue.push(std::move(move_info), (*move_info).*field_to_optimize);
  }

  if (verbose)
    std::cout << "Determined " << queue.size() << " points to move"
              << std::endl;

  size_t num_vertex_moved = 0;

  while (!queue.empty()) {
    auto move_info = queue.pop();
    if (move_info.has_value()) {
      move_and_update(std::move(move_info.value()), queue, field_to_optimize);
      num_vertex_moved++;
    } else
      break;
  }

  if (verbose)
    std::cout << "Actually moved " << num_vertex_moved << " points"
              << std::endl;

  if (verbose)
    std::cout << "Finished FOURTH pass" << std::endl;
  return num_vertex_moved > 0;
}

void TetMeshOptimizer::move_and_update(
    std::unique_ptr<Move_info> move_info,
    UpdatableQueue<Move_info, double, vertex> &queue,
    double Move_info::*field_to_optimize) {

  auto move_result = std::make_unique<Move_result>();

  move_vertex(std::move(move_info), move_result.get());

  for (vertex vertex_to_update : move_result->impacted_vertices) {
    std::unique_ptr<Move_info> move_info = get_move_info(vertex_to_update);
    if (move_info->is_good)
      queue.set(std::move(move_info), (*move_info).*field_to_optimize);
    else
      queue.remove(vertex_to_update);
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
          barycenter += vector3d(mesh_->vertices[act_vertex]);
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
  move_info->id = v;

  if (mesh_->isOnBoundary(v)) {
    move_info->is_good = false;
    return move_info;
  }

  std::vector<tetrahedra> incident_tetrahedras;
  mesh_->VT(v, incident_tetrahedras);

  pointType *barycenter;
  try {
    barycenter = get_barycenter(v, incident_tetrahedras);
    move_info->coord_to_move_to = barycenter;
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
        std::max(post_transformation_energy, get_quality_measure(tet_points));
  }

  move_info->is_good = pre_transformation_energy > post_transformation_energy;

  if (move_info->is_good) {
    move_info->delta = pre_transformation_energy - post_transformation_energy;
    move_info->pre_energy = pre_transformation_energy;
  }

  return move_info;
}

void TetMeshOptimizer::move_vertex(std::unique_ptr<Move_info> move_info,
                                   Move_result *move_result) {
  move_result->success = false;

  vertex v = move_info->vertex;
  std::vector<tetrahedra> incident_tetrahedras;

  mesh_->VTfull(v, incident_tetrahedras);
  mesh_->VV(v, move_result->impacted_vertices);

  mesh_->vertices[v] = move_info->coord_to_move_to;

  for (tetrahedra t : incident_tetrahedras)
    tetrahedras_energy[t] = get_quality_measure(t);

  move_result->success = true;
}

/// FIFTH PASS

bool TetMeshOptimizer::fifth_pass(
    double Split_tetrahedra_info::*field_to_optimize,
    UpdatableQueue<Split_tetrahedra_info, double, tetrahedra> &queue) {
  if (verbose)
    std::cout << "\nStarting FIFTH pass" << std::endl;

  for (tetrahedra t = 0; t < mesh_->numTets(); t++) {
    auto split_info = get_split_tetrahedra_info(t);
    if (split_info->is_good) {
      queue.push(std::move(split_info), (*split_info).*field_to_optimize);
    }
  }
  if (verbose)
    std::cout << "Determined " << queue.size() << " tetrahedras to split"
              << std::endl;

  size_t num_tetrahedra_splitted = 0;
  std::vector<tetrahedra> impacted_tets;

  while (!queue.empty()) {
    auto split_info = queue.pop();
    if (split_info.has_value()) {
      split_tetrahedra_and_update(std::move(split_info.value()), queue,
                                  field_to_optimize);
      num_tetrahedra_splitted++;
    } else
      break;
  }

  temp_remap_vertex.clear();

  if (verbose) {
    std::cout << "Actually splitted " << num_tetrahedra_splitted
              << " tets\nFinished FIFTH pass" << std::endl;
  }

  return num_tetrahedra_splitted != 0;
}

std::unique_ptr<TetMeshOptimizer::Split_tetrahedra_info>
TetMeshOptimizer::get_split_tetrahedra_info(tetrahedra t) {
  // ChronoAuto chrono(__func__);

  std::unique_ptr<Split_tetrahedra_info> split_info =
      std::make_unique<Split_tetrahedra_info>();
  if (mesh_->has_infinite_vertex(t)) {
    split_info->is_good = false;
    return split_info;
  }

  split_info->tetrahedra = t;
  split_info->id = t;

  std::vector<tetrahedra> incident_tetrahedras, incident_tetrahedras_opposite;
  std::unordered_set<vertex> vertices_edge, vertices_opposite_edge;
  int size_of_intersection;

  edge current_edge, best_edge, current_opposite_edge;

  double pre_transformation_energy = 0., post_transformation_energy = DBL_MAX,
         current_post_transformation_energy, current_pre_transformation_energy;

  pointType *split_point;

  int num_inc_tet = 0;

  for (int i = 1; i < 4; i++) {
    for (int j = 0; j < i; j++) {

      current_post_transformation_energy = 0.;
      current_pre_transformation_energy = 0.;

      current_edge = make_edge(mesh_->get_i_th_vertex_of_tetrahedra(t, i),
                               mesh_->get_i_th_vertex_of_tetrahedra(t, j));

      mesh_->oppositeTetEdgePair(t, current_edge, current_opposite_edge);

      if (mesh_->isOnBoundary(current_opposite_edge.first,
                              current_opposite_edge.second)) {
        continue;
      }

      split_point = get_split_point(current_edge);

      incident_tetrahedras.clear();
      mesh_->ETfull(current_edge.first, current_edge.second,
                    incident_tetrahedras);
      num_inc_tet = 0;

      for (tetrahedra neigh_tet : incident_tetrahedras) {
        if (!mesh_->has_infinite_vertex(neigh_tet) && neigh_tet != t &&
            mesh_->mark_tetrahedra[t] == DT_IN) {
          num_inc_tet++;
          std::vector<pointType *> tet_points_1 =
              mesh_->getTetPoints(neigh_tet);
          std::vector<pointType *> tet_points_2 =
              mesh_->getTetPoints(neigh_tet);

          tet_points_1[mesh_->get_index_of_vertex_in_tet(
              current_edge.first, neigh_tet)] = split_point;
          tet_points_2[mesh_->get_index_of_vertex_in_tet(
              current_edge.second, neigh_tet)] = split_point;

          current_post_transformation_energy =
              std::max(current_post_transformation_energy,
                       std::max(get_quality_measure(tet_points_1),
                                get_quality_measure(tet_points_2)));
          current_pre_transformation_energy = std::max(
              current_pre_transformation_energy, tetrahedras_energy[neigh_tet]);
        } else {
          if (neigh_tet == t) {
            current_pre_transformation_energy =
                std::max(current_pre_transformation_energy,
                         tetrahedras_energy[neigh_tet]);
          }
        }
      }

      if (num_inc_tet <= 1) {
        continue;
      }

      incident_tetrahedras_opposite.clear();

      mesh_->ETfull(current_opposite_edge.first, current_opposite_edge.second,
                    incident_tetrahedras_opposite);

      size_of_intersection = 0;

      vertices_edge.clear();
      vertices_opposite_edge.clear();

      for (tetrahedra t : incident_tetrahedras) {
        for (int i = 0; i < 4; i++) {
          if (!vertices_edge.contains(
                  mesh_->get_i_th_vertex_of_tetrahedra(t, i))) {
            vertices_edge.insert(mesh_->get_i_th_vertex_of_tetrahedra(t, i));
          }
        }
      }

      for (tetrahedra t : incident_tetrahedras_opposite) {
        for (int i = 0; i < 4; i++) {
          if (vertices_edge.contains(
                  mesh_->get_i_th_vertex_of_tetrahedra(t, i)) &&
              !vertices_opposite_edge.contains(
                  mesh_->get_i_th_vertex_of_tetrahedra(t, i))) {
            size_of_intersection++;
            vertices_opposite_edge.insert(
                mesh_->get_i_th_vertex_of_tetrahedra(t, i));
          }
        }
      }

      if (size_of_intersection != 4) {
        current_post_transformation_energy = DBL_MAX;
        continue;
      }

      num_inc_tet = 0;

      for (tetrahedra neigh_tet : incident_tetrahedras_opposite) {
        if (!mesh_->has_infinite_vertex(neigh_tet) && neigh_tet != t &&
            mesh_->mark_tetrahedra[t] == DT_IN) {
          num_inc_tet++;
          std::vector<pointType *> tet_points_1 =
              mesh_->getTetPoints(neigh_tet);
          std::vector<pointType *> tet_points_2 =
              mesh_->getTetPoints(neigh_tet);

          tet_points_1[mesh_->get_index_of_vertex_in_tet(
              current_opposite_edge.first, neigh_tet)] = split_point;
          tet_points_2[mesh_->get_index_of_vertex_in_tet(
              current_opposite_edge.second, neigh_tet)] = split_point;

          current_post_transformation_energy =
              std::max(current_post_transformation_energy,
                       std::max(get_quality_measure(tet_points_1),
                                get_quality_measure(tet_points_2)));
          current_pre_transformation_energy = std::max(
              current_pre_transformation_energy, tetrahedras_energy[neigh_tet]);
        } else {
          if (neigh_tet == t) {
            current_pre_transformation_energy =
                std::max(current_pre_transformation_energy,
                         tetrahedras_energy[neigh_tet]);
          }
        }
      }

      if (num_inc_tet <= 1)
        continue;

      if (current_post_transformation_energy <
              current_pre_transformation_energy &&
          current_post_transformation_energy < post_transformation_energy) {
        post_transformation_energy = current_post_transformation_energy;
        pre_transformation_energy = current_pre_transformation_energy;
        best_edge = current_edge;
      }
    }
  }

  split_info->is_good = pre_transformation_energy != DBL_MAX &&
                        post_transformation_energy < pre_transformation_energy;

  if (split_info->is_good) {
    split_info->delta = pre_transformation_energy - post_transformation_energy;
    split_info->pre_energy = pre_transformation_energy;
    split_info->edge_to_split = best_edge;
  }

  return split_info;
}

void TetMeshOptimizer::split_tetrahedra_and_update(
    std::unique_ptr<Split_tetrahedra_info> split_info,
    UpdatableQueue<Split_tetrahedra_info, double, tetrahedra> &queue,
    double Split_tetrahedra_info::*field_to_optimize) {

  // ChronoAuto chrono(__func__);
  auto split_result = std::make_unique<Split_tetrahedra_result>();
  split_tetrahedra(std::move(split_info), split_result.get());

  for (tetrahedra t : split_result->removed_tetrahedras) {
    queue.remove(t);
    if (t < mesh_->numTets()) {
      std::unique_ptr<Split_tetrahedra_info> split_info =
          get_split_tetrahedra_info(t);
      if (split_info->is_good) {
        queue.set(std::move(split_info), (*split_info).*field_to_optimize);
      }
    }
  }

  for (tetrahedra t : split_result->impacted_tetrahedras) {
    queue.remove(t);
    if (t < mesh_->numTets()) {
      std::unique_ptr<Split_tetrahedra_info> split_info =
          get_split_tetrahedra_info(t);
      if (split_info->is_good) {
        queue.set(std::move(split_info), (*split_info).*field_to_optimize);
      }
    }
  }
}

void TetMeshOptimizer::split_tetrahedra(
    std::unique_ptr<Split_tetrahedra_info> split_info,
    Split_tetrahedra_result *split_result) {
  // ChronoAuto chrono(__func__);

  split_result->success = false;

  tetrahedra tet_to_split = split_info->tetrahedra;
  edge edge_to_split = split_info->edge_to_split;

  auto info = get_split_tetrahedra_info(tet_to_split);

  if (!info->is_good) {
    throw std::runtime_error("Check before splitting is not good for tet " +
                             std::to_string(tet_to_split));
  }

  edge opposite_edge, opposite_edge_update;

  mesh_->oppositeTetEdgePair(tet_to_split, edge_to_split, opposite_edge);

  std::unique_ptr<Split_edge_result> split_edge_result =
      std::make_unique<Split_edge_result>();

  split_edge(std::move(std::make_unique<Split_edge_info>(edge_to_split)),
             split_edge_result.get());

  for (tetrahedra t : split_edge_result->impacted_tetrahedras) {
    if (mesh_->tetHasVertex(t, edge_to_split.first)) {
      mesh_->oppositeTetEdgePair(
          t, make_edge(edge_to_split.first, mesh_->numVertices() - 1),
          opposite_edge_update);
      if (opposite_edge_update.first != INFINITE_VERTEX)
        mesh_->VTfull(opposite_edge_update.first,
                      split_result->impacted_tetrahedras);
      if (opposite_edge_update.second != INFINITE_VERTEX)
        mesh_->VTfull(opposite_edge_update.second,
                      split_result->impacted_tetrahedras);
    } else {
      if (mesh_->tetHasVertex(t, edge_to_split.second)) {
        mesh_->oppositeTetEdgePair(
            t, make_edge(edge_to_split.second, mesh_->numVertices() - 1),
            opposite_edge_update);
        if (opposite_edge_update.first != INFINITE_VERTEX)
          mesh_->VTfull(opposite_edge_update.first,
                        split_result->impacted_tetrahedras);
        if (opposite_edge_update.second != INFINITE_VERTEX)
          mesh_->VTfull(opposite_edge_update.second,
                        split_result->impacted_tetrahedras);
      }
    }
  }

  split_edge_result = std::make_unique<Split_edge_result>();

  split_edge(std::move(std::make_unique<Split_edge_info>(opposite_edge)),
             split_edge_result.get());

  for (tetrahedra t : split_edge_result->impacted_tetrahedras) {
    if (mesh_->tetHasVertex(t, opposite_edge.first)) {
      mesh_->oppositeTetEdgePair(
          t, make_edge(opposite_edge.first, mesh_->numVertices() - 1),
          opposite_edge_update);
      if (opposite_edge_update.first != INFINITE_VERTEX)
        mesh_->VTfull(opposite_edge_update.first,
                      split_result->impacted_tetrahedras);
      if (opposite_edge_update.second != INFINITE_VERTEX)
        mesh_->VTfull(opposite_edge_update.second,
                      split_result->impacted_tetrahedras);
    } else {
      if (mesh_->tetHasVertex(t, opposite_edge.second)) {
        mesh_->oppositeTetEdgePair(
            t, make_edge(opposite_edge.second, mesh_->numVertices() - 1),
            opposite_edge_update);
        if (opposite_edge_update.first != INFINITE_VERTEX)
          mesh_->VTfull(opposite_edge_update.first,
                        split_result->impacted_tetrahedras);
        if (opposite_edge_update.second != INFINITE_VERTEX)
          mesh_->VTfull(opposite_edge_update.second,
                        split_result->impacted_tetrahedras);
      }
    }
  }

  edge edge_to_collapse =
      make_edge(mesh_->numVertices() - 2, mesh_->numVertices() - 1);

  if (!link_condition(edge_to_collapse)) {
    throw std::runtime_error("Collapsing the edge in splitting tetrahedra will "
                             "cause break of link condition");
  }

  std::unique_ptr<Collapse_info> collapse_info =
      std::make_unique<Collapse_info>(edge_to_size_t(edge_to_collapse),
                                      edge_to_collapse,
                                      mesh_->numVertices() - 2);

  auto collapse_result = std::make_unique<Collapse_result>();

  collapse(std::move(collapse_info), collapse_result.get());

  split_result->removed_tetrahedras.reserve(
      split_result->removed_tetrahedras.size() +
      distance(collapse_result->removed_tetrahedras.begin(),
               collapse_result->removed_tetrahedras.end()));
  split_result->removed_tetrahedras.insert(
      split_result->removed_tetrahedras.end(),
      collapse_result->removed_tetrahedras.begin(),
      collapse_result->removed_tetrahedras.end());

  split_result->impacted_tetrahedras.reserve(
      split_result->impacted_tetrahedras.size() +
      distance(collapse_result->impacted_tetrahedras.begin(),
               collapse_result->impacted_tetrahedras.end()));
  split_result->impacted_tetrahedras.insert(
      split_result->impacted_tetrahedras.end(),
      collapse_result->impacted_tetrahedras.begin(),
      collapse_result->impacted_tetrahedras.end());

  remove_duplicates(split_result->impacted_tetrahedras);
  remove_duplicates(split_result->removed_tetrahedras);

  split_result->success = true;
}

void TetMeshOptimizer::get_edges_from_tetrahedras(
    std::vector<std::pair<uint32_t, uint32_t>> &all_edges,
    std::vector<tetrahedra> &tets) const {
  for (tetrahedra t : tets) {
    const uint32_t *tn = mesh_->tet_node.data() + (t << 2);
    if (tn[3] == INFINITE_VERTEX)
      continue;
    for (int i = 0; i < 4; i++) {
      for (int j = i + 1; j < 4; j++)
        if (tn[i] < tn[j])
          all_edges.push_back(std::pair<uint32_t, uint32_t>(tn[i], tn[j]));
        else
          all_edges.push_back(std::pair<uint32_t, uint32_t>(tn[j], tn[i]));
    }
  }
  std::sort(all_edges.begin(), all_edges.end());
  all_edges.erase(std::unique(all_edges.begin(), all_edges.end()),
                  all_edges.end());
}
void TetMeshOptimizer::get_incident_edges_from_tetrahedras(
    std::vector<std::pair<uint32_t, uint32_t>> &all_edges,
    std::vector<tetrahedra> &tets) const {
  std::vector<vertex> neighbour;
  std::unordered_set<vertex> explored_vertex;
  for (tetrahedra t : tets) {
    const uint32_t *tn = mesh_->tet_node.data() + (t << 2);
    if (tn[3] == INFINITE_VERTEX)
      continue;
    for (int i = 0; i < 4; i++) {
      if (!explored_vertex.contains(tn[i])) {
        neighbour.clear();
        mesh_->VV(tn[i], neighbour);
        for (vertex n : neighbour) {
          if (tn[i] < n)
            all_edges.push_back(std::pair<uint32_t, uint32_t>(tn[i], n));
          else
            all_edges.push_back(std::pair<uint32_t, uint32_t>(n, tn[i]));
        }
        explored_vertex.insert(tn[i]);
      }
    }
  }
  std::sort(all_edges.begin(), all_edges.end());
  all_edges.erase(std::unique(all_edges.begin(), all_edges.end()),
                  all_edges.end());
}

double TetMeshOptimizer::getTotalEnergy() {
  const uint32_t num_tets = mesh_->numTets();
  double total_energy = 0., current_tet_energy;
  for (uint32_t t = 0; t < num_tets; t++) {
    current_tet_energy = get_quality_measure(t);
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
      current_tet_energy = get_quality_measure(t);
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
      current_tet_energy = get_quality_measure(t);
      if (current_tet_energy != DBL_MAX) {
        num_real_tets++;
        mean_energy += current_tet_energy;
      }
    }
  }
  return mean_energy / num_real_tets;
}

std::string TetMeshOptimizer::get_energy_distribution() {
  std::ostringstream ss;
  for (tetrahedra t = 0; t < mesh_->numTets(); t++) {
    if (!mesh_->has_infinite_vertex(t) && mesh_->mark_tetrahedra[t] == DT_IN) {
      double rounded = std::ceil(tetrahedras_energy[t] * 100.0) / 100.0;
      ss << std::fixed << std::setprecision(2) << rounded << " ";
    }
  }
  return ss.str();
}
void TetMeshOptimizer::remove_tetrahedra(tetrahedra t) {
  tetrahedra last_tetrahedra = mesh_->numTets() - 1;
  if (t != last_tetrahedra) {
    mesh_->mark_tetrahedra[t] = mesh_->mark_tetrahedra[last_tetrahedra];
    tetrahedras_energy[t] = tetrahedras_energy[last_tetrahedra];
    for (int i = 0; i < 4; i++) {
      mesh_->tet_node[mesh_->get_i_th_corner_of_tetrahedra(t, i)] =
          mesh_->tet_node[mesh_->get_i_th_corner_of_tetrahedra(last_tetrahedra,
                                                               i)];
      mesh_->setMutualNeighbors(
          mesh_->tet_neigh[mesh_->get_i_th_corner_of_tetrahedra(last_tetrahedra,
                                                                i)],
          mesh_->get_i_th_corner_of_tetrahedra(t, i));
      if (mesh_->get_i_th_vertex_of_tetrahedra(last_tetrahedra, i) !=
              INFINITE_VERTEX &&
          mesh_->inc_tet[mesh_->get_i_th_vertex_of_tetrahedra(
              last_tetrahedra, i)] == last_tetrahedra)
        mesh_->inc_tet[mesh_->get_i_th_vertex_of_tetrahedra(last_tetrahedra,
                                                            i)] = t;
    }
    temp_remap_tetrahedra[last_tetrahedra] = t;
  }
  for (int i = 0; i < 4; i++)
    mesh_->tet_node.pop_back();
  for (int i = 0; i < 4; i++)
    mesh_->tet_neigh.pop_back();
  mesh_->mark_tetrahedra.pop_back();
  tetrahedras_energy.pop_back();
}

void TetMeshOptimizer::remove_vertex(vertex v) {
  vertex last_vertex = mesh_->numVertices() - 1;
  // std::cout << "Removing vertex " << v << " and remapping " << last_vertex
  //           << " to it" << std::endl;
  if (v != last_vertex) {
    mesh_->vertices[v] = mesh_->vertices[last_vertex];

    mesh_->marked_vertex[v] = mesh_->marked_vertex[last_vertex];

    mesh_->inc_tet[v] = mesh_->inc_tet[last_vertex];

    std::vector<tetrahedra> incident_tetrahedras;
    mesh_->VTfull(last_vertex, incident_tetrahedras);
    for (tetrahedra t : incident_tetrahedras)
      mesh_->tet_node[mesh_->get_corner_in_tet(t, last_vertex)] = v;

    temp_remap_vertex[last_vertex] = v;
  }

  mesh_->vertices.pop_back();
  mesh_->marked_vertex.pop_back();
  mesh_->inc_tet.pop_back();
}
/// VISUALISATION

void TetMeshOptimizer::register_tetrahedrisation(
    string mesh_name,
    const std::vector<std::vector<tetrahedra>> &highlighted_tetrahedras) {
  double mean_energy = getMeanEnergy();
  std::vector<std::array<double, 3>> coords;
  std::vector<std::array<size_t, 4>> tet_indices, tet_indices_out,
      inf_energy_tet_indices;
  std::vector<std::array<size_t, 4>> very_bad_tet_indices;
  std::vector<double> tet_energy, inf_tet_energy, out_tet_energy;
  std::vector<double> very_bad_tet_energy;
  double current_tet_energy;
  for (pointType *p : mesh_->vertices) {
    vector3d v(p);
    coords.push_back({v.c[0], v.c[1], v.c[2]});
  }
  for (corner i = 0; i < mesh_->tet_node.size(); i += 4) {

    if (i + 3 < mesh_->tet_node.size() &&
        !mesh_->has_infinite_vertex(
            mesh_->get_tetrahedra_index_from_corner(i))) {
      current_tet_energy =
          get_quality_measure(mesh_->vertices[mesh_->tet_node[i]],
                              mesh_->vertices[mesh_->tet_node[i + 1]],
                              mesh_->vertices[mesh_->tet_node[i + 2]],
                              mesh_->vertices[mesh_->tet_node[i + 3]]);
      if (current_tet_energy == DBL_MAX) {
        inf_energy_tet_indices.push_back(
            {mesh_->tet_node[i], mesh_->tet_node[i + 1], mesh_->tet_node[i + 2],
             mesh_->tet_node[i + 3]});
        inf_tet_energy.push_back(1e10);
      } else if (mesh_->mark_tetrahedra[i >> 2] == DT_OUT) {
        tet_indices_out.push_back({mesh_->tet_node[i], mesh_->tet_node[i + 1],
                                   mesh_->tet_node[i + 2],
                                   mesh_->tet_node[i + 3]});
        out_tet_energy.push_back(1e10);
      } else if (current_tet_energy > 2000) {
        very_bad_tet_indices.push_back(
            {mesh_->tet_node[i], mesh_->tet_node[i + 1], mesh_->tet_node[i + 2],
             mesh_->tet_node[i + 3]});
        very_bad_tet_energy.push_back(current_tet_energy);
      } else {
        tet_energy.push_back(current_tet_energy);
        tet_indices.push_back({mesh_->tet_node[i], mesh_->tet_node[i + 1],
                               mesh_->tet_node[i + 2], mesh_->tet_node[i + 3]});
      }
    }
  }

  int group = 1;
  std::vector<std::array<size_t, 4>> highlight_tets_indices;

  for (auto tets : highlighted_tetrahedras) {
    highlight_tets_indices.clear();
    for (tetrahedra t : tets) {
      highlight_tets_indices.push_back(
          {mesh_->get_i_th_vertex_of_tetrahedra(t, 0),
           mesh_->get_i_th_vertex_of_tetrahedra(t, 1),
           mesh_->get_i_th_vertex_of_tetrahedra(t, 2),
           mesh_->get_i_th_vertex_of_tetrahedra(t, 3)});
    }
    polyscope::registerTetMesh(mesh_name + " highlights " +
                                   std::to_string(group),
                               coords, highlight_tets_indices);
    group++;
  }

  polyscope::registerTetMesh(mesh_name, coords, tet_indices);
  polyscope::registerTetMesh(mesh_name + " very bad", coords,
                             very_bad_tet_indices);
  polyscope::registerTetMesh(mesh_name + " inf", coords,
                             inf_energy_tet_indices);
  polyscope::registerTetMesh(mesh_name + " out", coords, tet_indices_out);
  polyscope::getVolumeMesh(mesh_name)->addCellScalarQuantity(
      "Tetrahedra energy", tet_energy);
  polyscope::getVolumeMesh(mesh_name)->setEnabled(false);
  polyscope::getVolumeMesh(mesh_name + " very bad")->setEnabled(false);
  polyscope::getVolumeMesh(mesh_name + " very bad")
      ->addCellScalarQuantity("Bad tetrahedra energy", very_bad_tet_energy);
  polyscope::getVolumeMesh(mesh_name + " out")->setEnabled(false);
  polyscope::getVolumeMesh(mesh_name + " out")
      ->addCellScalarQuantity("Bad tetrahedra energy", out_tet_energy);
  polyscope::getVolumeMesh(mesh_name + " inf")->setEnabled(false);
  polyscope::getVolumeMesh(mesh_name + " inf")
      ->addCellScalarQuantity("Bad tetrahedra energy", inf_tet_energy);
}
