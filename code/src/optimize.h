#pragma once

#include "numeric_wrapper.h"
#include "quality_measure.h"
#include "random_map.h"
#include "tetmesh.h"
#include "updatable_priority_queue.h"
#include "updatable_queue_template.h"
#include "vector_3d.h"
#include <cfloat>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

inline std::unordered_map<std::string, long long> tempsCumules;

using vertex = TetMesh::vertex;
using edge = TetMesh::edge;
using corner = TetMesh::corner;
using tetrahedra = TetMesh::tetrahedra;

class ChronoAuto {
  std::string nom;
  std::chrono::high_resolution_clock::time_point debut;

public:
  ChronoAuto(const std::string &nomFonction) : nom(nomFonction) {
    debut = std::chrono::high_resolution_clock::now();
  }

  ~ChronoAuto() {
    auto fin = std::chrono::high_resolution_clock::now();
    auto duree =
        std::chrono::duration_cast<std::chrono::microseconds>(fin - debut)
            .count();
    if (tempsCumules.contains(nom))
      tempsCumules[nom] += duree;
    else
      tempsCumules[nom] = duree;
  }
};

class TetMeshOptimizer {
public:
  bool verbose = false;
  bool random = false;
  bool log = false;
  std::ofstream *time_log;
  std::ofstream *mean_energy_log;
  std::ofstream *max_energy_log;

  void set_mesh(TetMesh &mesh) { mesh_ = &mesh; }
  void set_quality_measure(
      std::function<double(const pointType *, const pointType *,
                           const pointType *, const pointType *)>
          quality_measure) {
    quality_measure_ = quality_measure;
  };
  bool optimize();
  void get_all_tets_energy();

  double getTotalEnergy();
  double getMaxEnergy();
  double getMeanEnergy();
  std::string get_energy_distribution();

  void register_tetrahedrisation(
      string mesh_name,
      const std::vector<std::vector<tetrahedra>> &higlighted_tetrahedras = {});

  TetMeshOptimizer() {};
  TetMeshOptimizer(TetMesh &mesh) { mesh_ = &mesh; };
  ~TetMeshOptimizer() {};

private:
  void log_message(std::string message);
  void log_energy();

  TetMesh *mesh_;
  std::function<double(const pointType *, const pointType *, const pointType *,
                       const pointType *)>
      quality_measure_ = energy_dirichlet;

  double get_quality_measure(const pointType *p1, const pointType *p2,
                             const pointType *p3, const pointType *p4);
  double get_quality_measure(const tetrahedra t);
  double get_quality_measure(const std::vector<pointType *> &tet_points);

  bool optimize_only_DT_IN = true;

  std::vector<double>
      tetrahedras_energy; // Contains the energy for each tetrahedra t necessary
                          // to compute optimization pass, it's initialized with
                          // get_all_tets_energy()
  std::map<TetMesh::vertex, TetMesh::vertex>
      temp_remap_vertex; // Temporary map for remapping when removing a vertex
  std::map<TetMesh::tetrahedra, TetMesh::tetrahedra>
      temp_remap_tetrahedra; // Temporary map for remapping when removing a
                             // tetrahedra

  struct Op_info {
    bool is_good;
    double delta;
    double pre_energy;
    int prioritize = 0;
  };

  struct Op_result {
    bool success;
  };

  struct Split_edge_info : public Op_info {
    TetMesh::edge edge;
    size_t id;
    pointType *split_point = nullptr;
    vertex split_vertex = INFINITE_VERTEX;

    Split_edge_info(TetMesh::edge _e) : edge(_e) {}
    Split_edge_info() = default;
  };

  struct Split_edge_result : public Op_result {
    std::vector<edge> impacted_edges;
    std::vector<tetrahedra> impacted_tetrahedras;
  };

  struct Collapse_info : public Op_info {
    size_t id;
    TetMesh::edge edge;
    vertex collapse_vertex;

    Collapse_info(size_t _id, TetMesh::edge _e, vertex _v)
        : id(_id), edge(_e), collapse_vertex(_v) {}

    Collapse_info() = default;
  };

  struct Collapse_result : public Op_result {
    std::vector<edge> impacted_edges;
    std::vector<edge> removed_edges;
    std::vector<tetrahedra> removed_tetrahedras;
    std::vector<tetrahedra> impacted_tetrahedras;
    std::unordered_set<edge, TetMesh::edge_hash> set_impacted_edges;
  };

  struct Swap_edge_info : public Op_info {
    size_t id;
    TetMesh::edge edge;
    TetMesh::vertex collapse_vertex;
  };

  struct Swap_edge_result : public Op_result {
    std::vector<edge> impacted_edges;
    std::vector<edge> removed_edges;
  };

  struct Swap_face_info : public Op_info {
    TetMesh::corner id;
    TetMesh::corner face;
  };

  struct Swap_face_result : public Op_result {
    std::vector<corner> impacted_faces;
  };

  struct Move_info : public Op_info {
    TetMesh::vertex id;
    TetMesh::vertex vertex;
    pointType *coord_to_move_to;
  };

  struct Move_result : public Op_result {
    std::vector<vertex> impacted_vertices;
  };

  struct Split_tetrahedra_info : public Op_info {
    TetMesh::tetrahedra tetrahedra;
    TetMesh::tetrahedra id;
    TetMesh::edge edge_to_split;
  };

  struct Split_tetrahedra_result : public Op_result {
    std::vector<tetrahedra> impacted_tetrahedras;
    std::vector<tetrahedra> removed_tetrahedras;
  };

  /// Default fields to optimize

  double Split_edge_info::*optim_field_edge_split = &Split_edge_info::delta;
  double Collapse_info::*optim_field_collapse = &Collapse_info::delta;
  double Swap_edge_info::*optim_field_edge_swap = &Swap_edge_info::delta;
  double Swap_face_info::*optim_field_face_swap = &Swap_face_info::delta;
  double Move_info::*optim_field_vertex_move = &Move_info::delta;
  double Split_tetrahedra_info::*optim_field_tet_split =
      &Split_tetrahedra_info::delta;

  void mark_vertices(std::vector<vertex> &vertices_to_mark, unsigned char mark);

  // FIRST PASS related functions:

  // Execute first pass (refining) of optimization process as described in
  // sec 3.2 of tetwild MAX
  bool first_pass(double Split_edge_info::*field_to_optimize,
                  UpdatableQueue<Split_edge_info, double, size_t> &queue);

  // Returns the maximum of the energy of the two tetrahedras that will produce
  // the split of the edge on the tetrahedra t that is adjacent to the edge to
  // split
  double get_energy_from_splitting(TetMesh::tetrahedra t,
                                   TetMesh::edge edge_to_split,
                                   pointType *potential_split_point);

  pointType *get_split_point(TetMesh::edge e, double t = 0.5);
  // Returns if it's valid and worth to split the edge
  std::unique_ptr<Split_edge_info> get_split_edge_info(TetMesh::edge e);

  void
  split_edge_and_update(std::unique_ptr<Split_edge_info> split_info,
                        UpdatableQueue<Split_edge_info, double, size_t> &queue,
                        double Split_edge_info::*field_to_optimize);

  // Split the edge edge_to_split with the vertex split_vertex in the middle
  void split_edge(std::unique_ptr<Split_edge_info> split_info,
                  Split_edge_result *split_result);

  // SECOND PASS related functions:

  // Execute second pass (coarsening) of optimization process as described in
  // sec 3.2 of tetwild MAX
  bool second_pass(double Collapse_info::*field_to_optimize,
                   UpdatableQueue<Collapse_info, double, size_t> &queue);

  // Remap and edge for the collapsing pass
  void remap_edge(TetMesh::edge &edge) {
    while (temp_remap_vertex.contains(edge.first)) {
      edge.first = temp_remap_vertex[edge.first];
    }
    while (temp_remap_vertex.contains(edge.second)) {
      edge.second = temp_remap_vertex[edge.second];
    }
  }
  void remap_vertex(TetMesh::vertex &vertex) {
    while (temp_remap_vertex.contains(vertex)) {
      vertex = temp_remap_vertex[vertex];
    }
  }

  void remap_tetrahedra(TetMesh::tetrahedra &tet) {
    while (temp_remap_tetrahedra.contains(tet)) {
      tet = temp_remap_tetrahedra[tet];
    }
  }

  double get_energy_from_collapsing(
      TetMesh::edge e,
      std::vector<TetMesh::tetrahedra> &incident_tetrahedras_v2,
      bool debug = false);

  // Returns a pair, the first is whether it's valid and worth to collapse the
  // edge, the second is on which end-vertices it should be collapsed (1 if on
  // edge.first, 2 if on edge.second)
  std::unique_ptr<Collapse_info> get_collapse_info(TetMesh::edge &edge,
                                                   bool debug = false);

  // Returns if the link condition is valid to collapse an edge
  // For e := v1--v2, it returns if one_ring(v1) \intersected one_ring(v2) =
  // one_ring(e) where one_ring(vertex) is the neighbour vertices of the vertex
  // and one_ring(edge) is the vertices that share a face with v1 and v2 (e)
  bool link_condition(TetMesh::edge e);
  void collapse_and_update(std::unique_ptr<Collapse_info> collapse_info,
                           UpdatableQueue<Collapse_info, double, size_t> &queue,
                           double Collapse_info::*field_to_optimize);

  // Collapse an edge onto its first endpoint
  void collapse(std::unique_ptr<Collapse_info> collapse_info,
                Collapse_result *collapse_result, bool debug = false);

  // THIRD PASS related functions:

  // Execute third pass (swapping) of optimization process as described in
  // sec 3.2 of tetwild MAX
  bool third_pass(double Swap_edge_info::*field_to_optimize_edge,
                  double Swap_face_info::*field_to_optimize_face);

  // Execute the swap of the faces, here the only swap implemented is the 2-3
  // swap
  bool third_pass_face(double Swap_face_info::*field_to_optimize,
                       UpdatableQueue<Swap_face_info, double, corner> &queue);
  void
  swap_face_and_update(std::unique_ptr<Swap_face_info> swap_info,
                       UpdatableQueue<Swap_face_info, double, corner> &queue,
                       double Swap_face_info::*field_to_optimize);

  // Returns if it's valid and worth to swap a face (the face is identified by
  // one of the two corner which is opposed to it)
  std::unique_ptr<Swap_face_info> get_swap_face_info(TetMesh::corner face);

  // Returns the maximum energy of the 3 new tetrahedras created by 2-3 swap
  // on the face identified by one of its opposite corner
  double get_energy_from_swapping_face(TetMesh::corner face);

  // 2-3 swap
  void swap_face(std::unique_ptr<Swap_face_info> swap_info,
                 Swap_face_result *swap_result);

  // Execute the swap of the edges, it works by first splitting an edge then
  // collapse one of the edge created by the split vertex
  bool third_pass_edge(double Swap_edge_info::*field_to_optimize,
                       UpdatableQueue<Swap_edge_info, double, size_t> &queue);

  void
  swap_edge_and_update(std::unique_ptr<Swap_edge_info> swap_info,
                       UpdatableQueue<Swap_edge_info, double, size_t> &queue,
                       double Swap_edge_info::*field_to_optimize);

  void swap_edge(std::unique_ptr<Swap_edge_info> swap_info,
                 Swap_edge_result *swap_result);

  // Returns a pair where the first is wether it's valid and worth to swap the
  // edge, the second is on which vertex we should collapse the edge created
  // by the splitting
  std::unique_ptr<Swap_edge_info> get_swap_edge_info(TetMesh::edge e,
                                                     bool verbose = false);

  // FOURTH PASS related functions:

  // Execute fourth pass (smoothing) of optimization process as described in
  // sec 3.2 of tetwild MAX
  bool fourth_pass(double Move_info::*field_to_optimize,
                   UpdatableQueue<Move_info, double, vertex> &queue);
  void move_and_update(std::unique_ptr<Move_info> move_info,
                       UpdatableQueue<Move_info, double, vertex> &queue,
                       double Move_info::*field_to_optimize);

  // Returns a pair where the first is wether it's valid and worth to move the
  // vertex to its center of mass, the second is the center of mass
  std::unique_ptr<Move_info> get_move_info(TetMesh::vertex v);

  pointType *
  get_barycenter(TetMesh::vertex v,
                 std::vector<TetMesh::tetrahedra> &incident_tetrahedras);
  void move_vertex(std::unique_ptr<Move_info> move_info,
                   Move_result *move_result);

  // FIFTH PASS related functions:
  bool
  fifth_pass(double Split_tetrahedra_info::*field_to_optimize,
             UpdatableQueue<Split_tetrahedra_info, double, tetrahedra> &queue);

  std::unique_ptr<Split_tetrahedra_info>
  get_split_tetrahedra_info(TetMesh::tetrahedra t);

  void split_tetrahedra_and_update(
      std::unique_ptr<Split_tetrahedra_info> split_info,
      UpdatableQueue<Split_tetrahedra_info, double, tetrahedra> &queue,
      double Split_tetrahedra_info::*field_to_optimize);

  void split_tetrahedra(std::unique_ptr<Split_tetrahedra_info> split_info,
                        Split_tetrahedra_result *split_result);

  void get_edges_from_tetrahedras(std::vector<TetMesh::edge> &all_edges,
                                  std::vector<TetMesh::tetrahedra> &tets) const;
  void get_incident_edges_from_tetrahedras(
      std::vector<TetMesh::edge> &all_edges,
      std::vector<TetMesh::tetrahedra> &tets) const;

  void remove_vertex(vertex v);
  void remove_tetrahedra(tetrahedra t);
};
