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

// ----------------- Macros --------------------

using vertex = TetMesh::vertex;
using edge = TetMesh::edge;
using corner = TetMesh::corner;
using tetrahedra = TetMesh::tetrahedra;

// ----------------- End macros ----------------

class TetMeshOptimizer {
public:
  bool verbose = false;
  bool random = false;
  bool log = false;

  std::ofstream *time_log;
  std::ofstream *mean_energy_log;
  std::ofstream *max_energy_log;

  // Set functions
  void set_mesh(TetMesh &mesh) { mesh_ = &mesh; }
  void set_quality_measure(
      std::function<double(const pointType *, const pointType *,
                           const pointType *, const pointType *)>
          quality_measure) {
    quality_measure_ = quality_measure;
  };
  void set_verbose(bool v) { verbose = v; }
  void set_random(bool r) { random = r; }
  void set_log(bool l) { log = l; }

  // Global optimization function
  bool optimize();

  // Fills tetrahedras_energy with the energy of each tetrahedron
  void get_all_tets_energy();

  // Return the average quality measure across all the tetrahedra of the model
  double get_max_energy();

  // Return the maximum quality measure across all the tetrahedra of the model
  double get_mean_energy();

  // Register the current state of the tetrahedrization on Polyscope to
  // visualize it
  void register_tetrahedrisation(
      string mesh_name,
      const std::vector<std::vector<tetrahedra>> &higlighted_tetrahedras = {});

  // Return a string with all the quality measure separated by ' '
  std::string get_energy_distribution();

  TetMeshOptimizer() {};
  TetMeshOptimizer(TetMesh &mesh) { mesh_ = &mesh; };
  ~TetMeshOptimizer() {};

private:
  // Log functions to either log a message or log the average and maximum energy
  void log_message(std::string message);
  void log_energy();

  // The mesh structure, detailed in tetmesh.h/cpp
  TetMesh *mesh_;

  // The quality measure function, it can be swapped by another one, the
  // function requires to be increasing with the "bad quality" of the
  // tetrahedron
  std::function<double(const pointType *, const pointType *, const pointType *,
                       const pointType *)>
      quality_measure_ = energy_dirichlet;

  // Several functions to access more easily to the quality measure
  double get_quality_measure(const pointType *p1, const pointType *p2,
                             const pointType *p3, const pointType *p4);
  double get_quality_measure(const tetrahedra t);
  double get_quality_measure(const std::vector<pointType *> &tet_points);

  // Parameter to only look for optimizing the tetrahedra inside the model and
  // not all the tetrahedra in the convex hull of the model
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

  // Template structure for the information of a simulation of an operation
  struct Op_info {
    bool is_good; // if set to true the operation is good to perform, so it may
                  // be added in the queue
    double delta; // difference between the local maximum quality measure before
                  // and after the operation
    double pre_energy; // local maximum quality measure before the operation
  };

  // Template structure for the result of the execution of an operation
  struct Op_result {
    bool success; // return if the operation succeeded
  };

  // ------------- Inherited struct -----------------
  //
  // First pass:

  struct Split_edge_info : public Op_info {
    TetMesh::edge edge; // edge to be splitted
    size_t id; // id of the edge (it is the concatenation of the two vertices
               // represented by an 32 bits integer concatenated in a 64 bits
               // integer)
    pointType *split_point =
        nullptr; // position of the split point, by default set to nullptr which
                 // is understood as putting the point at the midpoint of the
                 // edge
    vertex split_vertex = INFINITE_VERTEX; // id of the split vertex, by default
                                           // set to INFINITE_VERTEX which is
                                           // understood as adding a new vertex

    // Constructor
    Split_edge_info(TetMesh::edge _e) : edge(_e) {}
    Split_edge_info() = default;
  };

  struct Split_edge_result : public Op_result {
    std::vector<edge> impacted_edges;
    std::vector<tetrahedra> impacted_tetrahedras;
  };

  // Second pass:
  struct Collapse_info : public Op_info {
    TetMesh::edge edge;     // edge to be collapsed
    size_t id;              // same id as in split_edge_info
    vertex collapse_vertex; // the vertex which the edge is going to be
                            // collapsed on

    // Constructor
    Collapse_info(size_t _id, TetMesh::edge _e, vertex _v)
        : id(_id), edge(_e), collapse_vertex(_v) {}
    Collapse_info() = default;
  };

  struct Collapse_result : public Op_result {
    std::vector<edge> impacted_edges;
    std::vector<edge> removed_edges;
    std::vector<tetrahedra> removed_tetrahedras;
    std::vector<tetrahedra> impacted_tetrahedras;
  };

  // Third pass:
  struct Swap_edge_info : public Op_info {
    TetMesh::edge edge;              // edge to be swapped
    size_t id;                       // same id as in split_edge_info
    TetMesh::vertex collapse_vertex; // the vertex on which the new vertex
                                     // will be collapsed on
  };

  struct Swap_edge_result : public Op_result {
    std::vector<edge> impacted_edges;
    std::vector<edge> removed_edges;
  };

  struct Swap_face_info : public Op_info {
    TetMesh::corner face; // face to be swapped
    TetMesh::corner id;   // it is the same integer
  };

  struct Swap_face_result : public Op_result {
    std::vector<corner> impacted_faces;
  };

  // Fourth pass:

  struct Move_info : public Op_info {
    TetMesh::vertex vertex;      // vertex to be moved
    TetMesh::vertex id;          // it is the same integer
    pointType *coord_to_move_to; // coordinate to move the vertex
  };

  struct Move_result : public Op_result {
    std::vector<vertex> impacted_vertices;
  };

  // Fifth pass:

  struct Collapse_tetrahedra_info : public Op_info {
    TetMesh::tetrahedra tetrahedra; // tetrahedra to be collapsed
    TetMesh::tetrahedra id;         // it is the same interger
    TetMesh::edge edge_to_split;    // the edge where the first split will occur
  };

  struct Collapse_tetrahedra_result : public Op_result {
    std::vector<tetrahedra> impacted_tetrahedras;
    std::vector<tetrahedra> removed_tetrahedras;
  };

  /// Default fields to optimize

  double Split_edge_info::*optim_field_edge_split = &Split_edge_info::delta;
  double Collapse_info::*optim_field_collapse = &Collapse_info::delta;
  double Swap_edge_info::*optim_field_edge_swap = &Swap_edge_info::delta;
  double Swap_face_info::*optim_field_face_swap = &Swap_face_info::delta;
  double Move_info::*optim_field_vertex_move = &Move_info::delta;
  double Collapse_tetrahedra_info::*optim_field_tet_split =
      &Collapse_tetrahedra_info::delta;

  void mark_vertices(std::vector<vertex> &vertices_to_mark, unsigned char mark);

  // FIRST PASS related functions:

  // Execute first pass
  bool first_pass(double Split_edge_info::*field_to_optimize,
                  UpdatableQueue<Split_edge_info, double, size_t> &queue);

  // Returns the maximum of the energy of the two tetrahedras that will produce
  // the split of the edge on the tetrahedra t that is adjacent to the edge to
  // split
  double get_energy_from_splitting(TetMesh::tetrahedra t,
                                   TetMesh::edge edge_to_split,
                                   pointType *potential_split_point);

  // Return the split point for an edge e, by default the split point is located
  // at the midpoint of the edge so t is set by default to 0.5, otherwise it
  // returns: t*e_1 + (1-t)*e_2
  pointType *get_split_point(TetMesh::edge e, double t = 0.5);

  // Returns if it's valid and worth to split the edge
  std::unique_ptr<Split_edge_info> get_split_edge_info(TetMesh::edge e);

  // Method to do the split and then to update the queue
  void
  split_edge_and_update(std::unique_ptr<Split_edge_info> split_info,
                        UpdatableQueue<Split_edge_info, double, size_t> &queue,
                        double Split_edge_info::*field_to_optimize);

  // Split the edge of the split info and put the result in split_result
  void split_edge(std::unique_ptr<Split_edge_info> split_info,
                  Split_edge_result *split_result);

  // SECOND PASS related functions:

  // Execute second pass
  bool second_pass(double Collapse_info::*field_to_optimize,
                   UpdatableQueue<Collapse_info, double, size_t> &queue);

  // Remap an edge
  void remap_edge(TetMesh::edge &edge) {
    while (temp_remap_vertex.contains(edge.first)) {
      edge.first = temp_remap_vertex[edge.first];
    }
    while (temp_remap_vertex.contains(edge.second)) {
      edge.second = temp_remap_vertex[edge.second];
    }
  }

  // Remap a vertex
  void remap_vertex(TetMesh::vertex &vertex) {
    while (temp_remap_vertex.contains(vertex)) {
      vertex = temp_remap_vertex[vertex];
    }
  }

  // Remap a tetrahedron
  void remap_tetrahedra(TetMesh::tetrahedra &tet) {
    while (temp_remap_tetrahedra.contains(tet)) {
      tet = temp_remap_tetrahedra[tet];
    }
  }

  // Return the maximum local quality measure after a supposed collapse of e
  double get_energy_from_collapsing(
      TetMesh::edge e,
      std::vector<TetMesh::tetrahedra> &incident_tetrahedras_v2,
      bool debug = false);

  // Return a Collapse_info with all the information detailed in the struct
  std::unique_ptr<Collapse_info> get_collapse_info(TetMesh::edge &edge,
                                                   bool debug = false);

  // Returns if the link condition is valid to collapse an edge
  // For e := v1--v2, it returns if one_ring(v1) \intersected one_ring(v2) =
  // one_ring(e) where one_ring(vertex) is the neighbour vertices of the vertex
  // and one_ring(edge) is the vertices that share a face with v1 and v2 (e)
  bool link_condition(TetMesh::edge e);

  // Collapse the edge of the collapse info and update the queue
  void collapse_and_update(std::unique_ptr<Collapse_info> collapse_info,
                           UpdatableQueue<Collapse_info, double, size_t> &queue,
                           double Collapse_info::*field_to_optimize);

  // Collapse the edge of collapse info
  void collapse(std::unique_ptr<Collapse_info> collapse_info,
                Collapse_result *collapse_result, bool debug = false);

  // THIRD PASS related functions:

  // Execute third pass
  bool third_pass(double Swap_edge_info::*field_to_optimize_edge,
                  double Swap_face_info::*field_to_optimize_face);

  // Execute the swap of the faces, here the only swap implemented is the 2-3
  // swap
  bool third_pass_face(double Swap_face_info::*field_to_optimize,
                       UpdatableQueue<Swap_face_info, double, corner> &queue);

  // Swap the face and update the queue
  void
  swap_face_and_update(std::unique_ptr<Swap_face_info> swap_info,
                       UpdatableQueue<Swap_face_info, double, corner> &queue,
                       double Swap_face_info::*field_to_optimize);

  // Return the information of the simulation of swapping face
  std::unique_ptr<Swap_face_info> get_swap_face_info(TetMesh::corner face);

  // Returns the maximum energy of the 3 new tetrahedras created by 2-3 swap
  // on the face identified by one of its opposite corner
  double get_energy_from_swapping_face(TetMesh::corner face);

  // Execute the swap
  void swap_face(std::unique_ptr<Swap_face_info> swap_info,
                 Swap_face_result *swap_result);

  // Execute the swap of the edges, it works by first splitting an edge then
  // collapse one of the edge created by the split vertex
  bool third_pass_edge(double Swap_edge_info::*field_to_optimize,
                       UpdatableQueue<Swap_edge_info, double, size_t> &queue);

  // Swap the edge and update the queue
  void
  swap_edge_and_update(std::unique_ptr<Swap_edge_info> swap_info,
                       UpdatableQueue<Swap_edge_info, double, size_t> &queue,
                       double Swap_edge_info::*field_to_optimize);

  // Swap the edge
  void swap_edge(std::unique_ptr<Swap_edge_info> swap_info,
                 Swap_edge_result *swap_result);

  // Returns the information of the simulation of swapping e
  std::unique_ptr<Swap_edge_info> get_swap_edge_info(TetMesh::edge e,
                                                     bool verbose = false);

  // FOURTH PASS related functions:

  // Execute fourth pass
  bool fourth_pass(double Move_info::*field_to_optimize,
                   UpdatableQueue<Move_info, double, vertex> &queue);

  // Move the vertex and update the queue
  void move_and_update(std::unique_ptr<Move_info> move_info,
                       UpdatableQueue<Move_info, double, vertex> &queue,
                       double Move_info::*field_to_optimize);

  // Returns the information of the simulation of moving vertex v
  std::unique_ptr<Move_info> get_move_info(TetMesh::vertex v);

  // Returns the barycenter of the neighbour of v
  pointType *
  get_barycenter(TetMesh::vertex v,
                 std::vector<TetMesh::tetrahedra> &incident_tetrahedras);

  // Move the vertex
  void move_vertex(std::unique_ptr<Move_info> move_info,
                   Move_result *move_result);

  // FIFTH PASS related functions:

  // Execute fifth pass
  bool fifth_pass(
      double Collapse_tetrahedra_info::*field_to_optimize,
      UpdatableQueue<Collapse_tetrahedra_info, double, tetrahedra> &queue);

  // Returns the infomation of the simulation of collapsing tetrahedron t
  std::unique_ptr<Collapse_tetrahedra_info>
  get_collapse_tetrahedra_info(TetMesh::tetrahedra t);

  // Collapse the tetrahedron and update the queue
  void collapse_tetrahedra_and_update(
      std::unique_ptr<Collapse_tetrahedra_info> split_info,
      UpdatableQueue<Collapse_tetrahedra_info, double, tetrahedra> &queue,
      double Collapse_tetrahedra_info::*field_to_optimize);

  // Collapse the tetrahedron
  void collapse_tetrahedra(std::unique_ptr<Collapse_tetrahedra_info> split_info,
                           Collapse_tetrahedra_result *split_result);

  // Fills all_edges with the edges of the tetrahedra contained in tets
  void get_edges_from_tetrahedras(std::vector<TetMesh::edge> &all_edges,
                                  std::vector<TetMesh::tetrahedra> &tets) const;

  // Fills all_edges with all the incident edges of all the vertices of the
  // tetrahedra contained in tets
  void get_incident_edges_from_tetrahedras(
      std::vector<TetMesh::edge> &all_edges,
      std::vector<TetMesh::tetrahedra> &tets) const;

  // Remove the vertex v of the mesh structure and update it accordingly
  void remove_vertex(vertex v);
  // Remove the tetrahedron of the mesh structure and update it accordingly
  void remove_tetrahedra(tetrahedra t);
};
