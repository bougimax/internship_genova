#include "delaunay.h"
#include "numeric_wrapper.h"
#include <vector>

class TetMeshOptimizer {
public:
  void set_mesh(TetMesh &mesh) { mesh_ = mesh; }
  void optimize();

private:
  TetMesh &mesh_;
  bool verbose = false;
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

  struct Split_edge_info : public Op_info {
    TetMesh::edge edge;
    std::map<TetMesh::tetrahedra, double> energy_per_tet;
    pointType *split_point;
  };
  struct Split_tetrahedra_info : public Op_info {
    TetMesh::tetrahedra tetrahedra;
    TetMesh::edge edge_1;
    pointType *split_point_1;
    TetMesh::edge edge_2;
    pointType *split_point_2;
  };
  struct Collapse_info : public Op_info {
    TetMesh::edge edge;
  };
  struct Swap_edge_info : public Op_info {
    TetMesh::edge edge;
    TetMesh::vertex collapse_vertex;
  };
  struct Swap_face_info : public Op_info {
    std::pair<TetMesh::corner, TetMesh::corner> face;
  };
  struct Move_info : public Op_info {
    TetMesh::vertex vertex;
    pointType *barycenter;
  };

  // FIRST PASS related functions:

  // Execute first pass (refining) of optimization process as described in
  // sec 3.2 of tetwild MAX
  void first_pass();

  // Returns the maximum of the energy of the two tetrahedras that will produce
  // the split of the edge on the tetrahedra t that is adjacent to the edge to
  // split
  double get_energy_from_splitting(TetMesh::tetrahedra t,
                                   TetMesh::edge edge_to_split,
                                   pointType *potential_split_point);

  pointType *get_split_point(TetMesh::edge e);
  // Returns if it's valid and worth to split the edge
  std::unique_ptr<Split_edge_info> get_split_edge_info(TetMesh::edge e);

  void split_edge_and_update(
      TetMesh::edge e,
      better_priority_queue::updatable_priority_queue<size_t, double> &queue);

  // Split the edge edge_to_split with the vertex split_vertex in the middle
  void split_edge(TetMesh::edge edge_to_split, TetMesh::vertex split_vertex,
                  std::vector<TetMesh::tetrahedra> &impacted_tetrahedras);

  // SECOND PASS related functions:

  // Execute second pass (coarsening) of optimization process as described in
  // sec 3.2 of tetwild MAX
  void second_pass();

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
      std::vector<TetMesh::tetrahedra> &incident_tetrahedras_v2);

  // Returns a pair, the first is whether it's valid and worth to collapse the
  // edge, the second is on which end-vertices it should be collapsed (1 if on
  // edge.first, 2 if on edge.second)
  std::unique_ptr<Collapse_info> get_collapse_info(TetMesh::edge &edge);

  // Returns if the link condition is valid to collapse an edge
  // For e := v1--v2, it returns if one_ring(v1) \intersected one_ring(v2) =
  // one_ring(e) where one_ring(vertex) is the neighbour vertices of the vertex
  // and one_ring(edge) is the vertices that share a face with v1 and v2 (e)
  bool link_condition(TetMesh::edge e);
  void collapse_on_v1_and_update(
      TetMesh::edge edge_to_collapse,
      better_priority_queue::updatable_priority_queue<size_t, double> &queue);

  // Collapse an edge onto its first endpoint
  bool collapse_on_v1(TetMesh::edge e,
                      std::vector<TetMesh::tetrahedra> &impacted_tetrahedras,
                      std::vector<TetMesh::edge> &removed_edges);

  // THIRD PASS related functions:

  // Execute third pass (swapping) of optimization process as described in
  // sec 3.2 of tetwild MAX
  void third_pass();

  // Execute the swap of the faces, here the only swap implemented is the 2-3
  // swap
  void third_pass_face();
  void swap_face_and_update(
      corner face_to_swap,
      better_priority_queue::updatable_priority_queue<corner, double> &queue);

  // Returns if it's valid and worth to swap a face (the face is identified by
  // one of the two corner which is opposed to it)
  std::unique_ptr<Swap_face_info> get_swap_face_info(corner face);

  // Returns the maximum energy of the 3 new tetrahedras created by 2-3 swap
  // on the face identified by one of its opposite corner
  double get_energy_from_swapping_face(corner face);

  // 2-3 swap
  bool swap_face(corner face, std::vector<corner> &impacted_faces,
                 bool prevent_inversion = true, double min_energy = DBL_MAX);

  // Execute the swap of the edges, it works by first splitting an edge then
  // collapse one of the edge created by the split vertex
  void third_pass_edge();

  void
  swap_edge_and_update(TetMesh::edge edge_to_swap,
                       TetMesh::vertex collapse_vertex,
                       better_priority_queue::updatable_priority_queue<
                           size_t, std::pair<double, TetMesh::vertex>> &queue);

  void swap_edge(TetMesh::edge edge_to_swap, TetMesh::vertex collapse_vertex,
                 std::vector<TetMesh::edge> &impacted_edges);

  // Returns a pair where the first is wether it's valid and worth to swap the
  // edge, the second is on which vertex we should collapse the edge created
  // by the splitting
  std::unique_ptr<Swap_edge_info> get_swap_edge_info(TetMesh::edge e);

  // FOURTH PASS related functions:

  // Execute fourth pass (smoothing) of optimization process as described in
  // sec 3.2 of tetwild MAX
  void fourth_pass();
  void move_and_update(
      TetMesh::vertex v,
      better_priority_queue::updatable_priority_queue<TetMesh::vertex, double>
          &queue);

  // Returns a pair where the first is wether it's valid and worth to move the
  // vertex to its center of mass, the second is the center of mass
  std::unique_ptr<Move_info> get_move_info(TetMesh::vertex v);

  pointType *
  get_barycenter(TetMesh::vertex v,
                 std::vector<TetMesh::tetrahedra> &incident_tetrahedras);
  void move_vertex(TetMesh::vertex v, pointType *coord_to_move);

  // FIFTH PASS related functions:
  void fifth_pass();

  std::unique_ptr<TetMesh::Split_tetrahedra_info>
  get_split_tetrahedra_info(TetMesh::tetrahedra t);

  bool try_to_split_tetrahedra(std::unique_ptr<Split_tetrahedra_info> info);

  void get_edges_from_tetrahedras(std::vector < std::pair<TetMesh::edge> &
                                      all_edges,
                                  std::vector<TetMesh::tetrahedra> &tets) const;

  double getTetEnergy(uint64_t t) const;

  double getTotalEnergy();
  double getMaxEnergy();
  double getMeanEnergy();

  void get_all_tets_energy();
};
}
