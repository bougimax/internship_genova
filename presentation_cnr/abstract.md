# Constrained Delaunay tetrahedrisation refinement

The talk will focus on my internship at CNR IMATI.

Recently, Lorenzo Diazzi, Marco Attene, and colleagues published an article on Constrained Delaunay Tetrahedrization (CDT). Up to that point, they had developed a robust tetrahedrization pipeline using geometric predicates. However, they still lacked a refinement process for the tetrahedral output—something that other well-known software in the domain like TetGen or CGAL already offers.

CDTs have a wide range of applications. They are commonly used in physics-based simulations such as the finite element method, solutions to partial differential equations, or ray tracing. However, these simulations require high-quality meshes, which raw CDTs often fail to provide. Therefore, a refinement process is essential to improve overall mesh quality and ensure the robustness and accuracy of the simulations.

The objective of this internship was to design and validate a greedy algorithm for refining Constrained Delaunay Tetrahedrizations.

In this talk, I will present CDTs and the algorithm we designed.
