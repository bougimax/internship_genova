import polyscope as ps
import numpy as np
import sys

ps.init()


def parse_off_file(filepath):
    index_line = 0
    num_vertices, num_faces, num_edges = None, None, None
    vertices, faces = [], []
    is_off = False
    with open(filepath, "r") as obj:
        for line in obj.readlines():
            if (len(line) > 0 and line[0] == "#") or len(line) == 0:
                continue
            if "OFF" in line:
                is_off = True
                continue
            if index_line == 0:
                num_vertices, num_faces, num_edges = [int(x) for x in line.split()]
            elif index_line <= num_vertices:
                # We are parsing a vertex
                vertices.append([float(x) for x in line.split()])
            elif num_vertices < index_line <= num_vertices + num_faces:
                # We are parsing a face
                faces.append([int(x) for x in line.split()[1:]])
            index_line += 1
    return vertices, faces


def parse_off_tet_file(filepath):
    index_line = 0
    num_vertices, num_inner_tets, num_outer_tets = None, None, None
    vertices, inner_faces, outer_faces = [], [], []
    with open(filepath, "r") as obj:
        for line in obj.readlines():
            if (len(line) > 0 and line[0] == "#") or len(line) == 0:
                continue
            if index_line == 0:
                num_vertices = int(line.split()[0])
            elif index_line == 1:
                num_inner_tets = int(line.split()[0])
            elif index_line == 2:
                num_outer_tets = int(line.split()[0])
            elif index_line <= num_vertices + 2:
                # We are parsing a vertex
                vertices.append([float(x) for x in line.split()])
            elif num_vertices + 2 < index_line <= num_vertices + num_inner_tets + 2:
                # We are parsing an inner tetrahedra
                tetrahedra_vertices = [int(x) for x in line.split()[1:]]
                for i in range(4):
                    inner_faces.append(
                        tetrahedra_vertices[:i] + tetrahedra_vertices[i + 1 :]
                    )
            elif (
                num_vertices + num_inner_tets + 2
                < index_line
                <= num_vertices + num_inner_tets + num_outer_tets + 2
            ):
                # We are parsing an outer tetrahedra
                tetrahedra_vertices = [int(x) for x in line.split()[1:]]
                for i in range(4):
                    outer_faces.append(
                        tetrahedra_vertices[:i] + tetrahedra_vertices[i + 1 :]
                    )
            index_line += 1
    return vertices, inner_faces, outer_faces


INPUT_FILEPATH = sys.argv[1]

input_vertices, input_faces = parse_off_file(INPUT_FILEPATH)
input = ps.register_surface_mesh(
    "Input",
    np.array(input_vertices),
    np.array(input_faces),
    enabled=False,
    edge_width=1,
)
output_vertices, output_inner_faces, output_outer_faces = parse_off_tet_file(
    INPUT_FILEPATH + ".tet"
)
output_inner = ps.register_surface_mesh(
    "Output inner",
    np.array(output_vertices),
    np.array(output_inner_faces),
    edge_width=1,
)
output_outer = ps.register_surface_mesh(
    "Output outer",
    np.array(output_vertices),
    np.array(output_outer_faces),
    enabled=False,
    edge_width=1,
)
ps.show()
