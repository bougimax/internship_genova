import polyscope as ps
import numpy as np
import sys

ps.init()
ps.set_transparency_mode("pretty")
ps.set_ground_plane_mode("none")


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
    vertices, inner_tets, outer_tets = [], [], []
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
                inner_tets.append([int(x) for x in line.split()[1:]])
            elif (
                num_vertices + num_inner_tets + 2
                < index_line
                <= num_vertices + num_inner_tets + num_outer_tets + 2
            ):
                # We are parsing an outer tetrahedra
                outer_tets.append([int(x) for x in line.split()[1:]])
            index_line += 1
    return vertices, inner_tets, outer_tets


def get_faces_from_tets(tets):
    faces = []
    for tet in tets:
        for i in range(4):
            faces.append(tet[:i] + tet[i + 1 :])
    return faces


def parse_splitted_edges(filepath):
    splitted_edges = set()
    with open(filepath, "r") as f:
        for line in f.readlines():
            splitted_edges.add(tuple(sorted(tuple(map(int, line.split())))))
    return splitted_edges


def parse_splitted_tetrahedras(filepath):
    splitted_tets = []
    with open(filepath, "r") as f:
        for line in f.readlines():
            splitted_tets.append(list(map(int, line.split())))
    return splitted_tets


def get_edges_from_faces(faces):
    seen_edges = set()
    edges = []
    for f in faces:
        n_vertices_faces = len(f)
        for i in range(n_vertices_faces):
            if tuple(sorted([f[i], f[(i + 1) % n_vertices_faces]])) not in seen_edges:
                edges.append(tuple(sorted([f[i], f[(i + 1) % n_vertices_faces]])))
                seen_edges.add(tuple(sorted([f[i], f[(i + 1) % n_vertices_faces]])))
    return edges


def add_splitted_edges_only(filepath):
    remap_vertex = {}
    new_vertices = []
    num_v = 0
    splitted_edges = list(parse_splitted_edges(filepath))
    for v1, v2 in splitted_edges:
        if v1 not in remap_vertex:
            remap_vertex[v1] = num_v
            new_vertices.append(output_vertices[v1])
            num_v += 1
        if v2 not in remap_vertex:
            remap_vertex[v2] = num_v
            new_vertices.append(output_vertices[v2])
            num_v += 1
    ps.register_curve_network(
        "Splitted edges only",
        np.array(new_vertices),
        np.array(
            list(map(lambda l: list(map(lambda x: remap_vertex[x], l)), splitted_edges))
        ),
    )


INPUT_FILEPATH = sys.argv[1]

input_vertices, input_faces = parse_off_file(INPUT_FILEPATH)

input = ps.register_surface_mesh(
    "Input",
    np.array(input_vertices),
    np.array(input_faces),
    enabled=False,
    edge_width=1,
)
output_vertices, output_inner_tets, output_outer_tets = parse_off_tet_file(
    INPUT_FILEPATH + ".tet"
)

output_inner = ps.register_volume_mesh(
    "Output inner",
    np.array(output_vertices),
    np.array(output_inner_tets),
    edge_width=0,
    transparency=0.05,
    color=(1, 1, 1),
    enabled=False,
)

output_edges = get_edges_from_faces(get_faces_from_tets(output_inner_tets))

splitted_edges = list(
    parse_splitted_edges(INPUT_FILEPATH + ".splitted_edges"),
)

output_inner_edges = ps.register_curve_network(
    "Output inner cn", np.array(output_vertices), np.array(output_edges), radius=0.0001
)

output_splited_edges = ps.register_curve_network(
    "Splitted edges only",
    np.array(output_vertices),
    np.array(splitted_edges),
    radius=0.0001,
)

add_splitted_edges_only(INPUT_FILEPATH + ".splitted_edges")

output_inner_edges.add_color_quantity(
    "Splitted edges",
    np.array(
        [(1, 0, 0) if e in splitted_edges else (0.8, 0.8, 0.8) for e in output_edges]
    ),
    defined_on="edges",
    enabled=True,
)

splitted_tets_ind = parse_splitted_tetrahedras(INPUT_FILEPATH + ".splitted_tetrahedras")

splitted_tets = ps.register_volume_mesh(
    "Splitted tets", np.array(output_vertices), np.array(splitted_tets_ind)
)


ps.show()
