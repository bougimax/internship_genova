import sys


def convert_obj_to_off(input_filepath, output_filepath):
    with open(input_filepath, "r") as input_file:
        with open(output_filepath, "w") as output_file:
            output_file.write("OFF\n")
            vertices, faces = [], []
            for line in input_file.readlines():
                if line[0:2] == "v ":
                    vertices.append(list(map(float, line[2:].split())))
                elif line[0:2] == "f ":
                    current_face = list(
                        map(
                            lambda x: int(x) - 1,
                            map(lambda x: x.split("/")[0], line[2:].split()),
                        )
                    )
                    if len(current_face) == 3:
                        faces.append([3] + current_face)
                    elif len(current_face) == 4:
                        faces.append([3] + current_face[:3])
                        faces.append([3] + current_face[2:] + [current_face[0]])
            output_file.write(f"{len(vertices)} {len(faces)} 0\n")
            for x in vertices + faces:
                output_file.write(" ".join(map(str, x)) + "\n")


input_fp = sys.argv[1]
output_fp = input_fp[:-3] + "off"

convert_obj_to_off(input_fp, output_fp)
