import re


def remove_comments(code):
    # Remove // and /* */ comments
    code_without_comments = re.sub(r"//.*", "", code)
    code_without_comments = re.sub(
        r"/\*.*?\*/", "", code_without_comments, flags=re.DOTALL
    )
    return code_without_comments


def extract_methods(cpp_file, class_name):
    with open(cpp_file, "r", encoding="utf-8") as f:
        content = f.read()

    content = remove_comments(content)

    # Find the class block
    class_pattern = rf"class\s+{class_name}\s*{{(.*?)\n}};"
    match = re.search(class_pattern, content, re.DOTALL)

    if not match:
        print(f"Class '{class_name}' not found.")
        return []

    class_body = match.group(1)

    # Regex to detect declared methods
    method_pattern = r"(?:virtual\s+)?(?:static\s+)?[\w:<>\s*&]+?\s+(\w+)\s*\([^;]*\)\s*(?:const)?\s*;"
    methods = re.findall(method_pattern, class_body)

    forbidden_words = {"if", "for", "while", "return", "switch", "else"}
    methods = [m for m in methods if m not in forbidden_words]

    return methods


def prefix_method_calls(code, method_names):
    # Pattern : on cherche les methodes, mais on s'assure qu'il n'y a PAS mesh_-> ou :: juste avant
    # (?<!mesh_->) et (?<!::) sont des lookbehind négatifs
    pattern = (
        r"(?<!mesh_->)(?<!::)\b("
        + "|".join(re.escape(m) for m in method_names)
        + r")\s*"
    )

    def replacer(match):
        method = match.group(1)
        return f"mesh_->{method}"

    return re.sub(pattern, replacer, code)


FILEPATH_ORIGIN = "../src/tetmesh.h"
FILEPATH_TO_MODIFY = "delaunay.cpp"
FILEPATH_H = "delaunay.h"

with open(FILEPATH_TO_MODIFY, "r", encoding="utf-8") as f:
    contenu = f.read()

# Méthodes à remplacer (extraites avec le script précédent)
methodes = extract_methods(FILEPATH_ORIGIN, "TetMesh")
methodes_to_delete = extract_methods(FILEPATH_H, "TetMeshTetrahedrizer")
methodes = [x for x in methodes if x not in methodes_to_delete] + [
    "getTetNodes",
    "getTetNeighs",
    "resizeTets",
    "vertices",
    "inc_tet",
    "tet_node" "tet_neigh",
    "mark_tetrahedra",
    "marked_vertex",
    "has_outer_vertices",
]
# Appliquer le remplacement
code_modifie = prefix_method_calls(contenu, methodes)

# Sauvegarde (optionnelle)
with open(
    FILEPATH_TO_MODIFY.split(".")[0] + "_modified." + FILEPATH_TO_MODIFY.split(".")[1],
    "w",
    encoding="utf-8",
) as f:
    f.write(code_modifie)

print("Méthodes remplacées avec mesh_->")
