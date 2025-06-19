import gmsh
import meshio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from nodes import Node
from elements import Element
from solve import Solve

"""
def fixed_load_mesh_objects(geo_file="geo.geo", msh_file="mesh.msh", n_nodes=9):
    # Leer malla
    mesh = meshio.read(msh_file)

    # Crear nodos
    nodes = [Node(i + 1, x, y) for i, (x, y, _) in enumerate(mesh.points)]
    print(f"Se han creado {len(nodes)} nodos")

    # Obtener nodos de borde con etiquetas físicas
    boundary_nodes = {}
    gmsh.initialize()
    gmsh.open(msh_file)

    physicals = gmsh.model.getPhysicalGroups(1)
    name_map = {}
    for dim, tag in physicals:
        name = gmsh.model.getPhysicalName(dim, tag)
        name_map[tag] = name
        boundary_nodes[name] = set()

        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        for entity in entities:
            _, _, node_tags = gmsh.model.mesh.getElements(dim, entity)
            for nlist in node_tags:
                for node_id in nlist:
                    boundary_nodes[name].add(int(node_id))

    gmsh.finalize()

    # Asignar etiquetas de borde a los nodos
    for node in nodes:
        node.boundary_label = []
        for name, id_set in boundary_nodes.items():
            if node.id in id_set:
                print(name)
                node.boundary_label.append(name)

    # Crear elementos genéricos (filtrando solo elementos de tipo 'triangle' o 'quad')
    # Crear elementos Quad9
    elements = []
    for cell_block in mesh.cells:
        if cell_block.type == "quad9": #Ojjo aqui
            for i, node_ids in enumerate(cell_block.data):
                # Convertir los nodos de base 0 a base 1
                node_ids = [int(id) + 1 for id in node_ids]  # +1 para pasar a base 1
                elements.append(Element(i + 1, node_ids, n_nodes))  # Crear el elemento Quad9

    return nodes, elements
"""


def fixed_load_mesh_objects(geo_file="geo.geo", msh_file="mesh.msh", n_nodes=7):
    nodes = []
    # Coordenadas en rejilla 3x3 para subdividir el cuadrado [0,1]x[0,1]
    # y nodos intermedios
    coords = [
        # Esquinas (vértices)
        (1, 0.0, 0.0),   # A
        (2, 1.0, 0.0),   # B
        (3, 0.5, 0.5),   # C
        (4, 0.0, 1),   # D
        (5, 1.0, 1.0),   # E
       
    ]

    # Crear nodos
    added_ids = set()
    for nid, x, y in coords:
        if nid not in added_ids:
            node = Node(nid, x, y)
            if x in [0.0, 1.0] or y in [0.0, 1.0]:
                node.boundary_label.append("Dirichlet")
            nodes.append(node)
            added_ids.add(nid)

    # Crear elementos Quad7: [N1, N2, N3, N4, N5, N6, N7]
    elements = []

    # Elemento 1: inferior izquierdo
    elements.append(Element(1, [1, 2, 3, 4], n_nodes))
    elements.append(Element(2, [2, 5, 4, 3], n_nodes))

   
    return nodes, elements



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_solution_3d(nodes):
    """
    Esta función plotea la solución FEM y la solución analítica en 3D usando los nodos.
    La solución analítica es node.u y la solución FEM es node.u_fem.
    """
    # Extraer las coordenadas y las soluciones de desplazamiento para X, Y, Z
    x_coords = [node.x for node in nodes]
    y_coords = [node.y for node in nodes]
    z_fem = [node.u_fem for node in nodes]  # Solución FEM
    z_analytical = [node.u for node in nodes]  # Solución analítica

    # Configurar la figura para la visualización 3D
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Graficar la solución FEM
    ax.scatter(x_coords, y_coords, z_fem, c='b', marker='o', label='Solución FEM')

    # Graficar la solución analítica
    ax.scatter(x_coords, y_coords, z_analytical, c='r', marker='^', label='Solución Analítica')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Solucion (u)')

    ax.legend()
    ax.set_title('Comparación de Soluciones FEM y Analítica')

    plt.savefig("GRAFICOS/Validacion_FEM/solucion_comparativa.png")



def main(alpha):   
    Estructure = None
    nodes = None
    elements = None

    geo_file = "GMSH_FILES/validacion_fem.geo"
    mesh_file = "GMSH_FILES/validacion_fem.msh"

    # Genero la malla
    nodes, elements = fixed_load_mesh_objects(geo_file=geo_file, msh_file=mesh_file, n_nodes=4)

    # Obtengo la solución numérica por nodo
    for node in nodes:
        node.solve_u(alpha)

    # Resuelvo la estructura
    Estructure = Solve(nodes, elements, alpha)
    Estructure.solve_matrix()

    # Solución analítica
    solucion_analitica = Estructure.semi_norm_H1_0(alpha)

    print(f"Solución analítica: {solucion_analitica}")

    # Solución FEM
    solucion_fem = Estructure.femm_solution()
    print(f"Solución FEM: {solucion_fem}")

    error = np.abs(solucion_analitica - solucion_fem)

    # Graficar las soluciones
    plot_solution_3d(nodes)  # Llamada a la función de graficado

    # Guardar en .txt
    with open("VALIDACION_FEM/resultados.txt", "a") as f:
        f.write(f"N = {1}, R = {2}, alpha = {alpha}\n")
        f.write(f"Error: {error:.6e}\n")
        f.write("-" * 40 + "\n")

    print("Resultados guardados en resultados.txt")


if __name__ == "__main__":
    open("VALIDACION_FEM/resultados.txt", "w").close()
    alpha = 3
    
    main(alpha)
