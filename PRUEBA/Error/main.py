import gmsh
import meshio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import os

from nodes import Node
from elements import Element
from solve import Solve


def fixed_load_mesh_objects(geo_file="geo.geo", msh_file="mesh.msh", n_nodes=4):
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
        print(f"Tipo de celda: {cell_block.type}")
        if cell_block.type == "quad": #Ojjo aqui
            for i, node_ids in enumerate(cell_block.data):
                # Convertir los nodos de base 0 a base 1
                node_ids = [int(id) + 1 for id in node_ids]  # +1 para pasar a base 1
                elements.append(Element(i + 1, node_ids, n_nodes))  # Crear el elemento Quad9

    print(elements)
    return nodes, elements

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

    plt.savefig("GRAFICOS/Validacion_FEM/solucion_comparativa_2.png")


def main(alpha):   

    Estructure = None
    nodes = None
    elements = None

    geo_file = "LST/geo.geo"
    mesh_file = "GMSH_FILES/Test2.msh"

    #En primero lugar modifico el archivo geo
    #modificar_geo(geo_file, geo_file, N, R)

    #Genero la malla
    nodes, elements = fixed_load_mesh_objects(geo_file=geo_file, msh_file=mesh_file)

    #Obtengo la solucion numerica por nodo
    for node in nodes:
        node.solve_u(alpha)

    #Resuelvo la estructura
    Estructure = Solve(nodes, elements, alpha)

    Estructure.solve_matrix()

    #errores = error(nodes)
    solucion_analitica = Estructure.semi_norm_H1_0(alpha)
    print(f"Solución analítica: {solucion_analitica}")
    solucion_fem = Estructure.femm_solution()
    print(f"Solución FEM: {solucion_fem}")


    plot_solution_3d(nodes)

    # Guardar en .txt

    

    

if __name__ == "__main__":

    main(alpha=3)
   
       
        
    

    
    
    
    
