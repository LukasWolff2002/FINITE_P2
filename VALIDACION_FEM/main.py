import gmsh
import meshio
import numpy as np
import matplotlib.pyplot as plt

from nodes import Node
from elements import Element
from solve import Solve


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
                node.boundary_label.append(name)

 

    # Crear elementos genéricos (filtrando solo elementos de tipo 'triangle' o 'quad')
    # Crear elementos Quad9
    elements = []
    for cell_block in mesh.cells:
        if cell_block.type == "quad9":
            for i, node_ids in enumerate(cell_block.data):
                # Convertir los nodos de base 0 a base 1
                node_ids = [int(id) + 1 for id in node_ids]  # +1 para pasar a base 1
                elements.append(Element(i + 1, node_ids, n_nodes))  # Crear el elemento Quad9


    return nodes, elements

def main(alpha):   

    Estructure = None
    nodes = None
    elements = None

    geo_file = "GMSH_FILES/validacion_fem.geo"
    mesh_file = "GMSH_FILES/validacion_fem.msh"

    #Genero la malla
    nodes, elements = fixed_load_mesh_objects(geo_file=geo_file, msh_file=mesh_file, n_nodes=9)


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

    error = np.abs(solucion_analitica - solucion_fem)

    # Guardar en .txt

    with open("VALIDACION_FEM/resultados.txt", "a") as f:
        f.write(f"N = {1}, R = {2}, alpha = {alpha}\n")
        f.write(f"Error: {error:.6e}\n")
        f.write("-" * 40 + "\n")

    print("Resultados guardados en resultados.txt")

    

if __name__ == "__main__":

    open("VALIDACION_FEM/resultados.txt", "w").close()
    alpha = 5
    
    main(alpha)
   
       
        
    

    
    
    
    
