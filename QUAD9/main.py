import os
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import meshio
import matplotlib.tri as mtri
import re

from nodes import Node
from material import Material
from membrane import Membrane
from Quad2D import Quad9
from solve import Solve
from graph import plot_results
from collections import defaultdict

def make_nodes_groups(output_file, restriccion_x = None, restriccion_y = None, restriccion_xy = None):
    mesh = meshio.read(output_file)
    
    tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}
    grupos = defaultdict(dict)  # nombre_grupo: {id_nodo: Node}

    # Procesar elementos tipo quad9
    for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
        if cell_block.type != "quad9":
            continue
        for quad, tag in zip(cell_block.data, phys_tags):
            nombre = tag_to_name.get(tag, str(tag))
            for node_id in quad:
                x, y = mesh.points[node_id][:2]
                if node_id not in grupos[nombre]:
                    grupos[nombre][node_id] = Node(node_id + 1, [x, y])

    # Procesar líneas tipo line3 para condiciones de borde
    for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
        if cell_block.type != "line3":
            continue
        for line, tag in zip(cell_block.data, phys_tags):
            nombre = tag_to_name.get(tag, str(tag))
            for node_id in line:
                x, y = mesh.points[node_id][:2]
                restrain = [0, 0]
                if restriccion_x is not None and nombre in restriccion_x:
                    restrain = [1, 0]
                if restriccion_y is not None and nombre in restriccion_y:
                    restrain = [0, 1]
                if restriccion_xy is not None and nombre in restriccion_xy:
                    restrain = [1, 1]
                if node_id not in grupos[nombre]:
                    grupos[nombre][node_id] = Node(node_id + 1, [x, y], restrain=restrain)
                else:
                    grupos[nombre][node_id].restrain = restrain  # Actualiza si ya existe

    # Convertir a listas
    grupos_final = {nombre: list(nodos.values()) for nombre, nodos in grupos.items()}

    # Visualizar (si está disponible)
    #Node.plot_nodes_por_grupo(grupos_final, title, show_ids=False, save=False)

    return grupos_final, mesh

def make_sections(grupos, thickness_dict, E, nu, gamma):
    
    sections = {}

    for group in thickness_dict:
        material = Material(E, nu, gamma)
        sections[group] = Membrane(thickness_dict[group], material)

    nodes_dict = {}
    for group in grupos:
        for node in grupos[group]:
            nodes_dict[node.index] = node

    return sections, nodes_dict

def make_quad9_elements(mesh, sections, nodes_dict):
    quads = mesh.cells_dict.get('quad9', [])
    tags = mesh.cell_data_dict["gmsh:physical"].get("quad9", [])
    elements = []
    used_nodes = set()
    nodos_faltantes = []
    errores_jacobiano = []

    for i in range(len(tags)):
        section_tag = str(tags[i])
        if section_tag not in sections:
            print(f"⚠️ Tag físico {section_tag} no tiene sección asociada. Elemento {i + 1} omitido.")
            continue

        section = sections[section_tag]
        node_ids = quads[i]

        try:
            nodos = [nodes_dict[node_id + 1] for node_id in node_ids]
        except KeyError as e:
            nodos_faltantes.append(node_ids)
            print(f"❌ Nodo no encontrado en nodes_dict: {e}")
            continue

        for nodo in nodos:
            used_nodes.add(nodo)

        # Intentamos crear el elemento y capturamos errores de Jacobiano
        try:
            element = Quad9(i + 1, nodos, section)
            elements.append(element)
        except ValueError as ve:
            print(f"❌ Error en el elemento {i + 1} con Jacobiano no positivo:")
            print(f"   Nodos: {[n.index for n in nodos]}")
            print(f"   Coordenadas:")
            for j, n in enumerate(nodos):
                print(f"     Nodo local {j}: ID {n.index}, coord = {n.coord}")
            errores_jacobiano.append(i + 1)
            continue

    if nodos_faltantes:
        print(f"❌ Se omitieron {len(nodos_faltantes)} elementos por nodos faltantes.")
    if errores_jacobiano:
        print(f"⚠️ Se omitieron {len(errores_jacobiano)} elementos por Jacobiano negativo.")

    return elements, list(used_nodes)

def apply_distributed_force(grupo_nodos, fuerza_total_x, estructura):
    nodos = grupo_nodos
    n = len(nodos)
    if n < 2:
        print("Se requieren al menos dos nodos para aplicar fuerza distribuida.")
        return

    # Calcular posiciones acumuladas según distancia entre nodos (longitud sobre la curva)
    posiciones = [0.0]
    for i in range(1, n):
        dx = nodos[i].coord[0] - nodos[i-1].coord[0]
        dy = nodos[i].coord[1] - nodos[i-1].coord[1]
        distancia = np.sqrt(dx**2 + dy**2)
        posiciones.append(posiciones[-1] + distancia)
    total_longitud = posiciones[-1]

    # Inicializar fuerzas nodales
    nodal_forces = {}

    # Aplicar fuerza proporcional al tramo entre posiciones adyacentes
    for i in range(n):
        if i == 0:
            # Primer nodo: mitad de la diferencia con siguiente nodo
            fuerza = (posiciones[1] - posiciones[0]) / total_longitud * fuerza_total_x * 0.5
        elif i == n-1:
            # Último nodo: mitad de la diferencia con nodo anterior
            fuerza = (posiciones[-1] - posiciones[-2]) / total_longitud * fuerza_total_x * 0.5
        else:
            # Nodo interno: mitad de tramo anterior + mitad de tramo siguiente
            fuerza = ((posiciones[i] - posiciones[i-1]) + (posiciones[i+1] - posiciones[i])) / total_longitud * fuerza_total_x * 0.5
        nodal_forces[nodos[i].index] = fuerza

    # Aplicar fuerzas en X
    for node in nodos:
        fx = nodal_forces[node.index]
        dof_x, dof_y = node.dofs
        estructura.apply_force(dof_x, fx)
        estructura.apply_force(dof_y, 0.0)
        #print(f"Nodo {node.index} ← Fx = {fx:.3f} N, Fy = 0.000 N, coordenadas y = {node.coord[1]:.3f}")

def apply_self_weight(elements, rho, estructura):
    """
    Aplica peso propio a cada elemento Quad2D como fuerza puntual centrada interpolada.
    """
    g = 9.81
    P = 0.0

    for element in elements:
        centroid = element.get_centroid()
        area = element.A
        t = element.thickness
        peso = area * t * rho * g
        P += peso

        f_local = element.apply_point_body_force(
            x=centroid[0], y=centroid[1], force_vector=[0, -peso]
        )

        for idx_local, dof_global in enumerate(element.index):
            estructura.apply_force(dof_global, f_local[idx_local])

    print(f"✅ Peso total aplicado: {P:.3f} N")
    return P

def main(title, mesh_file, self_weight=True): 


    E = 210e3  # MPa
    nu = 0.3
    rho = 7800e-9 # kg/mm³

    thickness_dict = {"1": 200}

    restriccion_x = None #De otra forma hacer listas
    restriccion_y = None
    restriccion_xy = ["Restriccion XY"]

    grupos, mesh = make_nodes_groups(mesh_file, restriccion_x, restriccion_y, restriccion_xy)
    sections, nodes_dict = make_sections(grupos, thickness_dict, E, nu, rho)
    elements, used_nodes = make_quad9_elements(mesh, sections, nodes_dict)

    estructure = Solve(used_nodes, elements)

    if self_weight:
        
        # Aplicar peso propio a los elementos
        Peso = apply_self_weight(elements, rho, estructure)

    nodos_fuerza = grupos["Fuerza"]
    apply_distributed_force(nodos_fuerza, fuerza_total_x=-1200000, estructura=estructure)

 

    desplazamientos = estructure.solve()


    plot_results(
        estructure,
        elements,
        title=title,
        def_scale=1,
        force_scale=1e-4,
        reaction_scale=1e-2,
        sigma_y_tension=250, 
        sigma_y_compression=250
    )


if __name__ == "__main__":

    mesh_file = "GMSH_FILES/Quad9.msh"
    main(title="Test_Quad9/Resultados", mesh_file=mesh_file, self_weight=True)
