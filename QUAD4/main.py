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
from Quad2D import Quad2D
from solve import Solve
from graph import plot_results

def make_nodes_groups(output_file, restriccion_x = None, restriccion_y = None, restriccion_xy = None):
    
    mesh = meshio.read(output_file)
    
    # Traducimos los tags físicos a nombres
    tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}
    grupos = {}


    # Elementos tipo "quad" → para el dominio estructural
    for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
        if cell_block.type != "quad":
            continue
        
        for quad, tag in zip(cell_block.data, phys_tags):
            nombre = tag_to_name.get(tag, f"{tag}")

            if nombre not in grupos:
                grupos[nombre] = []
              
            for node_id in quad:
                x, y = mesh.points[node_id][:2]
                grupos[nombre].append(Node(node_id + 1, [x, y]))

    # Elementos tipo "line" → condiciones de borde
    for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
        if cell_block.type != "line":
            continue
        for line, tag in zip(cell_block.data, phys_tags):
            nombre = tag_to_name.get(tag, f"{tag}")
            if nombre not in grupos:
                grupos[nombre] = []
            for node_id in line:
                x, y = mesh.points[node_id][:2]
                restrain = [0, 0]
                if restriccion_x is not None and nombre in restriccion_x:
                    restrain = [1, 0]
                if restriccion_y is not None and nombre in restriccion_y:
                    restrain = [0, 1]
                if restriccion_xy is not None and nombre in restriccion_xy:
                    restrain = [1, 1]
                grupos[nombre].append(Node(node_id + 1, [x, y], restrain=restrain))

    # Eliminar nodos duplicados por grupo (según id)
    for nombre in grupos:
        nodos_unicos = {}
        for nodo in grupos[nombre]:
            nodos_unicos[nodo.index] = nodo
        grupos[nombre] = list(nodos_unicos.values())

    # Visualización opcional
    #Node.plot_nodes_por_grupo(grupos, title, show_ids=False, save=False)
    print(grupos)
    return grupos, mesh

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

def make_quad2d_elements(mesh, sections, nodes_dict):
    quads = mesh.cells_dict.get('quad', [])
    tags = mesh.cell_data_dict["gmsh:physical"].get("quad", [])
    elements = []
    used_nodes = set()
    nodos_faltantes = []

    for i in range(len(tags)):
        section_tag = str(tags[i])
        if section_tag not in sections:
            print(f"⚠️ Tag físico {section_tag} no tiene sección asociada. Elemento {i+1} omitido.")
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

        #print(nodos)

        element = Quad2D(i + 1, nodos, section)
        elements.append(element)

    if nodos_faltantes:
        print(f"❌ Se omitieron {len(nodos_faltantes)} elementos por nodos faltantes.")
    
    return elements, list(used_nodes)

def plot_all_elements(elements, title, show_ids=True):
    all_x = []
    all_y = []

    # Recopilar coordenadas de todos los nodos
    for elem in elements:
        coords = elem.xy  # accede directamente a las coordenadas
        coords = np.vstack([coords, coords[0]])  # cerrar el polígono
        all_x.extend(coords[:, 0])
        all_y.extend(coords[:, 1])

    # Márgenes y límites
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 8
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))

    for elem in elements:
        coords = elem.xy
        coords = np.vstack([coords, coords[0]])  # cerrar el polígono

        ax.plot(coords[:, 0], coords[:, 1], 'k-', linewidth=1)

        if show_ids:
            for nodo, (x, y) in zip(elem.node_list, coords[:-1]):
                ax.text(x, y, f'N{nodo.index}', color='black', fontsize=6, ha='center', va='center')

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("All Quad2D elements")
    ax.grid(True)

    plt.show()
    
def apply_distributed_force_x(grupo_nodos, fuerza_total_x, estructura):
    """
    Aplica una fuerza distribuida vertical (por ejemplo, peso) sobre una línea formada por nodos.
    La fuerza se reparte proporcionalmente a la longitud de los tramos y se descompone en x e y.
    """

    nodos = grupo_nodos
    n = len(nodos)
    if n < 2:
        print("Se requieren al menos dos nodos para aplicar fuerza distribuida.")
        return

    # Paso 1: calcular longitud total
    longitudes = []
    total_length = 0
    for i in range(n - 1):
        dx = nodos[i+1].coord[0] - nodos[i].coord[0]
        dy = nodos[i+1].coord[1] - nodos[i].coord[1]
        L = np.sqrt(dx**2 + dy**2)
        longitudes.append(L)
        total_length += L

    q_lineal = fuerza_total_x / total_length  # fuerza por metro

    # Paso 2: inicializar diccionario de fuerzas
    nodal_forces = {node.index: np.array([0.0, 0.0]) for node in nodos}

    for i in range(n - 1):
        ni = nodos[i]
        nj = nodos[i + 1]
        xi, yi = ni.coord
        xj, yj = nj.coord

        dx = xj - xi
        dy = yj - yi
        L = longitudes[i]

        # Vector perpendicular hacia "abajo"
        
        vy = dy / L
        nx = -vy
      

        Fi = q_lineal * L
        fx = Fi * nx
       
        nodal_forces[ni.index] += np.array([fx / 2, 0])
        nodal_forces[nj.index] += np.array([fx / 2, 0])

    # Paso 3: aplicar fuerzas al sistema
    for node in nodos:
        fx, fy = nodal_forces[node.index]
        dof_x, dof_y = node.dofs
        #estructura.apply_force(dof_x, fx)
        estructura.apply_force(dof_x, fx)
        #print(f"Nodo {node.index} ← Fx = {fx:.3f} N, Fy = {fy:.3f} N")

def apply_distributed_force_y(grupo_nodos, fuerza_total_y, estructura):
    """
    Aplica una fuerza distribuida vertical (por ejemplo, peso) sobre una línea formada por nodos.
    La fuerza se reparte proporcionalmente a la longitud de los tramos y se descompone en x e y.
    """

    nodos = grupo_nodos
    n = len(nodos)
    if n < 2:
        print("Se requieren al menos dos nodos para aplicar fuerza distribuida.")
        return

    # Paso 1: calcular longitud total
    longitudes = []
    total_length = 0
    for i in range(n - 1):
        dx = nodos[i+1].coord[0] - nodos[i].coord[0]
        dy = nodos[i+1].coord[1] - nodos[i].coord[1]
        L = np.sqrt(dx**2 + dy**2)
        longitudes.append(L)
        total_length += L

    q_lineal = fuerza_total_y / total_length  # fuerza por metro

    # Paso 2: inicializar diccionario de fuerzas
    nodal_forces = {node.index: np.array([0.0, 0.0]) for node in nodos}

    for i in range(n - 1):
        ni = nodos[i]
        nj = nodos[i + 1]
        xi, yi = ni.coord
        xj, yj = nj.coord

        dx = xj - xi
        dy = yj - yi
        L = longitudes[i]

        # Vector perpendicular hacia "abajo"
        
        vx = dy / L
        ny = -vx
      

        Fi = q_lineal * L
        fy = Fi * ny
       
        nodal_forces[ni.index] += np.array([0, fy / 2])
        nodal_forces[nj.index] += np.array([0, fy / 2])

    # Paso 3: aplicar fuerzas al sistema
    for node in nodos:
        fx, fy = nodal_forces[node.index]
        dof_x, dof_y = node.dofs
        #estructura.apply_force(dof_x, fx)
        estructura.apply_force(dof_y, fy)
        #print(f"Nodo {node.index} ← Fx = {fx:.3f} N, Fy = {fy:.3f} N")


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
    sections, nodes_dict = make_sections(grupos, thickness_dict=thickness_dict,E=E, nu=nu, gamma=rho)
    elements, used_nodes = make_quad2d_elements(mesh, sections, nodes_dict)

    #plot_all_elements(elements, "All Quad2D elements", show_ids=True)

    estructure = Solve(used_nodes, elements)

    if self_weight:
        
        # Aplicar peso propio a los elementos
        apply_self_weight(elements, rho, estructure)

        pass

    nodos_fuerza = grupos["Fuerza"]
    apply_distributed_force_y(nodos_fuerza, fuerza_total_y=1200000, estructura=estructure)

    estructure.solve()

    plot_results(
        estructure,
        elements,
        title=title,
        def_scale=1,
        force_scale=1,
        reaction_scale=1e-2,
        sigma_y_tension=250, 
        sigma_y_compression=250
    )
    

if __name__ == "__main__":

    mesh_file = "GMSH_FILES/Quad4.msh"
    main(title="Test_Quad4/Resultados", mesh_file=mesh_file, self_weight=True)
