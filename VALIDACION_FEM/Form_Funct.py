import numpy as np
import sympy as sp

# Definimos las variables simbólicas
xi, eta = sp.symbols('xi eta')

# Funciones de forma para un elemento cuadrilátero de 9 nodos (Quad9)
def N1(xi, eta):
    return 0.25 * xi * eta * (xi - 1) * (eta - 1)

def N2(xi, eta):
    return 0.25 * xi * eta * (xi + 1) * (eta - 1)

def N3(xi, eta):
    return 0.25 * xi * eta * (xi + 1) * (eta + 1)

def N4(xi, eta):
    return 0.25 * xi * eta * (xi - 1) * (eta + 1)

def N5(xi, eta):
    return 0.5 * (1 - xi ** 2) * eta * (eta - 1)

def N6(xi, eta):
    return 0.5 * xi * (xi + 1) * (1 - eta ** 2)

def N7(xi, eta):
    return 0.5 * (1 - xi ** 2) * eta * (eta + 1)

def N8(xi, eta):
    return 0.5 * xi * (xi - 1) * (1 - eta ** 2)

def N9(xi, eta):
    return (1 - xi ** 2) * (1 - eta ** 2)

# Función para calcular derivadas algebraicas
def calculate_derivatives():
    N_funcs = [N1(xi, eta), N2(xi, eta), N3(xi, eta), N4(xi, eta),
               N5(xi, eta), N6(xi, eta), N7(xi, eta), N8(xi, eta), N9(xi, eta)]
    
    # Derivadas respecto a xi (zeta)
    dN_dxi = [sp.diff(N, xi) for N in N_funcs]
    
    # Derivadas respecto a eta
    dN_deta = [sp.diff(N, eta) for N in N_funcs]

    print("Derivadas respecto a xi (zeta):")
    for i, dN in enumerate(dN_dxi):
        print(f"dN{i + 1}/dξ = {dN}")
    
    print("\nDerivadas respecto a eta:")
    for i, dN in enumerate(dN_deta):
        print(f"dN{i + 1}/dη = {dN}")

# Coordenadas de los nodos en el plano (xi, eta)
node_coords = {
    1: (-1, -1),
    2: ( 1, -1),
    3: ( 1,  1),
    4: (-1,  1),
    5: ( 0, -1),
    6: ( 1,  0),
    7: ( 0,  1),
    8: (-1,  0),
    9: ( 0,  0)
}

# Función para validar las funciones de forma para un nodo específico
def validate_shape_function(N_funcs, node_index):
    xi, eta = node_coords[node_index + 1]  # Obtener las coordenadas del nodo
    
    # Validar la función de forma del nodo específico
    value = N_funcs[node_index](xi, eta)
    print(f"En el nodo {node_index + 1}: N{node_index + 1} = {value}")
    
    # Validar que la función de forma en su nodo correspondiente sea 1
    if node_index + 1 == node_index + 1:
        assert np.isclose(value, 1), f"Error: N{node_index + 1} en nodo {node_index + 1} no es 1."
    else:
        # Validar que la función de forma en los otros nodos sea 0
        assert np.isclose(value, 0), f"Error: N{node_index + 1} en nodo {node_index + 1} no es 0."

    # Validar que las demás funciones sean 0 en el nodo actual
    for i in range(9):
        if i != node_index:  # Si no es el mismo nodo
            value_other = N_funcs[i](xi, eta)
            print(f"En el nodo {node_index + 1}: N{i + 1} = {value_other}")
            assert np.isclose(value_other, 0), f"Error: N{i + 1} en nodo {node_index + 1} no es 0."

# Función para validar todas las funciones de forma para cada nodo
def validate_all_shape_functions():
    N_funcs = [N1, N2, N3, N4, N5, N6, N7, N8, N9]
    for node_index in range(9):
        print(f"\n--- Validando la función de forma para el nodo {node_index + 1} ---")
        validate_shape_function(N_funcs, node_index)

# Llamar a la función para validar las funciones de forma para todos los nodos
validate_all_shape_functions()

# Llamar a la función para calcular las derivadas simbólicas de las funciones de forma
calculate_derivatives()
