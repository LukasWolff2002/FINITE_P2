import numpy as np
import sympy as sp

class Element:
    def __init__(self, id, node_ids, n_nodes):
        self.id = id
        self.node_ids = node_ids  # Lista de nodos
        self.num_nodes = n_nodes  # Número de nodos en el elemento
        self.K = None  # Se calculará con get_stiffness_matrix

    def shape_functions(self, xi, eta):
        # Funciones de forma para un elemento cuadrilátero de 9 nodos
        N1 = 0.25 * xi * eta * (xi - 1) * (eta - 1)
        N2 = 0.25 * xi * eta * (xi + 1) * (eta - 1)
        N3 = 0.25 * xi * eta * (xi + 1) * (eta + 1)
        N4 = 0.25 * xi * eta * (xi - 1) * (eta + 1)
        N5 = 0.5 * (1 - xi ** 2) * eta * (eta - 1)
        N6 = 0.5 * xi * (xi + 1) * (1 - eta ** 2)
        N7 = 0.5 * (1 - xi ** 2) * eta * (eta + 1)
        N8 = 0.5 * xi * (xi - 1) * (1 - eta ** 2)
        N9 = (1 - xi ** 2) * (1 - eta ** 2)

        # Generar la lista de funciones de forma
        N = np.array([N1, N2, N3, N4, N5, N6, N7, N8, N9])
        
        return N

    def shape_function_derivatives(self, xi, eta):
        # Derivadas de las funciones de forma respecto a xi y eta
        N = self.shape_functions(xi, eta)

        # Las derivadas simbólicas con sympy
        xi, eta = sp.symbols('xi eta')
        N_sym = np.array([
            0.25 * xi * eta * (xi - 1) * (eta - 1),
            0.25 * xi * eta * (xi + 1) * (eta - 1),
            0.25 * xi * eta * (xi + 1) * (eta + 1),
            0.25 * xi * eta * (xi - 1) * (eta + 1),
            0.5 * (1 - xi ** 2) * eta * (eta - 1),
            0.5 * xi * (xi + 1) * (1 - eta ** 2),
            0.5 * (1 - xi ** 2) * eta * (eta + 1),
            0.5 * xi * (xi - 1) * (1 - eta ** 2),
            (1 - xi ** 2) * (1 - eta ** 2)
        ])

        dN_dxi = np.array([sp.diff(N_i, xi) for N_i in N_sym])
        dN_deta = np.array([sp.diff(N_i, eta) for N_i in N_sym])

        # Convertir de sympy a numpy para la evaluación numérica
        dN_dxi_func = np.array([sp.lambdify((xi, eta), d, "numpy") for d in dN_dxi])
        dN_deta_func = np.array([sp.lambdify((xi, eta), d, "numpy") for d in dN_deta])

        return dN_dxi_func, dN_deta_func

    def get_B_matrix(self, nodes, xi, eta):
        coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in self.node_ids])  # Coordenadas de nodos
      
        dN_dxi_func, dN_deta_func = self.shape_function_derivatives(xi, eta)

        J = np.zeros((2, 2))
        for i in range(self.num_nodes):  
            J[0, 0] += dN_dxi_func[i](xi, eta) * coords[i, 0]
            J[0, 1] += dN_dxi_func[i](xi, eta) * coords[i, 1]
            J[1, 0] += dN_deta_func[i](xi, eta) * coords[i, 0]
            J[1, 1] += dN_deta_func[i](xi, eta) * coords[i, 1]

        detJ = np.linalg.det(J)
        if abs(detJ) < 1e-12:
            return np.zeros((2, self.num_nodes)), 0.0

        Jinv = np.linalg.inv(J)

        dN_dx = Jinv[0, 0] * np.array([dN_dxi_func[i](xi, eta) for i in range(self.num_nodes)]) + \
                Jinv[0, 1] * np.array([dN_deta_func[i](xi, eta) for i in range(self.num_nodes)])
        dN_dy = Jinv[1, 0] * np.array([dN_dxi_func[i](xi, eta) for i in range(self.num_nodes)]) + \
                Jinv[1, 1] * np.array([dN_deta_func[i](xi, eta) for i in range(self.num_nodes)])

        B = np.vstack((dN_dx, dN_dy))

        return B, detJ

    def get_stiffness_matrix(self, nodes):
        # Puntos de Gauss
        gauss_points = [
            (1/6, 1/6, 1/6),
            (2/3, 1/6, 1/6),
            (1/6, 2/3, 1/6)
        ]

        K = np.zeros((self.num_nodes, self.num_nodes))  # Matriz de rigidez general

        for xi, eta, w in gauss_points:
            B, detJ = self.get_B_matrix(nodes, xi, eta)
            if detJ == 0:
                continue
            K += w * detJ * (B.T @ B)

        self.K = K
        return self.K

    def get_load_vector(self, nodes, alpha):
        coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in self.node_ids])

        gauss_points = [
            (1/6, 1/6, 1/6),
            (2/3, 1/6, 1/6),
            (1/6, 2/3, 1/6)
        ]

        f_local = np.zeros(self.num_nodes)  # Vector de carga local

        for xi, eta, w in gauss_points:
            N = self.shape_functions(xi, eta)
            x = np.dot(N, coords[:, 0])
            y = np.dot(N, coords[:, 1])
            r2 = x**2 + y**2

            f_val = 0.0
            if not (r2 == 0 and alpha < 1):
                f_val = -alpha**2 * r2**(alpha / 2 - 1)

            _, detJ = self.get_B_matrix(nodes, xi, eta)
            f_local += f_val * w * detJ * N

        return f_local
