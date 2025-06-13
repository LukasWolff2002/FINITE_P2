import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

class Solve:
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.elements = elements
        self.ndof = max(dof for node in nodes for dof in node.dofs)
        self.K_global = lil_matrix((self.ndof + 1, self.ndof + 1))
        self.f_global = lil_matrix((self.ndof + 1, 1))
        self.u_global = np.zeros((self.ndof + 1, 1))  # sigue siendo denso

    def assemble(self):
        for elem in self.elements:
            ke = elem.Kg
            idx = elem.calculate_indices()

            if ke.shape != (len(idx), len(idx)):
                raise ValueError(f"❌ Dimensión inconsistente: Ke {ke.shape}, idx {len(idx)}")

            for i in range(len(idx)):
                for j in range(len(idx)):
                    self.K_global[idx[i], idx[j]] += ke[i, j]

    def apply_force(self, dof_index, value):
        self.f_global[dof_index, 0] += value

    def apply_forces_vector(self, force_vector):
        for i, val in enumerate(force_vector):
            self.f_global[i, 0] += val

    def apply_boundary_conditions(self):
        self.fixed_dofs = []
        self.free_dofs = []

        for node in self.nodes:
            for dof_val, dof_idx in zip(node.restrain, node.dofs):
                if dof_val == 1:
                    self.fixed_dofs.append(dof_idx)
                else:
                    self.free_dofs.append(dof_idx)

        self.fixed_dofs = np.array(self.fixed_dofs)
        self.free_dofs = np.array(self.free_dofs)

        for dof in self.fixed_dofs:
            self.K_global[dof, :] = 0
            self.K_global[:, dof] = 0
            self.K_global[dof, dof] = 1
            self.f_global[dof, 0] = 0

    def check_zero_rows(self):
        return np.where(~self.K_global.toarray().any(axis=1))[0]

    def solve(self):
        self.assemble()

        # Guardar copias dispersas
        self.K_original = self.K_global.copy().tocsr()
        self.f_original = self.f_global.copy().tocsc()

        self.apply_boundary_conditions()

        used_dofs = sorted(set(dof for node in self.nodes for dof in node.dofs))

        K_reduced = self.K_global[used_dofs, :][:, used_dofs].tocsr()
        f_reduced = self.f_global[used_dofs].toarray()

        rowsums = np.sum(np.abs(K_reduced.toarray()), axis=1)
        zero_rows = np.where(rowsums == 0)[0]
        if len(zero_rows) > 0:
            raise ValueError("Sistema subdeterminado: nodos con DOF sin rigidez.")

        u_reduced = spsolve(K_reduced, f_reduced).reshape(-1, 1)

        self.u_global = np.zeros_like(self.u_global)
        self.u_global[used_dofs] = u_reduced

        return self.u_global

    def get_displacement_at_node(self, node_index):
        node = self.nodes[node_index - 1]
        ux = self.u_global[node.dofs[0], 0]
        uy = self.u_global[node.dofs[1], 0]
        return ux, uy

    def compute_reactions(self):
        R_total = self.K_original @ self.u_global - self.f_original.toarray()

        self.reactions = np.zeros_like(self.u_global)
        self.reactions[self.fixed_dofs] = R_total[self.fixed_dofs]

        return self.reactions

    def print_applied_forces(self):
        f = self.f_original.toarray() if hasattr(self, 'f_original') else self.f_global.toarray()
        for node in self.nodes:
            dof_x, dof_y = node.dofs
            fx = f[dof_x][0] if dof_x < len(f) else 0.0
            fy = f[dof_y][0] if dof_y < len(f) else 0.0
            print(f"Nodo {node.index}: Fx = {fx:.2f}, Fy = {fy:.2f}")

    def print_summary(self):
        for node in self.nodes:
            ux, uy = self.get_displacement_at_node(node.index)
            print(f"Nodo {node.index}: ux = {ux:.6e}, uy = {uy:.6e}")
