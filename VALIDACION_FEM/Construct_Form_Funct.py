import sympy as sp

# Definir variables simbólicas
xi, eta = sp.symbols('xi eta')

# Monomios simbólicos (base polinómica)
monomios = [sp.S(1), xi, eta, xi * eta]

# Coordenadas de los 7 nodos: 4 vértices y 3 centros de borde (ajusta si es necesario)
nodos = [
    (-1, -1),  # N1
    ( 1, -1),  # N2
    ( 1,  1),  # N3
    (-1,  1),  # N4
   

]

# Construir matriz del sistema A
A = sp.Matrix([[m.subs({xi: xj, eta: yj}) for m in monomios] for (xj, yj) in nodos])

# Resolver funciones de forma
Ns = []
for i in range(len(nodos)):
    b = sp.Matrix([1 if j == i else 0 for j in range(len(nodos))])
    coef = A.LUsolve(b)
    Ni = sum(c * m for c, m in zip(coef, monomios))
    Ni = sp.simplify(Ni)
    Ns.append(Ni)
    print(f"N{i+1}(xi, eta) = {Ni}\n")

# Validar suma
N_sum = sum(Ns)
N_sum_simplified = sp.simplify(N_sum)
N_sum_simplified = int(N_sum_simplified)  # Convertir a entero para evitar problemas de comparación
print("Suma de funciones de forma:", N_sum_simplified)
assert N_sum_simplified == 1, f"Error: La suma de las funciones de forma no es 1, es {N_sum_simplified}"
