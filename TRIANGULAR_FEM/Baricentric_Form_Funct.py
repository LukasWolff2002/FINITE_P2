import sympy as sp

def funciones_de_forma_triangular():
    # Definir coordenadas baricéntricas
    xi1, xi2, xi3 = sp.symbols('xi1 xi2 xi3')

    # Asegurarse de que la suma de las coordenadas baricéntricas sea 1
    suma_baricentrica = xi1 + xi2 + xi3 - 1
    assert sp.simplify(suma_baricentrica) == 0, f"Error: Las coordenadas no suman 1, es {suma_baricentrica}"

    # Definir las coordenadas físicas de los vértices (en este caso, un triángulo genérico)
    # Coordenadas de los vértices del triángulo
    x1, y1 = sp.symbols('x1 y1')
    x2, y2 = sp.symbols('x2 y2')
    x3, y3 = sp.symbols('x3 y3')

    # Relación entre coordenadas físicas y coordenadas baricéntricas
    x = x1 * xi1 + x2 * xi2 + x3 * xi3
    y = y1 * xi1 + y2 * xi2 + y3 * xi3

    # Imprimir las coordenadas físicas en función de las baricéntricas
    print(f"Coordenadas físicas: x = {sp.simplify(x)}, y = {sp.simplify(y)}")

    # Definir las funciones de forma en función de las coordenadas baricéntricas
    N1 = xi1
    N2 = xi2
    N3 = xi3

    # Imprimir las funciones de forma
    print(f"\nFunciones de forma:")
    print(f"N1(xi1, xi2, xi3) = {sp.simplify(N1)}")
    print(f"N2(xi1, xi2, xi3) = {sp.simplify(N2)}")
    print(f"N3(xi1, xi2, xi3) = {sp.simplify(N3)}")

    # Validar la suma de las funciones de forma
    N_sum = N1 + N2 + N3
    N_sum_simplified = sp.simplify(N_sum)
    print("\nSuma de funciones de forma:", N_sum_simplified)

    # Verificar que la suma de las funciones de forma es 1
    assert N_sum_simplified == 1, f"Error: La suma de las funciones de forma no es 1, es {N_sum_simplified}"

# Llamar a la función para un elemento triangular
funciones_de_forma_triangular()
