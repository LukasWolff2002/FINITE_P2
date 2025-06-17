SetFactory("OpenCASCADE");

// Puntos del cuadrado
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

// Líneas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Superficie
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Número de divisiones en las curvas (para refinamiento)
n = 10;
r = 1.1;

// Refinamiento hacia (0,0) con progresión de la malla
Transfinite Curve{1, 2} = n Using Progression r;  // Curva entre (1,2)
Transfinite Curve{4, 3} = n Using Progression 1/r;  // Curva entre (4,3)

// Refinamiento de la superficie del cuadrado
Transfinite Surface {1} = {2, 3, 4, 1};  // Cuadrado refinado

// Asignación de superficies físicas
Physical Surface("Dominio") = {1};

// Asignación de condiciones de contorno (Dirichlet) para cada borde
Physical Line("Dirichlet 1") = {1};  // Lado 1 (de 1 a 2)
Physical Line("Dirichlet 2") = {2};  // Lado 2 (de 2 a 3)
Physical Line("Dirichlet 3") = {3};  // Lado 3 (de 3 a 4)
Physical Line("Dirichlet 4") = {4};  // Lado 4 (de 4 a 1)
