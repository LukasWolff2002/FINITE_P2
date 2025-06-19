SetFactory("OpenCASCADE");

// Definición de puntos
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1 - Sqrt(2)/2, 0, 1};
Point(4) = {1 - Sqrt(2)/2, 1, 0, 1};
Point(5) = {0, 1, 0, 1};
Point(6) = {1, 1, 0, 1};  // Centro del círculo (1, 1)


// Crear líneas entre los puntos
Line(1) = {1, 2};
Line(2) = {2, 3};

Line(4) = {4, 5};  // Línea recta entre puntos 4 y 5
Line(5) = {5, 1};  // Línea recta entre puntos 5 y 1


// Crear un arco (círculo que conecta el punto 4 y 3 con centro en 6)
Circle(6) = {4, 6, 3};  // Círculo que conecta 4, 6 (centro), y 3

// Definir el número de divisiones en cada línea (usando transfinite en las líneas rectas)
Transfinite Curve {1, 2, 3, 4, 5} = 30 Using Progression 1.1;  // Dividir las líneas rectas en 30 segmentos
Transfinite Curve {6} = 30 Using Progression 1.1;  // Dividir el arco en 30 segmentos

// Crear un Loop cerrado y la superficie
Line Loop(1) = {1, 2, -6, 4, 5};  // Usamos las líneas para crear un loop cerrado, incluyendo el arco
Plane Surface(1) = {1};  // Definir la superficie cerrada (triangular)
//+
Transfinite Curve {1} = 5 ;
Transfinite Curve {2} = 5 ;
Transfinite Curve {4} = 2 ;
Transfinite Curve {5} = 6 ;
Transfinite Curve {6} = 4 ;
