SetFactory("OpenCASCADE");

// Geometría
Point(1) = {0,0,0,1};
Point(2) = {2000,0,0,1};
Point(3) = {2000,200,0,1};
Point(4) = {0,200,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Mallado estructurado
Transfinite Curve{1, 3} = 200 Using Progression 1;
Transfinite Curve{2, 4} = 40 Using Progression 1;
Transfinite Surface{1};
Recombine Surface{1};  // ← importante para tener quads

// Etiquetas físicas para FEM
Physical Surface(1) = {1};
Physical Line("Restriccion XY") = {4};
Physical Line("Fuerza") = {2};

// Malla automática
Mesh 2;
