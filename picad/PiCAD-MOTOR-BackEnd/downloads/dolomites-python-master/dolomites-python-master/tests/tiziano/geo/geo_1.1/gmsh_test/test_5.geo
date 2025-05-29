// Simple geometry to test Gmsh mesher capability ------------------------------

// Question: observe the behaviour of the mesher with another twisted curve loop?
// -----------------------------------------------------------------------------

lc = 0.5;
// -----------------------------------------------------------------------------
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Point(5) = { 2,  1, 0, lc};
Point(6) = { 2, -1, 0, lc};
Point(7) = {-1, -1, 0, lc};
Point(8) = {-1,  2, 0, lc};
Point(9) = { 1,  2, 0, lc};

// -----------------------------------------------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {3, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 3};

// -----------------------------------------------------------------------------
// from inside to outside
Curve Loop(1) = {1, 2, -10, -9, -8, -7, -6, -5, 3, 4};
// from outside to inside
// Curve Loop(1) = {8, 9, 10, -2, -1, -4, -3, 5, 6, 7};
// -----------------------------------------------------------------------------
Plane Surface(1) = {1};
// -----------------------------------------------------------------------------
Physical Surface(3) = {1};

// -----------------------------------------------------------------------------

Mesh 2;  // do 2D mesh
Save "test_5.msh";
// -----------------------------------------------------------------------------
// --> Answer: the mesh is generated as wanted (i.e. the hole is recognized).
//             However, when changing the defintion of the twisted curve loop,
//             the orthonormal orientation of the 2D mesh elements seems to change.
