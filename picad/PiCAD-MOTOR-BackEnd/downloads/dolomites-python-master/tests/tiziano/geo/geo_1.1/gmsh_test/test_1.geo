// Simple geometry to test Gmsh mesher capability ------------------------------

// Question: can Gmsh mesher do the the mesh of a plane surface with a hole
// defined using a twisted curve loop?
// -----------------------------------------------------------------------------

lc = 0.2;
// -----------------------------------------------------------------------------
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Point(5) = {2, 1, 0, lc};
Point(6) = {2, 2, 0, lc};
Point(7) = {1, 2, 0, lc};

Point(8) = {-1, -1, 0, lc};
Point(9) = {3, -1, 0, lc};
Point(10) = {3, 3, 0, lc};
Point(11) = {-1, 3, 0, lc};
// -----------------------------------------------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {3, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 3};

Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 8};
// -----------------------------------------------------------------------------
Curve Loop(1) = {1, 2, -8, -7, -6, -5, 3, 4};
// Curve Loop(1) = -1*{1, 2, -8, -7, -6, -5, 3, 4};
// Curve Loop(2) = {-9, -10, -11, -12};
Curve Loop(2) = -1*{-9, -10, -11, -12};
// -----------------------------------------------------------------------------
Plane Surface(1) = {2, 1};
// -----------------------------------------------------------------------------
Physical Surface(3) = {1};

Physical Line(2) = {1, 2, 3, 4, 5, 6, 7, 8};
// -----------------------------------------------------------------------------
Mesh 2;  // 2D mesh
Save "test_1.msh";
// -----------------------------------------------------------------------------
// --> Answer: yes, when the twisted loop describes a hole, it can.
//             The surface is meshed as wanted.
//             However, a special care should be taken in deciding the curve loops
//             direction: the lines/points sorting along the loop determines the
//             sign of the normals of the 2D mesh elements inside the surface.
//             (+ <==> counterclockwise) ; (- <==> clockwise)
//             Negative normals may cause troubles in FE formulations.          
