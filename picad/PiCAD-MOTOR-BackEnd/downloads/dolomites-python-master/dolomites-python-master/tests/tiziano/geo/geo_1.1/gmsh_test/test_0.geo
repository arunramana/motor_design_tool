// Simple geometry to test Gmsh mesher capability ------------------------------

// Question: can Gmsh mesher do the the mesh of a plane surface with a
// physical line inside?
// -----------------------------------------------------------------------------

lc = 0.1;
// -----------------------------------------------------------------------------
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
// -----------------------------------------------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {2, 4};
// -----------------------------------------------------------------------------
Curve Loop(1) = {1, 2, 3, 4};
// -----------------------------------------------------------------------------
Plane Surface(1) = {1};
// -----------------------------------------------------------------------------
Physical Surface(3) = {1};

Physical Line(2) = {1, 2, 3, 4};
Physical Line(5) = {5};
// -----------------------------------------------------------------------------
Mesh 2;  // 2D mesh
Save "test_0.msh";
// -----------------------------------------------------------------------------
// --> Answer: yes it can, but the nodes along the line do not coindice with 2D
//             mesh elements node.
// ==> Is this an issue for some types of FE formulations?
