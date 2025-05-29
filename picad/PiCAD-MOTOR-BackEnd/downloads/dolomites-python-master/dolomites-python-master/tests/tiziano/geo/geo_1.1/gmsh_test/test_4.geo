// Simple geometry to test Gmsh mesher capability ------------------------------

// Question: How nodes acquire physical tags?
// -----------------------------------------------------------------------------

lc = 0.5;
// -----------------------------------------------------------------------------
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Point(5) = {2, 1, 0, lc};
Point(6) = {2, 2, 0, lc};
Point(7) = {1, 2, 0, lc};

// -----------------------------------------------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {3, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 3};

// -----------------------------------------------------------------------------
Curve Loop(1) = {1, 2, 3, 4};
// -----------------------------------------------------------------------------
Plane Surface(1) = {1};
// -----------------------------------------------------------------------------
Physical Surface(3) = {1};

Physical Line(2) = {1, 2, 3, 4};
Physical Line(7) = {5, 6, 7, 8};
// -----------------------------------------------------------------------------
// The default msh format version is the latest stable version.
// Use the command below to change msh format.\

Mesh.MshFileVersion = 2.2;  // <-- choose .msh format version
// Mesh.MshFileVersion = 4.1;  // <-- choose .msh format version (default 4.1, 18-6-2022)
Mesh 2;  // do 2D mesh
// Save "test_4_v22.msh";  // version 2.2
// Save "test_4_v41.msh";  // version 4.1
// -----------------------------------------------------------------------------
// --> Answer: nodes ph tags are not specified in any .msh file format (until 4.1).

//             Looking at the geometry and mesh in Gmsh GUI, we can see that
//             in format version 2.2 (test_4_v22.msh) every node acquires the
//             ph tag of one of the mesh element (2D or 1D) which is part of.

//             The same thing happens in msh format version 4.1 (test_4_v41.msh),
//             except for those nodes that coincide with points defined in the
//             geo model. All these nodes have a 0 ph tag.

//             Note that when using the command <Mesh 2> in a geo script the parser
//             generates the mesh with the last stable version (now 18-6-2022, v. 4.1).
//             The option Mesh.MshFileVersion changes the format only of the output
//             msh file, i.e. the option acts only when saving the actual mesh in a file.
//             You can note this behavior by openin this geo file with command
//             <Mesh 2D> uncommented and Mesh.MshFileVersion set to 2.2.
//             All nodes that coincide with geo points have still a 0 ph tag.  
