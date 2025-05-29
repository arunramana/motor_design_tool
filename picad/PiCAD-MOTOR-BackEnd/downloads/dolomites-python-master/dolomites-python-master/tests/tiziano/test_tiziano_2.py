# %% ---------------------------------------------------------------------------
# Test tiziano mesh and save mesh capabilities
# ------------------------------------------------------------------------------
# %% Import tizinao from dolomites and other useful modules
from dolomites import tiziano

# %% Init a tiziano drawing instance
draw = tiziano.drawing()

# %% Load test_0 into current drawing
draw.load('test_0.tiz')
draw.plot()

# %% Simplify drawing
draw.remove_point(0.6)
draw.labels.clear()
draw.holes.clear()

# %% Add other points and lines to the drawing
p0 = draw.add_point(-1, -1)
p1 = draw.add_point(-1, 2)
p2 = draw.add_point(2, 2)
p3 = draw.add_point(2, -1)

l1 = draw.add_line(p0, p1, 20)[0]
l2 = draw.add_line(p1, p2, 20)[0]
l3 = draw.add_line(p2, p3, 20)[0]
l4 = draw.add_line(p3, p0, 20)[0]

# %% Add labels and holes
lab = draw.add_label(0, 0, 100, 1)
hol = draw.add_hole(1.5, 1.5)

# %% Visulize drawing with matplotlib
draw.plot()

# %% Call triangle mesher
mesh = draw.mesh_triangle()

# %% print mesh entities
print('Triangle mesh nodes ---------------------------------------------------')
print(len(draw.triangle_points))
print(draw.triangle_points)
print('Triangle mesh edges ---------------------------------------------------')
print(len(draw.triangle_edges))
print(draw.triangle_edges)
print('Triangle mesh elements ------------------------------------------------')
print(len(draw.triangle_triangles))
print(draw.triangle_triangles)

# %% Add other nodes and edges to tiziano drawing mesh lists

draw.triangle_edges.extend([((1, 4), 20), ((4, 5), 20), ((5, 6), 20), ((6, 3), 20)])

# %% Save mesh in gmsh format
draw.save_mesh('test_2.msh')
