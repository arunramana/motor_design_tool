# %% ---------------------------------------------------------------------------
# Test tiziano basic functionalities
# ------------------------------------------------------------------------------
# %% Import tizinao from dolomites and other useful modules
import math
from dolomites import tiziano

# %% Init a tiziano drawing instance
draw = tiziano.drawing()

# %% Add points to the drawing
p1 = draw.add_point(1, 1)
p2 = draw.add_point(-1, 1)
p3 = draw.add_point(-1, -1)
p4 = draw.add_point(1, -1)

print('Points in the drawing: ------------------------------------------------')
print(draw.points)
draw.plot()  # visualize drawing with matplotlib

# %% Add lines to the drawing
l1 = draw.add_line(p1, p2, 10)[0]
l2 = draw.add_line(p2, p3, 10)[0]
l3 = draw.add_line(p3, p4, 10)[0]
l4 = draw.add_line(p4, p1, 10)[0]

draw.plot()  # visualize drawing with matplotlib

# %% Lines intersection
_ = draw.add_line(p1, p3)
_ = draw.add_line(p2, p4)

print('Lines in the drawing: -------------------------------------------------')
print(draw.lines)
draw.plot()  # visualize drawing with matplotlib

# %% Add discretized arcs to the drawing (with intersections)
p5 = draw.add_point(0.5, 0)
p6 = draw.add_point(-0.5, 0)
arc_pts1, arc_lns1 = draw.add_arc(p5, p6, math.pi, math.pi/18, 10)
arc_pts2, arc_lns2 = draw.add_arc(p5, p6, -math.pi, math.pi/18, 10)

# print actions for debugging
print('Discretized arcs in the drawing: --------------------------------------')
print('Points: =======================')
print(len(arc_pts1 + arc_pts2[1:-1]))
print('Lines: ========================')
print(len(arc_lns1 + arc_lns2))

draw.plot()  # visualize drawing with matplotlib

# %% Remove central point (and connected lines)
draw.remove_point(0.1)

draw.plot()  # visualize drawing with matplotlib

# %% Add labels and holes
lab1 = draw.add_label(0.7, 0, 101, 0.1)
lab2 = draw.add_label(0, 0.7, 102, 0.1)
lab3 = draw.add_label(-0.7, 0, 103, 0.1)
lab4 = draw.add_label(0, -0.7, 104, 0.1)
hol = draw.add_hole(0, 0)
draw.plot()  # visualize drawing with matplotlib

# %% Save drawing as tiziano file
draw.save('test_0.tiz')

# %% Call triangle to mesh the drawing
draw.mesh_triangle()

# %% Save mesh in gmsh format
draw.save_mesh('test_0.msh')
