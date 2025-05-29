# %% Import useful modules
import numpy as np
from dolomites import tiziano

# %% Init a new tiziano drawing
draw = tiziano.drawing()

# %% Draw a simple trial drawing (a square)
# add points
p1 = draw.add_point(1, 1)
p2 = draw.add_point(-1, 1)
p3 = draw.add_point(-1, -1)
p4 = draw.add_point(1, -1)

# add lines
lns = []
lns.extend(draw.add_line(p1, p2))
lns.extend(draw.add_line(p2, p3))
lns.extend(draw.add_line(p3, p4))
lns.extend(draw.add_line(p4, p1))

# add label
lab = draw.add_label(0, 0, 100, 2)

# %% Plot picture
draw.plot()

# %% Call triangle mesher
mesh = draw.mesh_triangle()

# %% Print mesh_triangle output
print(draw.triangle_points)
print(draw.triangle_edges)
print(draw.triangle_triangles)

# %% Compute area of square region
area_sq = draw.area([100])
print('Square area computed with tiziano: ', area_sq)
print('Direct square ares: ', lns[0].manhattan_length()*lns[1].manhattan_length())

# %% markdown
# Computaton of areas of regions with only straight edges is very accurate.

# %% Draw another simple trial drawing (a circle surrounded by a square)
# add points
p01 = draw.add_point(0.5, 0)
p02 = draw.add_point(-0.5, 0)

# add arcs (discretized)
_, _ = draw.add_arc(p01, p02, np.pi, np.pi/30)
_, _ = draw.add_arc(p01, p02, -np.pi, np.pi/30)

# modify labels
draw.labels.clear()
lab0 = draw.add_label(0, 0, 1000, 0.05)
lab1 = draw.add_label(0, 0.8, 100, 0.2)

# %% Plot picture
draw.plot()

# %% Call triangle mesher
mesh = draw.mesh_triangle()

# %% Print mesh_triangle output
print(draw.triangle_points)
print(draw.triangle_edges)
print(draw.triangle_triangles)

# %% Compute area of circle region
area_circ = draw.area([1000])
print('Circle area computed with tiziano: ', area_circ)
print('Direct circle area: ', np.pi*0.5**2)

# %%
area_tot = draw.area([1000, 100])
print('Total area computed with tiziano: ', area_tot)
print('Direct total area: ', 2*2)

# %% markdown
# Accuracy of area computaton in regions with rounded edges depends heavily on  
# how arc edges are discretized.  
# A finer discretizaion leads to a more accurate result.
