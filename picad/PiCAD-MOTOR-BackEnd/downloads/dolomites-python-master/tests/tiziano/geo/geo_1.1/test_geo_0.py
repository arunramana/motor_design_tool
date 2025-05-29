# %% markdown
### Simple geometry to test tiziano geo module's functionality
# %% Import useful modules
import numpy as np
from dolomites import tiziano
import geo

# %% Init a new tiziano drawing instance
draw = tiziano.drawing()

# %% Draw a simple geometry (2 attached squares)
# points coordinates
x = [0, 2, 4, 4, 2, 0]
y = [0, 0, 0, 2, 2, 2]
# add points
pts = []
for ii, xx in enumerate(x):
    pts.append(draw.add_point(xx, y[ii], coherence=False))

# add lines
lns = []
for ii, pp in enumerate(pts[:-1]):
    lns.extend(draw.add_line(pp, pts[ii+1], ph=5, coherence=False))
lns.extend(draw.add_line(pts[-1], pts[0], ph=5, coherence=False))
lns.extend(draw.add_line(pts[1], pts[4], coherence=False))
lns.extend(draw.add_line(pts[2], pts[4], coherence=False))
lns.extend(draw.add_line(pts[0], pts[4], ph=8, coherence=False))

# %% Add a third separated square and modify one line ph tag
lns[2].ph_tag = 2

p1 = draw.add_point(5, 0, coherence=False)
p2 = draw.add_point(7, 0, coherence=False)
p3 = draw.add_point(7, 2, coherence=False)
p4 = draw.add_point(5, 2, coherence=False)

pts.extend([p1, p2, p3, p4])

l1 = draw.add_line(p1, p2, ph=4, coherence=False)
l2 = draw.add_line(p2, p3, coherence=False)
l3 = draw.add_line(p3, p4, ph=4, coherence=False)
l4 = draw.add_line(p1, p4, ph=3, coherence=False)

lns.extend(l1 + l2 + l3 + l4)

# %% Add labels
xl = [0.8, 1.2, 2.2, 3.2, 6]
yl = [1, 1, 1, 1, 1]

labs = []
for ii, xx in enumerate(xl):
    labs.append(draw.add_label(xx, yl[ii], (ii+1)*10, 2))
labs[1].ph_tag = 10
labs[3].ph_tag = 30
labs[4].ph_tag = 10

# %%
draw.plot()

# %% Call triangle mesher and save two different mesh models
# (1)
for ll in labs:  # decrease mesh trangles max are to see some difference
    ll.area = 0.5
draw.mesh_triangle()  # default flags='qpzaeA' inside dolomites/tiziano

# %%
# (2)
draw.mesh_triangle(flags='pzeA')

# %%
print('--> ', draw.triangle_points)
print('--> ', draw.triangle_edges)
print('--> ', draw.triangle_triangles)

# %%
# draw.save_mesh('example_mesh_1.msh')
# draw.save_mesh('example_mesh_2.msh')

# %% Save geo
geo.save_geo(draw, 'test_geo_0', 0.5, call_triangle=False)
