# %% markdown
### Trial geometry to test tiziano geo module's capabilities

# %% Import useful modules
import numpy as np
from dolomites import tiziano
import geo

# %% Init a new tiziano drawing instance
draw = tiziano.drawing()

# %% Draw geometry in current tiziano instance
x1 = [-3,  3, 3, -3]
y1 = [-3, -3, 3,  3]

x2 = [0, 0, -1, -1]
y2 = [0, 1,  1,  0]

pts = []
lns = []

# external square
for ii, xx in enumerate(x1):
    pts.append(draw.add_point(xx, y1[ii], coherence=False))
for ii, pt in enumerate(pts[:-1]):
    lns.extend(draw.add_line(pt, pts[ii+1], ph=2, coherence=False))
lns.extend(draw.add_line(pts[-1], pts[0], ph=2, coherence=False))
draw.add_label(0, -2, 8, 0.5);

# first internal square
for ii, xx in enumerate(x2):
    pts.append(draw.add_point(xx, y2[ii], coherence=False))
for ii, pt in enumerate(pts[4:-1], 4):
    lns.extend(draw.add_line(pt, pts[ii+1], ph=4, coherence=False))
lns.extend(draw.add_line(pts[-1], pts[4], ph=4, coherence=False))
hol = draw.add_hole(-0.5, 0.5)

# second internal square
draw.mirror(lns[4:]+[hol], 0, 0, 1, 1, coherence=False);

# %%
draw.plot()

# %% Save geo
# draw.mesh_triangle()
geo.save_geo(draw, 'test_geo_2', 0.2, call_triangle=True, gmsh_mesh=None)
    
