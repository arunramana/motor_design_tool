# %% markdown
### Trial geometry to test tiziano geo module's capabilities

# %% Import useful modules
import numpy as np
from dolomites import tiziano
import geo

# %% Init a new tiziano drawing instance
draw = tiziano.drawing()

# %% Draw geometry in current tiziano instance
x = [ 0, 1, 0,  0, 8, 0]
y = [-1, 0, 1, -6, 0, 6]

pts = []
for ii, xx in enumerate(x):
    pts.append(draw.add_point(xx, y[ii], coherence=False))

lns = []
lns.extend(draw.add_line(pts[0], pts[1], coherence=False))
lns.extend(draw.add_line(pts[1], pts[2], coherence=False))
lns.extend(draw.add_line(pts[3], pts[4], ph=7, coherence=False))
lns.extend(draw.add_line(pts[4], pts[5], ph=7, coherence=False))

_, _, _, _ = draw.mirror(lns, 0, 0, 0, 1, coherence=False)

to_move = draw.select_rect(-1.1, -1.1, 2.2, 2.2, entity=1)

draw.copy_traslate(to_move, 2, 0, N=2, coherence=False);
draw.copy_traslate(to_move, -2, 0, N=2, coherence=False);
draw.copy_traslate(to_move, 0, 2, coherence=False);
draw.copy_traslate(to_move, 0, -2, coherence=False);

pt = draw.add_point(0, 3, coherence=False)
ln = draw.add_line(pt, pts[-1], ph=5, coherence=False)
# draw.mirror(ln, 0, 0, 1, 0, coherence=False);
draw.clear_selected()

# %% Add labels and holes
hol = draw.add_hole(-4, 0)
draw.copy_traslate([hol], 4, 0, N=2, coherence=False);
lab1 = draw.add_label(2, 0, 10, 0.2)
draw.copy_rotate([lab1], 90, deg=True, N=3, coherence=False);
draw.add_label( 6, 0, 8, 0.2);
draw.add_label(-6, 0, 8, 0.2);

# %%
draw.plot()

# %% Save geo
# draw.mesh_triangle()
geo.save_geo(draw, 'test_geo_1', 0.2, call_triangle=True, gmsh_mesh=2.2)
