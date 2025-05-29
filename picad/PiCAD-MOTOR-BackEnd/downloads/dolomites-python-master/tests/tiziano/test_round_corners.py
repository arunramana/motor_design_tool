# %% markdown
## Test tiziano round_corner methods capability

# %% Import useful modules
import math as mt
from dolomites import tiziano

# %% markdown
### Test round_corner method with corner_arc_solver1 (**solver N.1**)  

# Find corner arc parameters given two tiziano lines that make a corner.

# %% Init a new tiziano drawing
draw1 = tiziano.drawing()

# %% Draw rounding corners trial schema
x = [5, 2, 2, 1, 0]
y = [0, 2, 3, 3, 6]

pts = []

for ii, xx in enumerate(x):
    pts.append(draw1.add_point(xx, y[ii]))

lns = []
for ii, pp in enumerate(pts[:-1], 1):
    lns.extend(draw1.add_line(pp, pts[ii], 10))

print(len(pts))
print(len(lns))
draw1.plot()

# %% Mirroring actions to get final scheme
# coherence False to speed up drawing
pts1, lns1, _, _ = draw1.mirror(lns, 0, 0, 0, 1, entity=1, coherence=False)
pts.extend(pts1)
lns.extend(lns1)
pts2, lns2, _, _ = draw1.mirror(lns, 0, 0, 1, 0, entity=1, coherence=False)
pts.extend(pts2)
lns.extend(lns2)

print(len(pts))
print(len(draw1.points))
print(len(lns))
print(len(draw1.lines))
draw1.plot()

# %% Intersecting lines
l1 = draw1.add_line(pts[0], pts[5])
l2 = draw1.add_line(pts[4], pts[12])
l3 = draw1.add_line(pts[1], pts[13])
l4 = draw1.add_line(pts[6], pts[9])
draw1.plot()

# %% Trying tiziano round_corner capability (1)
# Corner arc solver N.1 --------------------------------------------------------
# Find rounding arc parameters (starting and ending points, arc_angle) given
# two tiziano lines that make a corner
ptsr1, lnsr1 = draw1.round_corner1(lns[0], lns[1], 0.2, 20*mt.pi/180, ph=10, coherence=False)
ptsr2, lnsr2 = draw1.round_corner1(lns[2], lns[1], 0.5, 20*mt.pi/180, ph=10, coherence=False)
ptsr3, lnsr3 = draw1.round_corner1(lns[2], lns[3], 0.3, 20*mt.pi/180, ph=10, coherence=False)
ptsr4, lnsr4 = draw1.round_corner1(lns[7], lns[3], 0.4, 20*mt.pi/180, ph=10)
ptsr5, lnsr5 = draw1.round_corner1(lns[7], lns[6], 0.3, 20*mt.pi/180, ph=10, coherence=False)
ptsr6, lnsr6 = draw1.round_corner1(lns[5], lns[6], 0.5, 20*mt.pi/180, ph=10, coherence=False)
ptsr7, lnsr7 = draw1.round_corner1(lns[4], lns[5], 0.2, 20*mt.pi/180, ph=10, coherence=False)
ptsr8, lnsr8 = draw1.round_corner1(lns[4], lns[12], 0.4, 20*mt.pi/180, ph=10)
ptsr9, lnsr9 = draw1.round_corner1(lns[12], lns[13], 0.2, 20*mt.pi/180, ph=10, coherence=False)
ptsr10, lnsr10 = draw1.round_corner1(lns[14], lns[13], 0.5, 20*mt.pi/180, ph=10, coherence=False)
ptsr11, lnsr11 = draw1.round_corner1(lns[14], lns[15], 0.3, 20*mt.pi/180, ph=10, coherence=False)
ptsr12, lnsr12 = draw1.round_corner1(lns[11], lns[15], 0.4, 20*mt.pi/180, ph=10)
ptsr13, lnsr13 = draw1.round_corner1(lns[11], lns[10], 0.3, 20*mt.pi/180, ph=10, coherence=False)
ptsr14, lnsr14 = draw1.round_corner1(lns[9], lns[10], 0.5, 20*mt.pi/180, ph=10, coherence=False)
ptsr15, lnsr15 = draw1.round_corner1(lns[9], lns[8], 0.2, 20*mt.pi/180, ph=10, coherence=False)
ptsr16, lnsr16 = draw1.round_corner1(lns[8], lns[0], 0.4, 20*mt.pi/180, ph=10)
# %%
# draw1.clear_selected()
draw1.plot()

# %%
print(len(draw1.points))
print(len(draw1.lines))

# %% Verify that mesh is created properly
# Add label to call triangle
lab1 = draw1.add_label(0.5, 0.5, 1, 0.2)
lab2 = draw1.add_label(-0.5, 0.5, 2, 0.2)
lab3 = draw1.add_label(-0.5, -0.5, 3, 0.2)
lab4 = draw1.add_label(0.5, -0.5, 4, 0.2)

# %%
draw1.plot()

# %% Mesh
mesh = draw1.mesh_triangle()

# %% Save as gmsh
draw1.save_mesh('test_round_corners_1.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test round corners - mesh model 1](test_round_corners_1.png)

# %% markdown
### Test corner_arc_solver2 (**solver N.2**)  

# Find corner arc parameters given a tiziano line, a tiziano point and the arc radius.

# %% Init a new tiziano drawing
draw2 = tiziano.drawing()

# %% Draw solver 2 trial schema
x = [0, 0, 3, -3, 4, -1.5]
y = [0, 1, 3,  2, 1,  1.5]

pts = []
lns = []

for ii, xx in enumerate(x[:-1]):
    pts.append(draw2.add_point(xx, y[ii]))

lns.extend(draw2.add_line(pts[1], pts[2], 10))
l = draw2.add_line(pts[1], pts[3], 10)

pt, lns1 = draw2.add_point(x[-1], y[-1], get_lns=True)

pts.append(pt)
lns.extend(lns1)

draw2.plot()

# %% Mirroring action to get final scheme
# coherence disabled to speed up the operation
pts1, lns1, _, _ = draw2.mirror(lns + [pts[-2]], 0, 0, 1, 0, coherence=False)
pts.extend(pts1)
lns.extend(lns1)

draw2.plot()

# %% Draw final lines
lns.extend(draw2.add_line(pts[3], pts[9], 10, coherence=False))
lns.extend(draw2.add_line(pts[4], pts[10], 10, coherence=False))

draw2.plot()

# %% Trying tiziano corner arc solver capability (2)
# Corner arc solver N.2 --------------------------------------------------------
# Find rounding arc parameters (starting and ending points, arc_angle) given
# a tiziano line, a tiziano point and the arc radius
_, _ = draw2.round_corner2(lns[1], pts[0], 1.4, 20*mt.pi/180, ph=10, coherence=False)
_, _ = draw2.round_corner2(lns[0], pts[0], 1.8, 20*mt.pi/180, rp=pts[1], ph=10, coherence=False)
_, _ = draw2.round_corner2(lns[4], pts[0], 1.4, 20*mt.pi/180, ph=10, coherence=False)
_, _ = draw2.round_corner2(lns[3], pts[0], 1.8, 20*mt.pi/180, rp=pts[6], ph=10, coherence=False)

_, _ = draw2.round_corner2(lns[6], pts[5], 1.1, 20*mt.pi/180, rp=pts[3], ph=10, coherence=False)
_, _ = draw2.round_corner2(lns[0], pts[4], 1.6, 20*mt.pi/180, rp=pts[2], ph=10, coherence=False)

_, _ = draw2.round_corner2(lns[6], pts[8], 1.1, 20*mt.pi/180, rp=pts[9], ph=10, coherence=False)
_, _ = draw2.round_corner2(lns[3], pts[10], 1.6, 20*mt.pi/180, rp=pts[7], ph=10, coherence=False)

draw2.plot()

# %% Verify that mesh is created properly
# Add label to call triangle
lab1 = draw2.add_label(-2, 0, 1, 0.2)
lab2 = draw2.add_label(2, 0, 2, 0.2)

draw2.plot()

# %% Mesh
mesh = draw2.mesh_triangle()

# %% Save as gmsh
draw2.save_mesh('test_round_corners_2.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test round corners - mesh model 2](test_round_corners_2.png)
