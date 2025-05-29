# %% Import useful modules
import numpy as np
from dolomites import tiziano

# %% Init a new tiziano drawing instance
draw = tiziano.drawing()

# %% Simple trial drawing
# add points
p1 = draw.add_point(1, 1)
p2 = draw.add_point(-1, 1)
p3 = draw.add_point(-1, -1)
p4 = draw.add_point(1, -1)
p01 = draw.add_point(0.5, 0)
p02 = draw.add_point(-0.5, 0)

# add lines
lns = []
lns.extend(draw.add_line(p1, p2))
lns.extend(draw.add_line(p2, p3))
lns.extend(draw.add_line(p3, p4))
lns.extend(draw.add_line(p4, p1))

# add arcs (discretized)
arc1_pts, arc1_lns = draw.add_arc(p01, p02, np.pi, np.pi/10)
arc2_pts, arc2_lns = draw.add_arc(p01, p02, -np.pi, np.pi/10)

# add labels and holes
lab = draw.add_label(0, 0.8, 100, 0.2)
hol = draw.add_hole(0, 0)

# %% Plot picture
draw.plot()

# %% Move all elements
to_move_all = [p1, p2, p3, p4] + arc1_pts + arc2_pts + lns + arc1_lns + arc2_lns + [lab] + [hol]
new_pts, new_lns, new_labs, new_hols = draw.move_traslate(to_move_all, -1, -1,
                                                          entity=-1, coherence=True)

# %% Move only some elements
draw.clear(entity=2)
to_move = [p1]
# to_move = [lns[0], lns[2]]
new_pts, new_lns, new_labs, new_hols = draw.move_traslate(to_move, 1, 1,
                                                          entity=-1, coherence=True)
                        
# %% Print info to screen and Plot picture
print(len(new_pts))
print(len(new_lns))
print(len(new_labs))
print(len(new_hols))
print(len(draw.points))
print(len(draw.lines))
print(len(draw.labels))
print(len(draw.holes))

draw.select(new_pts+new_lns+new_labs+new_hols)
draw.plot()
draw.clear_selected()

# %% Another simple trial for traslations
draw1 = tiziano.drawing()

# %%
x = [-2, -1, 1, 2]
y = [ 0,  0, 0, 0]
pts = []
for ii, xx in enumerate(x):
    pts.append(draw1.add_point(xx, y[ii]))
lns = []
for jj, pp in enumerate(pts[:-1]):
    lns.extend(draw1.add_line(pp, pts[jj+1]))

# %%
po1 = draw1.add_point(-2, 2)
po2 = draw1.add_point(-2, 0.5)
lo  = draw1.add_line(po1, po2)[0]

out_pts, out_lns, _, _ = draw1.copy_traslate([lo], 0.5, 0, N=8, entity=1)

# %% Plot picture
print(len(out_pts))
print(len(out_lns))
draw1.select(out_pts + out_lns)
draw1.plot()
draw1.clear_selected()

# %% Move geom. objects to test tiziano move_traslate() method
to_move = out_lns
# to_move = pts[:]
new_pts, new_lns, _, _ = draw1.move_traslate(to_move, 0.5, -1,
                                            entity=1, coherence=True)

# %% Plot picture
print(len(new_pts))
print(len(new_lns))
print(len(draw1.points))
print(len(draw1.lines))
draw1.select(new_pts + new_lns)
draw1.plot()
draw1.clear_selected()

# %% Check for bugged returned lines and points
print(new_pts)
print(new_lns)

# %%
# --> FIXED on 20th March 2022 (immediately). ----------------------------------
# 20th March 2022
# BUG: There is a bug in move_traslate() method and most probabily in all
# transformation methods (copy_traslate, copy_rotate, mirror, ...)
# --> in some particular cases transformation methods return more lines/points
# than the ones which are actually newly created. This happens frequently when
# many coherence checks are done due to many intersections between points and
# lines.

# ==> This bug should not be an issue for meshing purposes as long as the bugged
#     methods return more geometry objects that the one actually newly created.
#     However, a particular care must be taken when using these methods, until
#     the bug is fixed.
# ------------------------------------------------------------------------------

# %% Test move rotate
draw2 = tiziano.drawing()

# %% Simple trial geometry (a triangle)
pt1 = draw2.add_point(-1, 0)
pt2 = draw2.add_point(1, 0)
pt3 = draw2.add_point(0, 2)

lns_t1 = draw2.add_line(pt1, pt2)[0]
lns_t2 = draw2.add_line(pt2, pt3)[0]
lns_t3 = draw2.add_line(pt3, pt1)[0]

pto1 = draw2.add_point(-1, 2.5)
pto2 = draw2.add_point(1, 2.5)
lns_to = draw2.add_line(pto1, pto2)[0]

lab_t = draw2.add_label(0, 0.5, 100, 0.5)

# %% Plot picture
draw2.plot()

# %% Move rotate some geometry objects
# to_rot = [pto1, pto2, lns_to]
# to_rot = [pto1, pto2, lns_to, lns_t1, lns_t2, lns_t3, lab_t]
to_rot = [lns_t1, lns_t2, lns_t3, lab_t]
new_pts, new_lns, new_labs, new_hols = draw2.move_rotate(to_rot, -180, deg=True,
                                                         xc=pt3.x, yc=pt3.y, entity=-1)

# %% Plot picture
print(len(new_pts))
print(len(new_lns))
print(len(draw2.points))
print(len(draw2.lines))
draw2.select(new_pts + new_lns)
draw2.plot()
draw2.clear_selected()
