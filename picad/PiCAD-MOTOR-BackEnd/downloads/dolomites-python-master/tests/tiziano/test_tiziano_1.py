# %% ---------------------------------------------------------------------------
# Test tiziano load capabilities
# ------------------------------------------------------------------------------
# %% Import tizinao from dolomites and other useful modules
import math
from dolomites import tiziano

# %% Init a tiziano drawing instance
draw = tiziano.drawing()

# %% Add points to the drawing
p1 = draw.add_point(0.2, 0)
p2 = draw.add_point(-0.2, 0)

# %% Add lines to the drawing
_, _ = draw.add_arc(p1, p2, math.pi, math.pi/10, 20)
_, _ = draw.add_arc(p1, p2, -math.pi, math.pi/10, 20)

draw.plot()  # visualize drawing with matplotlib

# %% Load test_0 inside current drawing without deleting it
draw.load('test_0.tiz', cl=False)

draw.plot()  # visualize drawing with matplotlib

# %% Fix labels and holes
draw.holes.clear()
lab = draw.add_label(0, 0, 105, 0.1)
hol = draw.add_hole(0.4, 0)

draw.plot()  # visualize drawing with matplotlib

# %% Save and Mesh
draw.save('test_1.tiz')
_ = draw.mesh_triangle()
draw.save_mesh('test_1.msh')
