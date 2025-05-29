# %%
from dolomites import tiziano
import geo

# %%
draw = tiziano.drawing()

# %%
draw.load('test_geo_machine.tiz')

# %%
draw.plot()

# %%
geo.save_geo(draw, 'test_geo_machine', 0.010, call_triangle=True, gmsh_mesh=None)
