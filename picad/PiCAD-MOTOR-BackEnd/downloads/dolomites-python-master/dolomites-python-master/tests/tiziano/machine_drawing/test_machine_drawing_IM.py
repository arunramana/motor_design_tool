# %% markdown
## Test machine_drawing module capabilities

### Induction Machine drawing - Example
#### Induction motor IEC 100.4-150x90 standard lamination geometry.

# %%
# Import useful modules
from dolomites import tiziano
import machine_drawing as draw_em

# %%
# Define unit length
mm = 1e-3  # millimeters

# %%
# Induction Motor geometry data ================================================
# -- From IEC 100.4-150x90 datasheet --

r = 1         # inner rotor, outer stator
g = 0.2*mm    # airgap thickness

# Stator ***********************************************************************
Ds = 90*mm     # stator inner diameter
Dse = 150*mm   # stator outer diameter
Qs = 36        # stator number of slots
# Stator slot (trapezoidal slots - rectangular teeth) ---------------------------
wt_s = 4.15*mm   # stator teeth width
wso = 2.3*mm     # slot opening width
hso = 0.65*mm    # slot opening height
hs1 = 0*mm       # slot first height after opening (not needed in this case)
R1_s = 2.05*mm   # rounding radius of corner near slot opening
hwed_s = 1*mm    # wedge height
hs = 16.31*mm    # slot total height
R2_s = 1.2*mm    # roudning radius of corner far from slot opening
# ******************************************************************************

# Rotor (squirrel cage) ********************************************************
Dre = 30*mm      # rotor inner diameter (shaft diameter)
Qr = 28          # squirrel cage number of slots
# Rotor slot (trapezoidal slots- rectangular teeth) ----------------------------
wt_r = 5.1*mm    # rotor teeth width
wro = 1.1*mm     # slot opening width
hro = 0*mm       # slot opening height
hr1 = 0*mm       # slot first height after opening (not needed in this case)
R1_r = 2.18*mm   # rounding radius of corner near airgap
hwed_r = 0*mm    # wedge height (no wedge for casted bar conductors)
hr = 17.38*mm    # slot total height
R2_r = 0.56*mm   # rounding radius of corner far from airgap
# ******************************************************************************

# %%
# Compute needed geometrical data not in datasheet
Dr = Ds - 2*g
alpha_s = 360/Qs
alpha_r = 360/Qr
wsi = draw_em.calc_w(Ds + 2*r*(hso+hs1), Qs, wt_s)
wse = draw_em.calc_w(Ds + 2*r*hs, Qs, wt_s)
wri = draw_em.calc_w(Dr - 2*r*(hro+hr1), Qr, wt_r)
wre = draw_em.calc_w(Dr - 2*r*hr, Qr, wt_r)
print('Rotor diamter toward airgap: Dr = ', Dr*1e3, ' mm')
print('Stator slot width after opening: wsi = ', wsi*1e3, ' mm')
print('Stator slot width far from opening: wse = ', wse*1e3, ' mm')
print('Rotor slot width after opening: wri = ', wri*1e3, ' mm')
print('Rotor slot width far from opening: wre = ', wre*1e3, ' mm')
print('Stator slot mechanical angle: alpha_s = ', alpha_s, ' deg')
print('Rotor slot mechanical angle: alpha_r = ', alpha_r, ' deg')

# %%
# Init a tiziano drawing instance for IM
draw1 = tiziano.drawing()

# %%
# Draw stator slot
slot_pts_s, slot_lns_s = draw_em.slot(draw1, r, Ds, alpha_s, wso, hso, wsi, hs1, hwed_s, wse, hs, R1=R1_s, sol1=2, R2=R2_s, sol2=1, ang_pos=alpha_s/2, coherence=False)
print(len(slot_pts_s))
print(len(slot_lns_s))
print(len(draw1.points))
print(len(draw1.lines))
draw1.select(slot_pts_s + slot_lns_s)
draw1.plot()
draw1.clear_selected()

# %%
# Draw rotor slot
slot_pts_r, slot_lns_r = draw_em.slot(draw1, -r, Dr, alpha_r, wro, hro, wri, hr1, hwed_r, wre, hr, R1=R1_r, sol1=2, R2=R2_r, sol2=1, ang_pos=0, coherence=False)
print(len(slot_pts_r))
print(len(slot_lns_r))
print(len(draw1.points))
print(len(draw1.lines))
draw1.select(slot_pts_r + slot_lns_r)
draw1.plot()
draw1.clear_selected()

# %%
# Draw stator using 'stator' machine_drawing function
stat_pts, stat_lns, stat_Az0lns = draw_em.stator(draw1, slot_pts_s, slot_lns_s, Dse, Qs, alpha_s, coherence=False)
print(len(stat_pts))
print(len(stat_lns))
print(len(draw1.points))
print(len(draw1.lines))
draw1.select(stat_pts + stat_lns)
draw1.plot()
draw1.clear_selected()

# %%
# Draw squirrel cage rotor using 'stator' machine_drawing function again
rot_pts, rot_lns, rot_Az0lns = draw_em.stator(draw1, slot_pts_r, slot_lns_r, Dre, Qr, alpha_r, coherence=False)
print(len(rot_pts))
print(len(rot_lns))
print(len(draw1.points))
print(len(draw1.lines))
draw1.select(rot_pts + rot_lns)
draw1.plot()
draw1.clear_selected()

# %%
# Divide airgap into 5 equal layers according to GetDP formulation
gap_pts, gap_lns, lns_int, lns_out = draw_em.airgap(draw1, g, Dr)
print(len(gap_pts))
print(len(gap_lns))
print(len(lns_int))
print(len(lns_out))
# draw1.select(lns_int)
# draw1.plot()
# draw1.clear_selected()
# draw1.select(lns_out)
# draw1.plot()
# draw1.clear_selected()

# %%
# Assign labels and physical tags to stator according to GetDP formulation
msh_areas_stat = [wt_s**2/4, min(wsi, wse)*(hs-hso-hwed_s)/16, wso*hso/6]
lab_stat = draw_em.ph_stator(draw1, r, stat_Az0lns, Ds, Dse, Qs, alpha_s, wso, hso, hwed_s, hs, msh_areas_stat, ang_pos=alpha_s/2)
# draw1.plot()

# %%
# Assign labels and physical tags to rotor according to GetDP formulation
msh_areas_rot = [wt_r**2/4, min(wri, wre)*(hr-hro-hwed_r)/16, wro*hro/6]
lab_rot = draw_em.ph_IMrotor(draw1, r, rot_Az0lns, Dr, Dre, Qr, alpha_r, wro, hro, hwed_r, hr, msh_areas_rot, ang_pos=0)
# draw1.plot()

# %%
# Assign labels and physical tags to airgap layers according to GetDP formulation
gap_msh_area = (g/5)**2/2
lab_gap = draw_em.ph_airgap(draw1, lns_int, lns_out, r, g, Dr, gap_msh_area)

# %%
# Plot final tiziano IM drawing
draw1.plot()

# %%
# Mesh IM drawing
mesh1 = draw1.mesh_triangle()

# %%
# Save mesh as gmsh file for GetDP solver
draw1.save_mesh('test_IM.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing IM - mesh model](test_IM.png)
