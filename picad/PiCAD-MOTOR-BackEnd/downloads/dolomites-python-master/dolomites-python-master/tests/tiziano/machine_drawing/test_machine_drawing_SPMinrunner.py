# %% markdown
## Test machine_drawing module capabilities

### Inrunner SPM Machine drawing - Example

# %%
# Import useful modules
from dolomites import tiziano
import machine_drawing as draw_em

# %%
# Define unit length
mm = 1e-3  # millimeters

# %%
# SPM Motor geometry data ======================================================

r = 1         # inner rotor, outer stator
g = 0.4*mm    # airgap thickness

# Stator ***********************************************************************
Ds = 103*mm     # stator inner diameter
Dse = 170*mm    # stator outer diameter
Qs = 36         # stator number of slots
# Stator slot (trapezoidal slots - rectangular teeth) ---------------------------
wt_s = 4.9*mm    # stator teeth width
wso = 2.5*mm     # slot opening width
hso = 0.64*mm    # slot opening height
hs1 = 0*mm       # slot first height after opening (not needed in this case)
R1_s = 2.25*mm   # rounding radius of corner near slot opening
hwed_s = 0*mm    # wedge height
hs = 17.3*mm     # slot total height
R2_s = 1.2*mm    # roudning radius of corner far from slot opening
# ******************************************************************************

# Rotor SPM ********************************************************************
p = 3                # number of pole pairs
Dre = 30*mm          # rotor inner diameter (shaft diameter)
typ_rot = 'norm'     # SPM type (normal or consequent-pole)
# Surface mounted Permanent Magnets --------------------------------------------
tm = 5*mm            # magnets thickness
tin = 0*mm           # magnets thicness buried inside rotor iron
alpha_tilt = 0       # magnets mechanical tilting angle in deg
typ_pm = 'radial'    # magnets shape (radial or square)
# ******************************************************************************

# %%
# Compute needed geometrical data not in datasheet
Dr = Ds - 2*g
alpha_s = 360/Qs
alpha_p = 360/(2*p)
alpha_m = alpha_p*2/3
wsi = draw_em.calc_w(Ds + 2*r*(hso+hs1), Qs, wt_s)
wse = draw_em.calc_w(Ds + 2*r*hs, Qs, wt_s)
print('Rotor diamter toward airgap (magnets thickness included): Dr = ', Dr*1e3, ' mm')
print('Stator slot width after opening: wsi = ', wsi*1e3, ' mm')
print('Stator slot width far from opening: wse = ', wse*1e3, ' mm')
print('Stator slot mechanical angle: alpha_s = ', alpha_s, ' deg')
print('Pole mechanical angle: alpha_p = ', alpha_p, ' deg')
print('Mechanical angle covered by magnet: alpha_m = ', alpha_m, ' deg')

# %%
# Init a tiziano drawing instance for SPM
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
# Assign labels and physical tags to stator according to GetDP formulation
msh_areas_stat = [wt_s**2/4, min(wsi, wse)*(hs-hso-hwed_s)/16, wso*hso/6]
lab_stat = draw_em.ph_stator(draw1, r, stat_Az0lns, Ds, Dse, Qs, alpha_s, wso, hso, hwed_s, hs, msh_areas_stat, ang_pos=alpha_s/2)
# draw1.plot()

# %%
# Save stator drawing as tiziano file
draw1.save('test_stator_36_170x103.tiz')

# %%
# Draw rotor SPM pole
pole_pts_r, pole_lns_r = draw_em.SPM_pole(draw1, r, Dr, alpha_p, alpha_m, tm, alpha_tilt=alpha_tilt, tin=tin, ang_pos=0, typ_pm=typ_pm, coherence=False)
print(len(pole_pts_r))
print(len(pole_lns_r))
print(len(draw1.points))
print(len(draw1.lines))
draw1.select(pole_pts_r + pole_lns_r)
draw1.plot()
draw1.clear_selected()

# %%
# Draw SPM rotor using 'SPM_rotor' machine_drawing function
rot_pts, rot_lns, rot_Az0lns = draw_em.SPM_rotor(draw1, pole_pts_r, pole_lns_r, Dre, p, alpha_p, ang_pos=0, typ_rot=typ_rot, coherence=False)
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
# Assign labels and physical tags to rotor according to GetDP formulation
msh_areas_rot = [abs(Dr/2 - Dre/2)**2/20, tm**2/4,  (tm-tin)**2/4]
lab_rot = draw_em.ph_SPMrotor(draw1, r, rot_Az0lns, Dr, Dre, p, alpha_p, alpha_m, tm, msh_areas_rot, alpha_tilt=alpha_tilt, tin=tin, ang_pos=0, typ_pm=typ_pm, typ_rot=typ_rot)
# draw1.plot()

# %%
# Assign labels and physical tags to airgap layers according to GetDP formulation
gap_msh_area = (g/5)**2/2
lab_gap = draw_em.ph_airgap(draw1, lns_int, lns_out, r, g, Dr, gap_msh_area)

# %%
# Plot final tiziano SPM drawing
draw1.plot()

# %%
# Mesh SPM drawing
mesh1 = draw1.mesh_triangle()

# %%
# Save mesh as gmsh file for GetDP solver
draw1.save_mesh('test_SPMinrunner_1.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing inrunner SPM - mesh model 1](test_SPMinrunner_1.png)

# %% markdown
### Inrunner SPM with **titled** and partially **buried magnets** - Example

# %%
# Change magnets' geometrical parameters
alpha_tilt = 5  # magnet tilting angle
tin = 3*mm      # magnet thickness inside rotor iron

# %%
# Init a new tiziano drawing instance
draw2 = tiziano.drawing()

# %%
# Load stator
draw2.load('test_stator_36_170x103.tiz')
draw2.plot()

# %%
# Draw new rotor SPM pole
pole_pts_r, pole_lns_r = draw_em.SPM_pole(draw2, r, Dr, alpha_p, alpha_m, tm, alpha_tilt=alpha_tilt, tin=tin, ang_pos=0, typ_pm=typ_pm, coherence=False)

# %%
# Draw SPM rotor using 'SPM_rotor' machine_drawing function
rot_pts, rot_lns, rot_Az0lns = draw_em.SPM_rotor(draw2, pole_pts_r, pole_lns_r, Dre, p, alpha_p, ang_pos=0, typ_rot=typ_rot, coherence=False)
draw2.plot()

# %%
# Divide airgap into 5 equal layers according to GetDP formulation
gap_pts, gap_lns, lns_int, lns_out = draw_em.airgap(draw2, g, Dr)

# %%
# Assign labels and physical tags to rotor according to GetDP formulation
msh_areas_rot = [abs(Dr/2 - Dre/2)**2/20, tm**2/4,  (tm-tin)**2/4]
lab_rot = draw_em.ph_SPMrotor(draw2, r, rot_Az0lns, Dr, Dre, p, alpha_p, alpha_m, tm, msh_areas_rot, alpha_tilt=alpha_tilt, tin=tin, ang_pos=0, typ_pm=typ_pm, typ_rot=typ_rot)

# %%
# Assign labels and physical tags to airgap layers according to GetDP formulation
gap_msh_area = (g/5)**2/2
lab_gap = draw_em.ph_airgap(draw2, lns_int, lns_out, r, g, Dr, gap_msh_area)

# %%
# Plot final tiziano SPM drawing
draw2.plot()

# %%
# Mesh SPM drawing
mesh2 = draw2.mesh_triangle()

# %%
# Save mesh as gmsh file for GetDP solver
draw2.save_mesh('test_SPMinrunner_2.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing inrunner SPM - mesh model 2](test_SPMinrunner_2.png)

# %% markdown
### Inrunner consequent-pole SPM with squared magnets - Example

# %%
# Change rotor and magnets' geometrical parameters
typ_rot = 'cons'   # consequent-pole rotor
typ_pm = 'square'  # squared magnets shape
alpha_tilt = 0     # magnet tilting angle
tin = 0*mm         # magnet thickness inside rotor iron

# %%
# Init a new tiziano drawing instance
draw3 = tiziano.drawing()

# %%
# Load stator
draw3.load('test_stator_36_170x103.tiz')
# draw3.plot()

# %%
# Draw rotor SPM pole PAIR
pole_pts_r, pole_lns_r = draw_em.SPM_pole(draw3, r, Dr, 2*alpha_p, alpha_m, tm, alpha_tilt=-alpha_p/2 + alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, coherence=False)
# draw3.plot()

# %%
# Draw SPM rotor using 'SPM_rotor' machine_drawing function
rot_pts, rot_lns, rot_Az0lns = draw_em.SPM_rotor(draw3, pole_pts_r, pole_lns_r, Dre, p, alpha_p, ang_pos=0, typ_rot=typ_rot, coherence=False)
# draw3.plot()

# %%
# Divide airgap into 5 equal layers according to GetDP formulation
gap_pts, gap_lns, lns_int, lns_out = draw_em.airgap(draw3, g, Dr)

# %%
# Assign labels and physical tags to rotor according to GetDP formulation
msh_areas_rot = [abs(Dr/2 - Dre/2)**2/20, tm**2/4,  (tm-tin)**2/4]
lab_rot = draw_em.ph_SPMrotor(draw3, r, rot_Az0lns, Dr, Dre, p, alpha_p, alpha_m, tm, msh_areas_rot, alpha_tilt=alpha_tilt, tin=tin, ang_pos=0, typ_pm=typ_pm, typ_rot=typ_rot)

# %%
# draw3.plot()

# %%
# Assign labels and physical tags to airgap layers according to GetDP formulation
gap_msh_area = (g/5)**2/2
lab_gap = draw_em.ph_airgap(draw3, lns_int, lns_out, r, g, Dr, gap_msh_area)

# %%
# Plot final tiziano SPM drawing
draw3.plot()

# %%
# Mesh SPM drawing
mesh3 = draw3.mesh_triangle()

# %%
# Save mesh as gmsh file for GetDP solver
draw3.save_mesh('test_SPMinrunner_3.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing inrunner SPM - mesh model 3](test_SPMinrunner_3.png)
