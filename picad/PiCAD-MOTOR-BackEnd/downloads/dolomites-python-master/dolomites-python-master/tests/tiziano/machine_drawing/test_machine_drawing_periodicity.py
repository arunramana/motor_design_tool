# %% markdown
## Test machine_drawing module capabilities

### Inrunner SPM Machine drawing with periodicity - Example

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
# Periodicity
per = 6
Q_sim = int(Qs/per)

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
# draw1.select(slot_pts_s + slot_lns_s)
# draw1.plot()
# draw1.clear_selected()

# %%
# Draw stator using 'stator' machine_drawing function
stat_pts, stat_lns, stat_Az0lns = draw_em.stator(draw1, slot_pts_s, slot_lns_s, Dse, Q_sim, alpha_s, ang_pos=alpha_s/2, coherence=False)
print(len(stat_pts))
print(len(stat_lns))
print(len(draw1.points))
print(len(draw1.lines))
# draw1.select(stat_pts + stat_lns)
# draw1.plot()
# draw1.clear_selected()

# %%
# Assign labels and physical tags to stator according to GetDP formulation
msh_areas_stat = [wt_s**2/4, min(wsi, wse)*(hs-hso-hwed_s)/16, wso*hso/6]
lab_stat = draw_em.ph_stator(draw1, r, stat_Az0lns, Ds, Dse, Q_sim, alpha_s, wso, hso, hwed_s, hs, msh_areas_stat, ang_pos=alpha_s/2)
draw1.plot()

# %%
# Save stator drawing as tiziano file
draw1.save('test_stator_36_170x103_per6.tiz')

# %%
# Draw rotor SPM pole
pole_pts_r, pole_lns_r = draw_em.SPM_pole(draw1, r, Dr, alpha_p, alpha_m, tm, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, coherence=False)
print(len(pole_pts_r))
print(len(pole_lns_r))
print(len(draw1.points))
print(len(draw1.lines))
# draw1.select(pole_pts_r + pole_lns_r)
# draw1.plot()
# draw1.clear_selected()

# %%
# Draw SPM rotor using 'SPM_rotor' machine_drawing function
rot_pts, rot_lns, rot_Az0lns = draw_em.SPM_rotor(draw1, pole_pts_r, pole_lns_r, Dre, p/per, alpha_p, ang_pos=alpha_p/2, typ_rot=typ_rot, coherence=False)
print(len(rot_pts))
print(len(rot_lns))
print(len(draw1.points))
print(len(draw1.lines))
# draw1.select(rot_pts + rot_lns)
# draw1.plot()
# draw1.clear_selected()

# %%
# Assign labels and physical tags to rotor according to GetDP formulation
msh_areas_rot = [abs(Dr/2 - Dre/2)**2/20, tm**2/4, (tm-tin)**2/4]
lab_rot = draw_em.ph_SPMrotor(draw1, r, rot_Az0lns, Dr, Dre, p/per, alpha_p, alpha_m, tm, msh_areas_rot, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, typ_rot=typ_rot)
draw1.plot()

# %%
# Divide airgap into 5 equal layers according to GetDP formulation
gap_pts, gap_lns, lns_int, lns_out = draw_em.airgap(draw1, g, Dr, per=6)
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
# Assign labels and physical tags to airgap layers according to GetDP formulation
gap_msh_area = (g/5)**2/2
lab_gap = draw_em.ph_airgap(draw1, lns_int, lns_out, r, g, Dr, gap_msh_area, ang_pos=alpha_p/2)

# %%
# Close stator with periodic boundaries
lns_stat_r, lns_stat_l = draw_em.close_per(draw1, Ds, Dse, g, [18, 19], 80, per, ang_pos=0)
print(len(lns_stat_r))
print(len(lns_stat_l))
# draw1.plot()

# %%
# Close rotor with periodic boundaries
lns_rot_r, lns_rot_l = draw_em.close_per(draw1, Dr, Dre, g, [16, 17], 80, per, tout=tm, N_air=40, ang_pos=0)
print(len(lns_rot_r))
print(len(lns_rot_l))
# draw1.plot()

# %%
# Plot final tiziano SPM drawing
draw1.plot()

# %%
# Mesh SPM drawing
mesh1 = draw1.mesh_triangle()

# %%
# Close Moving Band at rotor side for periodicity
draw_em.cloe_MB_per(draw1, per, 30)

# %%
# Save mesh as gmsh file for GetDP solver
draw1.save_mesh('test_SPMinrunner_per_1.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing periodicity - mesh model 1](test_SPMinrunner_per_1.png)

# %% markdown
### Inrunner SPM Machine with periodicity - **Another way to do it** -

# Three equal layers in the airgap, instead of the (default) 5.

# %%
# Init a new tiziano drawing instance
draw2 = tiziano.drawing()

# %%
# Draw rotor
pole_pts_r, pole_lns_r = draw_em.SPM_pole(draw2, r, Dr, alpha_p, alpha_m, tm, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, coherence=False)
print(len(pole_pts_r))
print(len(pole_lns_r))
print(len(draw2.points))
print(len(draw2.lines))
# draw1.select(pole_pts_r + pole_lns_r)
# draw1.plot()
# draw1.clear_selected()

rot_pts, rot_lns, rot_Az0lns = draw_em.SPM_rotor(draw2, pole_pts_r, pole_lns_r, Dre, p/per, alpha_p, ang_pos=alpha_p/2, typ_rot=typ_rot, coherence=False)
print(len(rot_pts))
print(len(rot_lns))
print(len(draw2.points))
print(len(draw2.lines))
# draw1.select(rot_pts + rot_lns)
# draw1.plot()
# draw1.clear_selected()

msh_areas_rot = [abs(Dr/2 - Dre/2)**2/20, tm**2/4, (tm-tin)**2/4]
lab_rot = draw_em.ph_SPMrotor(draw2, r, rot_Az0lns, Dr, Dre, p/per, alpha_p, alpha_m, tm, msh_areas_rot, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, typ_rot=typ_rot)
draw2.plot()

# %%
# Draw rotor side airgap layers only
rot_gap_pts, rot_gap_lns, rot_mb_lns = draw_em.airgap_oneside(draw2, r, g, Dr, side='rot', da=1.5, per=6, ang_pos=0, ng=3)

# %%
# and assign them labels and ph tags
gap_msh_area = (g)**2/2
rot_gap_lab = draw_em.ph_airgap_oneside(draw2, rot_mb_lns, r, g, Dr, gap_msh_area, side='rot', ang_pos=alpha_p/2, ng=3)
draw2.plot()

# %%
# Close rotor with periodic boundaries
lns_rot_r, lns_rot_l = draw_em.close_per(draw2, Dr, Dre, g, [16, 17], 80, per, tout=tm, N_air=40, ang_pos=0, ng=3)
draw2.plot()

# %%
# Load stator (without clearing all objects alreay in the drawing! (cl=False))
draw2.load('test_stator_36_170x103_per6.tiz', cl=False)
draw2.plot()

# %%
# Stator side airgap layers
stat_gap_pts, stat_gap_lns, stat_mb_lns = draw_em.airgap_oneside(draw2, r, g, Ds, side='stat', da=1.5, per=6, ang_pos=0, ng=3)
gap_msh_area = (g)**2/2
stat_gap_lab = draw_em.ph_airgap_oneside(draw2, stat_mb_lns, r, g, Ds, gap_msh_area , side='stat', ang_pos=alpha_p/2, ng=3)
# draw2.plot()

# %%
# Close stator with periodic boundaries
lns_stat_r, lns_stat_l = draw_em.close_per(draw2, Ds, Dse, g, [18, 19], 80, per, ang_pos=0, ng=3)
draw2.plot()

# %%
# Mesh with triangle
mesh = draw2.mesh_triangle()

# %%
# Close periodic Moving Band
draw_em.cloe_MB_per(draw2, per)

# %%
# Save mesh in gmsh format
draw2.save_mesh('test_SPMinrunner_per_2.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing periodicity - mesh model 2](test_SPMinrunner_per_2.png)

# %% markdown
### Inrunner SPM Machine with periodicity - **Tilted and buried magnets** -

# %%
# Change magnets' geometrical parameters
alpha_tilt = 5  # magnet tilting angle
tin = 3*mm      # magnet thickness inside rotor iron

# %%
# Change periodicity
per = 3
Q_sim = int(Qs/per)

# %%
# Init an new tiziano drawing instance
draw3 = tiziano.drawing()

# %%
# Draw stator slot
slot_pts_s, slot_lns_s = draw_em.slot(draw3, r, Ds, alpha_s, wso, hso, wsi, hs1, hwed_s, wse, hs, R1=R1_s, sol1=2, R2=R2_s, sol2=1, ang_pos=alpha_s/2)
print(len(slot_pts_s))
print(len(slot_lns_s))
print(len(draw3.points))
print(len(draw3.lines))
# draw1.select(slot_pts_s + slot_lns_s)
# draw1.plot()
# draw1.clear_selected()

# %%
# Draw stator using 'stator' machine_drawing function
stat_pts, stat_lns, stat_Az0lns = draw_em.stator(draw3, slot_pts_s, slot_lns_s, Dse, Q_sim, alpha_s, ang_pos=alpha_s/2)
print(len(stat_pts))
print(len(stat_lns))
print(len(draw3.points))
print(len(draw3.lines))
# draw1.select(stat_pts + stat_lns)
# draw1.plot()
# draw1.clear_selected()

# %%
# Assign labels and physical tags to stator according to GetDP formulation
msh_areas_stat = [wt_s**2/4, min(wsi, wse)*(hs-hso-hwed_s)/16, wso*hso/6]
lab_stat = draw_em.ph_stator(draw3, r, stat_Az0lns, Ds, Dse, Q_sim, alpha_s, wso, hso, hwed_s, hs, msh_areas_stat, ang_pos=alpha_s/2)
draw3.plot()

# %%
# Save stator drawing as tiziano file
draw3.save('test_stator_36_170x103_per3.tiz')

# %%
# Draw rotor SPM pole
pole_pts_r, pole_lns_r = draw_em.SPM_pole(draw3, r, Dr, alpha_p, alpha_m, tm, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, coherence=False)
print(len(pole_pts_r))
print(len(pole_lns_r))
print(len(draw3.points))
print(len(draw3.lines))
# draw1.select(pole_pts_r + pole_lns_r)
# draw1.plot()
# draw1.clear_selected()

# %%
# Draw SPM rotor using 'SPM_rotor' machine_drawing function
rot_pts, rot_lns, rot_Az0lns = draw_em.SPM_rotor(draw3, pole_pts_r, pole_lns_r, Dre, p/per, alpha_p, ang_pos=alpha_p/2)
print(len(rot_pts))
print(len(rot_lns))
print(len(draw3.points))
print(len(draw3.lines))
# draw1.select(rot_pts + rot_lns)
# draw1.plot()
# draw1.clear_selected()

# %%
# Assign labels and physical tags to rotor according to GetDP formulation
msh_areas_rot = [abs(Dr/2 - Dre/2)**2/20, tm**2/4, (tm-tin)**2/4]
lab_rot = draw_em.ph_SPMrotor(draw3, r, rot_Az0lns, Dr, Dre, p/per, alpha_p, alpha_m, tm, msh_areas_rot, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, typ_rot=typ_rot)
# draw3.plot()

# %%
# Divide airgap into 3 equal layers according to GetDP formulation
gap_pts, gap_lns, lns_int, lns_out = draw_em.airgap(draw3, g, Dr, per=3, ng=3)
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
# Assign labels and physical tags to airgap layers according to GetDP formulation
gap_msh_area = (g)**2/2
lab_gap = draw_em.ph_airgap(draw3, lns_int, lns_out, r, g, Dr, gap_msh_area, ang_pos=alpha_p/2, ng=3)

# %%
# Close stator with periodic boundaries
lns_stat_r, lns_stat_l = draw_em.close_per(draw3, Ds, Dse, g, [18, 19], 80, per, ang_pos=0, ng=3)
print(len(lns_stat_r))
print(len(lns_stat_l))
# draw1.plot()

# %%
# Close rotor with periodic boundaries
lns_rot_r, lns_rot_l = draw_em.close_per(draw3, Dr, Dre, g, [16, 17], 80, per, tout=tm-tin, N_air=40, ang_pos=0, ng=3)
print(len(lns_rot_r))
print(len(lns_rot_l))
# draw1.plot()

# %%
# Plot final tiziano SPM drawing
draw3.plot()

# %%
# Mesh SPM drawing
mesh1 = draw3.mesh_triangle()

# %%
# Close Moving Band at rotor side for periodicity
draw_em.cloe_MB_per(draw3, per)

# %%
# Save mesh as gmsh file for GetDP solver
draw3.save_mesh('test_SPMinrunner_per_3.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing periodicity - mesh model 3](test_SPMinrunner_per_3.png)

# %% markdown
### Outrunner SPM Machine with periodicity

# %%
# SPM machine geometry data ====================================================
# 10 poles - 12 slots fractional slots concentrated winding

r = -1        # inner stator, outer rotor
g = 0.4*mm    # airgap thickness

# Stator ***********************************************************************
Ds = 120*mm     # stator outer diameter
Dse = 36*mm     # stator inner diameter
Qs = 12         # stator number of slots
# Stator slot (trapezoidal slots - rectangular teeth) ---------------------------
wt_s = 8*mm     # stator teeth width
wso = 8*mm      # slot opening width
hso = 3.4*mm    # slot opening height
hs1 = 2.5*mm       # slot first height after opening (not needed in this case)
R1_s = 0*mm        # rounding radius of corner near slot opening
hwed_s = 2.5*mm    # wedge height
hs = 30*mm         # slot total height
R2_s = 0.2*mm      # rounding radius of corner far from slot opening
# ******************************************************************************

# Rotor SPM ********************************************************************
p = 5                # number of pole pairs
Dre = 150*mm         # rotor outer diameter (external machine diameter)
typ_rot = 'norm'     # SPM type (normal or consequent-pole)
# Surface mounted Permanent Magnets --------------------------------------------
tm = 4*mm            # magnets thickness
tin = 0*mm           # magnets thicness buried inside rotor iron
alpha_tilt = 0       # magnets mechanical tilting angle in deg
typ_pm = 'radial'    # magnets shape (radial or square)
# ******************************************************************************

# %%
# Compute needed geometrical data not in datasheet
Dr = Ds - 2*r*g
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
# Periodicity
per = 2
Q_sim = int(Qs/per)

# %%
# Init an new tiziano drawing instance
draw4 = tiziano.drawing()

# %%
# Draw stator slot
slot_pts_s, slot_lns_s = draw_em.slot(draw4, r, Ds, alpha_s, wso, hso, wsi, hs1, hwed_s, wse, hs, R1=R1_s, sol1=2, R2=R2_s, sol2=1, ang_pos=alpha_s/2)
print(len(slot_pts_s))
print(len(slot_lns_s))
print(len(draw4.points))
print(len(draw4.lines))
# draw1.select(slot_pts_s + slot_lns_s)
# draw1.plot()
# draw1.clear_selected()

# %%
# Draw stator using 'stator' machine_drawing function
stat_pts, stat_lns, stat_Az0lns = draw_em.stator(draw4, slot_pts_s, slot_lns_s, Dse, Q_sim, alpha_s, ang_pos=alpha_s/2)
print(len(stat_pts))
print(len(stat_lns))
print(len(draw4.points))
print(len(draw4.lines))
# draw1.select(stat_pts + stat_lns)
# draw1.plot()
# draw1.clear_selected()

# %%
# Assign labels and physical tags to stator according to GetDP formulation
msh_areas_stat = [wt_s**2/4, min(wsi, wse)*(hs-hso-hwed_s)/16, wso*hso/6]
lab_stat = draw_em.ph_stator(draw4, r, stat_Az0lns, Ds, Dse, Q_sim, alpha_s, wso, hso, hwed_s, hs, msh_areas_stat, ang_pos=alpha_s/2)
# draw4.plot()

# %%
# Draw rotor SPM pole
pole_pts_r, pole_lns_r = draw_em.SPM_pole(draw4, r, Dr, alpha_p, alpha_m, tm, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, coherence=False)
print(len(pole_pts_r))
print(len(pole_lns_r))
print(len(draw4.points))
print(len(draw4.lines))
# draw1.select(pole_pts_r + pole_lns_r)
# draw1.plot()
# draw1.clear_selected()

# %%
# Draw SPM rotor using 'SPM_rotor' machine_drawing function
rot_pts, rot_lns, rot_Az0lns = draw_em.SPM_rotor(draw4, pole_pts_r, pole_lns_r, Dre, p/per, alpha_p, ang_pos=alpha_p/2)
print(len(rot_pts))
print(len(rot_lns))
print(len(draw4.points))
print(len(draw4.lines))
# draw1.select(rot_pts + rot_lns)
# draw1.plot()
# draw1.clear_selected()


# %%
# Assign labels and physical tags to rotor according to GetDP formulation
msh_areas_rot = [abs(Dr/2 - Dre/2)**2/20, tm**2/4, (tm-tin)**2/4]
lab_rot = draw_em.ph_SPMrotor(draw4, r, rot_Az0lns, Dr, Dre, p/per, alpha_p, alpha_m, tm, msh_areas_rot, alpha_tilt=alpha_tilt, tin=tin, ang_pos=alpha_p/2, typ_pm=typ_pm, typ_rot=typ_rot)
# draw4.plot()

# %%
# Divide airgap into 5 equal layers according to GetDP formulation
gap_pts, gap_lns, lns_int, lns_out = draw_em.airgap(draw4, g, Ds, per=2)
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
# Assign labels and physical tags to airgap layers according to GetDP formulation
gap_msh_area = (g/5)**2/2
lab_gap = draw_em.ph_airgap(draw4, lns_int, lns_out, r, g, Ds, gap_msh_area, ang_pos=alpha_p/2)

# %%
# Close stator with periodic boundaries
lns_stat_r, lns_stat_l = draw_em.close_per(draw4, Ds, Dse, g, [18, 19], 140, per, ang_pos=0)
print(len(lns_stat_r))
print(len(lns_stat_l))
# draw1.plot()

# %%
# Close rotor with periodic boundaries
lns_rot_r, lns_rot_l = draw_em.close_per(draw4, Dr, Dre, g, [16, 17], 80, per, tout=tm, N_air=40, ang_pos=0)
print(len(lns_rot_r))
print(len(lns_rot_l))
# draw1.plot()

# %%
# Plot final tiziano SPM drawing
draw4.plot()

# %%
# Mesh SPM drawing
mesh1 = draw4.mesh_triangle()

# %%
# Close Moving Band at rotor side for periodicity
draw_em.cloe_MB_per(draw4, per)

# %%
# Save mesh as gmsh file for GetDP solver
draw4.save_mesh('test_SPMoutrunner_per_1.msh')

# %% markdown
# This image shows the mesh model just created.  
# It was taken from [Onelab](http://onelab.info/) Gmsh interface.
# ![test machine drawing periodicity - mesh model 4](test_SPMoutrunner_per_1.png)
