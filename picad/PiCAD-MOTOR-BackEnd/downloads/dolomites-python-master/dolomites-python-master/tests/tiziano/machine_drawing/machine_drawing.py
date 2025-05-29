# Import useful modules
import numpy as np

# ==============================================================================
# Function to draw a trapezoidal or rectangular 2D slot section
# of a radial flux electric machine using dolomites-tiziano
def slot(draw, r, Dg, alpha_s, wso, hso, wsi, hs1, hwed, wse, hs, R1=0, sol1=2, R2=0, sol2=1, ang_pos=0, close_typ=0, coherence=False):
    """Draw a trapezoidal or rectangular slot of a radial-flux Electric Machine with dolomites-tiziano.

    This function draws a 2D trapezoidal or rectangular slot profile.
    To get a rectangular slot set wsi=wse, otherwise a trapezoidal profile is created.

    It works for both closed and open (or half-open) slots.
    To get closed slots set wso=0 and hso>0, otherwise a (half) open slot is created.

    It works for both inrunner (inner rotor) and outrunner (outer rotor) radial-flux EMs.
    Use r parameter to choose this setting.

    Work for trapezoidal and rectangular slots with or without rounded corners.
    If needed, corners can be rounded using corner_arc_solver1 or corner_arc_solver2
    tiziano drawing methods.
    When rounding radii values too big or too small are provided, corners are not
    rounded.

    This function does not impose any label nor physical tag on lines and points.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to draw the slot profile.
    r : int
         1 --> airgap inside
               (outer stator, or outer rotor with slots, i.e. squirrel cage)
        -1 --> airgap outside
               (inner stator, or inner rotor with slots, i.e. squirrel cage)
    Dg : int or float
        Diameter toward airgap [Dg > 0].
    alpha_s : int or float
        Slot mechanical angle in deg [alpha_s > 0].
    wso : int or float
        Slot opening width [0 <= wso <= min(wsi, wse)].
    hso : int or float
        Slot opening height [0 <= hso <= hs].
    hs1 : int or float
        Slot first height after opening [0 <= hs1 <= hs-hso].
    wsi : int or float
        Slot width near airgap, i.e. internal width [wso <= wsi <= wse].
        Set wsi=wse to get a rectangular slot.
    hwed : int or float
        Slot wedge height [-hso <= hwed <= hs-hso].
        -hso    --> close slot at opening (fully filled slot for bar conductors)
          0     --> close slot just after opening.
        hs-hso  --> do not close the slot (empty slot, no conductors inside slot)
    wse : int or float
        Slot width far from airgap, i.e. external width [wse > 0].
        Set wse=wsi to get a rectangular slot.
    hs : int or float
        Slot total height [hs > 0].
    R1 : int or float, default 0
        Rounding arc radius near slot opening.
        0 for sharp corners (default).
    sol1 : int, default 2
        1 --> round corner near airgap using tiziano corner_arc_solver1 method
        2 --> round corner near airgap using tiziano corner_arc_solver2 method
    R2 : int or float, default 0
        Rounding arc radius far from slot opening.
        0 for sharp corners (default).
    sol2 : int, default 1
        1 --> round corner far from airgap using tiziano corner_arc_solver1 method
        2 --> round corner far from airgap using tiziano corner_arc_solver2 method
    ang_pos : int or float, default 0
        Mechanical angular position with respect to origin in deg which in to
        draw the slot.
    close_typ : int, default 0
        Type for slot closing toward airgap.
        0 --> close slot with a straight line toward airgap (default)
        1 --> close slot with a discretized arc toward airgap
    coherence : bool, default False
        True  --> check coherence when tiziano drawing methods are called
        False --> do not check coherence when tiziano drawing methods are called

    Returns
    -------
    pts : list of tiziano points
        Tiziano points added to the drawing that form the slot profile.
    lns : list of tiziano lines
        Tiziano lines added to the drawing that form the slot profile.
    """
    # compute radius toward air gap
    Rg = Dg/2
    # compute half slot angle in radians
    alpha_s2 = alpha_s/2*np.pi/180
    # compute half slot opening angle in radians
    alpha_so2 = np.arcsin(wso/Dg)

    # compute rotation matrix for specified angular position
    ang_pos = ang_pos*np.pi/180  # angular position in radians
    pos_cos = np.cos(ang_pos)
    pos_sin = np.sin(ang_pos)
    rot_mat = np.array([[pos_cos, -pos_sin],
                        [pos_sin,  pos_cos]])

    # COMPUTE COORDINATES of KEY POINTS
    x1 = Rg*np.cos(alpha_s2)
    y1 = Rg*np.sin(alpha_s2)

    x2 = Rg*np.cos(alpha_so2)
    y2 = wso/2

    x3 = x2 + r*hso
    y3 = y2

    x4 = x3 + r*hs1
    y4 = wsi/2

    x5 = x2 + r*hs
    y5 = wse/2

    x6 = x5
    y6 = 0

    x_aux = x3 + r*hwed
    y_aux = 1.01*max(y5, y4)

    # apply rotation to specified angular position
    xy = np.array([[x1, x2, x3, x4, x5, x6, x_aux],
                   [y1, y2, y3, y4, y5, y6, y_aux]])
    xy = rot_mat @ xy

    # ADD POINTS and LINES of HALF SLOT
    pts = []  # list for key points added to the drawing
    lns = []  # list for key lines added to the drawing
    for _, coord in enumerate(xy[:, :-1].T):
        pts.append(draw.add_point(coord[0], coord[1], coherence=coherence))

    # auxiliary point to close slot according to wedge
    # there is no need to check coherence here because this point will be removed
    p_aux = draw.add_point(xy[0, -1], xy[1, -1], coherence=False)

    for ii, pp in enumerate(pts[1:-1], 1):
        lns.extend(draw.add_line(pp, pts[ii+1]))

    # draw arc toward airgap (1.5 deg fixed discretization angle)
    arc_pts1, arc_lns1 = draw.add_arc(pts[1], pts[0], alpha_s2-alpha_so2,
                                      1.5*np.pi/180, coherence=coherence)

    # round corner near slot opening, if needed
    # (15 deg fixed for rounding arc discretization)
    if len(lns) < 3:
        arc_pts2 = []
        arc_lns2 = []
    else:
        if pts[-1] == pts[-2]:
            l1 = lns[-2]
            l2 = lns[-1]
        else:
            l1 = lns[-3]
            l2 = lns[-2]
        if sol1 == 1:
            arc_pts2, arc_lns2 = draw.round_corner1(l1, l2, R1, 15*np.pi/180,
                                                    coherence=coherence)
        else:
            arc_pts2, arc_lns2 = draw.round_corner2(l2, pts[2], R1, 15*np.pi/180,
                                                    rp=pts[3], coherence=coherence)

    # round corner far from slot opening, if needed
    # (20 deg fixed for rounding arc discretization)
    # if corner near airgap was rounded, input line changes for second rounding
    if pts[-1] == pts[-2]:
        arc_pts3 = []
        arc_lns3 = []
    else:
        l1 = lns[-2]
        if len(arc_lns2) != 0:
            if sol1 == 1:
                l1 = arc_lns2[int(-1.5 - r*0.5)]
            else:
                l1 = arc_lns2[-1]
        if sol2 == 1:
            arc_pts3, arc_lns3 = draw.round_corner1(l1, lns[-1], R2, 20*np.pi/180,
                                                    coherence=coherence)
        else:
            arc_pts3, arc_lns3 = draw.round_corner2(l1, pts[-1], R2, 20*np.pi/180,
                                                    rp=pts[-2], coherence=coherence)

    # remove useless geometry objects (points, lines) from lists
    # if corners rounding is successful
    if len(arc_lns2) != 0:
        if sol1 == 2:
            pts.pop(-4)
        pts.pop(-3)
        lns.pop(-3)
        lns.pop(-2)

    if len(arc_lns3) != 0:
        pts.pop(-2)
        if len(arc_lns2) == 0:
            lns.pop(-2)
        lns.pop(-1)
        if len(arc_lns2) != 0:
            arc_lns2.remove(l1)
        if sol2 == 2:
            pts.pop(-1)

    # join lists to get half slot profile points and lines
    pts.extend(arc_pts1[1:-1] + arc_pts2 + arc_pts3)
    lns.extend(arc_lns1 + arc_lns2 + arc_lns3)

    # mirror half section to get entire slot profile
    mirror_pts, mirror_lns, _, _ = draw.mirror(lns + [p_aux],
                                               0, 0, pos_cos, pos_sin,
                                               coherence=coherence)

    # join lists to get full slot profile points and lines
    pts.extend(mirror_pts[:-1])
    lns.extend(mirror_lns)

    # closing line --> close slot opening toward airgap
    if close_typ == 0:  # close slot with a straight line
        lns.extend(draw.add_line(pts[1], mirror_pts[0], coherence=coherence))
    else:               # close slot with a discretized arc
        arc_cl_pts, arc_cl_lns = draw.add_arc(mirror_pts[0], pts[1], 2*alpha_so2,
                                              1.5*np.pi/180, coherence=coherence)
        pts.extend(arc_cl_pts[1:-1])
        lns.extend(arc_cl_lns)

    # closing line --> slot wedge
    if hwed not in (-hso, hs-hso):
        # here coherence must be checked,
        # because this closing line is created by intersection
        cl, cl_pts = draw.add_line(p_aux, mirror_pts[-1], get_pts=True)

        # fix line and points ouput lists
        for pp in cl_pts:
            to_remove = [ll for ll in lns if ll.shortest_distance(pp) < 1e-10]
            for ll in to_remove:
                lns.remove(ll)

        lns.append(cl[-2])
        if len(cl) == 7:
            lns.extend([cl[0], cl[1], cl[3], cl[4]])

        pts.extend(cl_pts)

    # remove auxiliary poinst and connected lines
    draw.remove_point([p_aux, mirror_pts[-1]])

    # remove double point instances from output list, if any
    # points order in output list is changed using set command!!!
    pts = list(set(pts))  # !!!

    return pts, lns

# ==============================================================================
# Function to draw a 2D stator of a radial flux electric machine from a slot profile
# using dolomites-tiziano
def stator(draw, slot_pts, slot_lns, Dse, Qs, alpha_s, ang_pos=0, coherence=False):
    """Draw a stator or a squirrel cage rotor of a radial-flux Electric Machine with dolomites-tiziano.

    This function draws a 2D stator or squirrel cage rotor profile made by Qs equal slots.
    Work for both outer and inner stators/rotors of radial-flux EMs.

    The stator is drawn by copy rotation of a slot profile,
    and then it is closed using arcs when needed.
    (external diameter for outer stators, shaft and/or passing clamp for inner stators)

    This function does not impose any label nor physical tag on lines and points.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to draw the stator.
        A 2D slot profile should be already in the input drawing instance.
    slot_pts : list of tiziano points
        Tiziano points that form the slot profile.
        These points should be already in the input drawing instance.
    slot_lns : list of tiziano lines
        Tiziano lines that form the slot profile.
        These lines should be already in the input drawing instance.
    Dse : int or float
        Stator external diameter [Dse >= 0].
        "External" means not toward airgap.
        External diameter for outer stators [Dse > Dg].
        Shaft and/or passing clamp diameter for inner stators [0 < Dse < Dg].
        It can be also set to zero if one does not want to close the stator (outer stators)
        or when shaft/passing clamp is not necessary in the drawing (inner stators).
    Qs : int
        Number of stator equal slots in the drawing.
        Copy rotate Qs-1 slots because one is supposed already in the input
        drawing instance.
    alpha_s : int or float
        Slot mechanical angle in deg [alpha_s > 0].
    ang_pos : int or float, default 0
        Mechanical angular position with respect to origin in deg which in the
        first slot was drawn.
    coherence : bool, default False
        True  --> check coherence when tiziano drawing methods are called
        False --> do not check coherence when tiziano drawing methods are called

    Returns
    -------
    stator_pts : list of tiziano points
        Tiziano points added to the drawing that form the stator.
    stator_lns : list of tiziano lines
        Tiziano lines added to the drawing that form the stator.
    arc_lns1 + arc_lns2: list of tiziano lines
        Tiziano lines that close the stator.
        Useful to easily impose physical tags for boundary conditions in next steps.
        Empty if Dse=0.
    """
    # angular position of first slot in rad
    ang_pos = ang_pos*np.pi/180

    # take all lines that form a slot
    # and copy rotate them with their endpoins to get the stator
    pts, lns, _, _ = draw.copy_rotate(slot_lns, alpha_s, deg=True, N=Qs-1, entity=1, coherence=coherence)

    # stator angle to be drawn [rad]
    alpha_sim = alpha_s*Qs*np.pi/180

    # if an external diameter is provided,
    if Dse != 0:
        if alpha_sim == 2*np.pi:  # no periodicity ==> entire machine
            # close the stator with two 180 deg arcs
            p1 = draw.add_point(0, Dse/2, coherence=coherence)
            p2 = draw.add_point(0, -Dse/2, coherence=coherence)
            arc_pts1, arc_lns1 = draw.add_arc(p1, p2, np.pi, 6*np.pi/180, coherence=coherence)
            arc_pts2, arc_lns2 = draw.add_arc(p1, p2, -np.pi, 6*np.pi/180, coherence=coherence)
            # do not consider point already added
            arc_pts2 = arc_pts2[1:-1]
        else:                    # machine periodicity exploited
            # draw external stator profile with only one arc
            start_ang = ang_pos-alpha_s/2*np.pi/180            # [rad]
            end_ang = ang_pos-alpha_s/2*np.pi/180 + alpha_sim  # [rad]
            p1 = draw.add_point(Dse/2*np.cos(start_ang), Dse/2*np.sin(start_ang), coherence=coherence)
            p2 = draw.add_point(Dse/2*np.cos(end_ang), Dse/2*np.sin(end_ang), coherence=coherence)
            arc_pts1, arc_lns1 = draw.add_arc(p1, p2, alpha_sim, 6*np.pi/180, coherence=coherence)
            arc_pts2 = []
            arc_lns2 = []

    else:  # else, no closing arcs are added
        arc_pts1 = []
        arc_lns1 = []
        arc_pts2 = []
        arc_lns2 = []

    # join lists to get stator points and lines output lists
    stator_pts = slot_pts + pts + arc_pts1 + arc_pts2
    stator_lns = slot_lns + lns + arc_lns1 + arc_lns2

    return stator_pts, stator_lns, arc_lns1 + arc_lns2

# ==============================================================================
# Function to draw a Surface mounted Permanent Magnet section of a radial flux electric machine
# using dolomites-tiziano
def SPM_pole(draw, r, Dg, alpha_p, alpha_m, tm, alpha_tilt=0, tin=0, ang_pos=0, typ_pm='radial', coherence=False):
    """Draw Surface mounted Permanent Magnet section of a radial-flux Electric Machine using dolomites-tiziano.

    This function draws a 2D pole section of a radial-flux SPM rotor.
    It works for both inrunner (inner rotor) and outrunner (outer rotor) radial-flux EMs.
    Use r parameter to choose this setting.

    It works for both radial or squared permanent magnet blocks.
    Use typ keyword argument to choose this setting.

    It works for partially buried magnets inside rotor iron.
    Set tin keyword argument properly for buried magnets.

    It works also for tilted magnets and consequent pole topologies, with or without magnets.
    Set alpha_tilt to have a tilting angle and tm=0 to have no magnet.

    For consequent pole topologies an entire pole pair can be drawn with a single
    call to this function. (Suggested method!)
    Set alpha_p as mechanical angle covered an ENTIRE pole PAIR. (2*alpha_p)
    Set alpha_tilt = -alpha_p/2 + alpha_tilt.

    This function does not impose any label nor physical tag on lines and points.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to draw the pole profile.
    r : int
         1 --> airgap outside (inner rotor --> inrunner SPM machine)
        -1 --> airgap inside (outer rotor --> outrunner SPM machine)
    Dg : int or float
        Diameter toward airgap [Dg > 0].
    alpha_p: int or float
        Pole mechnical angle in deg [0 < alpha_p < 180].
        Pole pair mechnical angle in deg in case of consequent pole.
    alpha_m : int or float
        Mechanical angle in deg covered by magnet at the airgap. [0 < alpha_m < alpha_p].
    tm : int or float
        Magnet thickness [tm >= 0].
    alpha_tilt : int or float, default 0
        Magnet tilting angle [0 <= alpha_tilt <= (alpha_p - alpha_m)/2]
    tin : int or float, default 0
        Magnet thickness buried inside iron [0 <= tin <= tm].
        0  --> magnet all in the airgap (default)
        tm --> magnet totally buried inside rotor iron
    ang_pos : int or float, default 0
        Mechanical angular position with respect to origin in deg which in to
        draw the pole profile.
    typ_pm : str, default 'radial'
        Magnet block type, 'radial' (deafult) or 'square'.
    coherence : bool, default False
        True  --> check coherence when tiziano drawing methods are called
        False --> do not check coherence when tiziano drawing methods are called

    Returns
    -------
    pts : list of tiziano points
        Tiziano points added to the drawing that form the pole profile.
    lns : list of tiziano lines
        Tiziano lines added to the drawing that form the pole profile.
    """
    # compute radius toward air gap
    Rg = Dg/2
    # compute half pole angle (alpha_p) in radians
    alpha_p2 = alpha_p/2*np.pi/180
    # compute half magnet covering angle (alpha_m) in radians
    alpha_m2 = alpha_m/2*np.pi/180
    # compute magnet tilt angle in radians
    alpha_tilt = alpha_tilt*np.pi/180
    # compute rotation matrix for specified angular position
    ang_pos = ang_pos*np.pi/180  # angular position in radians
    pos_cos = np.cos(ang_pos)
    pos_sin = np.sin(ang_pos)
    rot_mat = np.array([[pos_cos, -pos_sin],
                        [pos_sin,  pos_cos]])
    # compute magnet tickness in the air gap
    tout = tm - tin

    # COMPUTE KEY POINTS COORDINATES
    x1 = (Rg - r*tout)*np.cos(alpha_p2)
    y1 = (Rg - r*tout)*np.sin(alpha_p2)

    x2 = (Rg - r*tout)*np.cos(alpha_m2 + alpha_tilt)
    y2 = (Rg - r*tout)*np.sin(alpha_m2 + alpha_tilt)

    x3 = (Rg - r*tm)*np.cos(alpha_m2 + alpha_tilt)
    y3 = (Rg - r*tm)*np.sin(alpha_m2 + alpha_tilt)

    x4 = Rg*np.cos(alpha_m2 + alpha_tilt)
    y4 = Rg*np.sin(alpha_m2 + alpha_tilt)

    x5 = Rg*np.cos(-alpha_m2 + alpha_tilt)
    y5 = Rg*np.sin(-alpha_m2 + alpha_tilt)

    x6 = (Rg - r*tm)*np.cos(-alpha_m2 + alpha_tilt)
    y6 = (Rg - r*tm)*np.sin(-alpha_m2 + alpha_tilt)

    x7 = (Rg - r*tout)*np.cos(-alpha_m2 + alpha_tilt)
    y7 = (Rg - r*tout)*np.sin(-alpha_m2 + alpha_tilt)

    x8 = (Rg - r*tout)*np.cos(-alpha_p2)
    y8 = (Rg - r*tout)*np.sin(-alpha_p2)

    x9 = Rg*np.cos(alpha_p2)
    y9 = Rg*np.sin(alpha_p2)

    x10 = Rg*np.cos(-alpha_p2)
    y10 = Rg*np.sin(-alpha_p2)

    # change key points coordinates if the magnet block is squared
    if typ_pm == 'square':
        x2 = x4 - r*tout*np.cos(alpha_tilt)
        y2 = y4 - r*tout*np.sin(alpha_tilt)
        x3 = x4 - r*tm*np.cos(alpha_tilt)
        y3 = y4 - r*tm*np.sin(alpha_tilt)
        x6 = x5 - r*tm*np.cos(alpha_tilt)
        y6 = y5 - r*tm*np.sin(alpha_tilt)
        x7 = x5 - r*tout*np.cos(alpha_tilt)
        y7 = y5 - r*tout*np.sin(alpha_tilt)

    # apply rotation to specified angular position
    xy = np.array([[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10],
                   [y1, y2, y3, y4, y5, y6, y7, y8, y9, y10]])
    xy = rot_mat @ xy

    # ADD KEY POINTS
    pts = []  # list for points
    for pp in range(10):
        pts.append(draw.add_point(xy[0, pp], xy[1, pp], coherence=coherence))

    # ADD LINES
    lns = []  # list for lines
    # this two lines must be added with coherence set to True (default)
    # in the case of magnet partially buried inside iron an intersection with
    # a point occurs
    lns.extend(draw.add_line(pts[2], pts[3]))
    lns.extend(draw.add_line(pts[4], pts[5]))

    # add discretized arcs (2 deg fixed step angle)
    arc1_pts, arc1_lns = draw.add_arc(pts[1], pts[0], alpha_p2 - alpha_m2 - alpha_tilt,
                                      2*np.pi/180, coherence=coherence)
    arc4_pts, arc4_lns = draw.add_arc(pts[7], pts[6], alpha_p2 - alpha_m2 + alpha_tilt,
                                      2*np.pi/180, coherence=coherence)
    if typ_pm == 'radial':
        arc2_pts, arc2_lns = draw.add_arc(pts[4], pts[3], 2*alpha_m2,
                                          2*np.pi/180, coherence=coherence)
        arc3_pts, arc3_lns = draw.add_arc(pts[5], pts[2], 2*alpha_m2,
                                          2*np.pi/180, coherence=coherence)
        arcx_pts = []
        arcx_lns = []
    # if the magnet is squared two arcs become straight lines ...
    if typ_pm == 'square':
        arc2_pts = []
        arc2_lns = draw.add_line(pts[4], pts[3], coherence=coherence)
        arc3_pts = []
        arc3_lns = draw.add_line(pts[5], pts[2], coherence=coherence)
        # ... and an arc is added to close an additional air regions in front of pm
        if r == 1:  # only for inner rotors
            arcx_pts, arcx_lns = draw.add_arc(pts[4], pts[3], 2*alpha_m2,
                                              2*np.pi/180, coherence=coherence)
        else:
            arcx_pts = []
            arcx_lns = []
    # the air regions near megnet inside airgap are closed with discretized arcs
    # if the magnet is not completely buried inside iron (or there is no magnet)
    if tout != 0:
        arc5_pts, arc5_lns = draw.add_arc(pts[3], pts[8], alpha_p2 - alpha_m2 - alpha_tilt,
                                          2*np.pi/180, coherence=coherence)
        arc6_pts, arc6_lns = draw.add_arc(pts[9], pts[4], alpha_p2 - alpha_m2 + alpha_tilt,
                                          2*np.pi/180, coherence=coherence)
    else:
        arc5_pts = []
        arc5_lns = []
        arc6_pts = []
        arc6_lns = []

    # add lines to output list
    lns.extend(arc1_lns + arc2_lns + arc3_lns + arc4_lns + arc5_lns + arc6_lns + arcx_lns)

    # add to output list arcs discretization points
    pts.extend(arc1_pts + arc2_pts + arc3_pts + arc4_pts + arc5_pts + arc6_pts + arcx_pts)
    # remove double point instances from output list, if any
    # points order in output list is changed using set command!!!
    pts = list(set(pts))  # !!!

    return pts, lns

# ==============================================================================
# Function to draw a Surface mounted Permanent Magnet rotor of a radial flux electric machine
# using dolomites-tiziano
def SPM_rotor(draw, input_pts, input_lns, Dre, p, alpha_p, ang_pos=0, typ_rot='norm', coherence=False):
    """Draw a rotor of a SPM radial-flux Electric Machine with dolomites-tiziano.

    This function draws a 2D rotor profile made by p equal pole pairs.
    Work for both outer and inner rotors of radial-flux SPM EMs.
    Work also for asymmetric tilted magnets and consequent pole topologies.

    The rotor is drawn by copy rotation of a pole pair,
    and then it is closed using arcs when needed.
    (external diameter for outer rotors, shaft for inner rotors)

    This function does not impose any label nor physical tag on lines and points.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to draw the stator.
        A 2D slot profile should be already in the input drawing instance.
    input_pts : list of tiziano points
        typ = 'norm'    --> points that form a single SPM pole profile
                            (normal or tilted asymmetric topologies)
        typ = 'cons'    --> points that form already a pole pair
                            (consequent pole topologies)
        These points should be already in the input drawing instance.
    input_lns : list of tiziano lines
        typ = 'norm' --> lines that form a single SPM pole profile
                         (normal or tilted asymmetric topologies)
        typ = 'cons' --> lines that form already a pole pair
                         (consequent pole topologies)
        These lines should be already in the input drawing instance.
    Dre : int or float
        Rotor external diameter [Dre >= 0].
        "External" means not toward airgap.
        External diameter for outer rotors.
        Shaft and/or passing clamp diameter for inner rotors.
        It can be also set to zero if one does not want to close the rotor (outer rotors)
        or when shaft/passing clamp is not necessary in the drawing (inner rotors).
    p : int or float
        Number of rotor equal pole pairs in the drawing.
        Copy rotate p-1 pole pairs, because one is supposed already in the input
        drawing instance for consequent pole topologies, or it is created by
        pole mirroring for other topologies.
        If you want to draw an odd number of poles, specify p as int +- 0.5.
        This last option is effective only for no tilted magnets (pole asymmetry
        is not considered).
    alpha_p : int or float
        Single pole mechanical angle in deg [0 < alpha_p < 180].
    ang_pos : int or float, default 0
        Mechanical angular position with respect to origin in deg which in
        the input single pole profile is drawn in the case of 'normal' type;
        angular position of first pole pair in the case of 'consequent-pole' type.
    typ_rot : str, default 'norm'
        typ_rot = 'norm' --> normal or tilted asymmetric topologies
        typ_rot = 'cons' --> consequent pole topologies
    coherence : bool, default False
        True  --> check coherence when tiziano drawing methods are called
        False --> do not check coherence when tiziano drawing methods are called

    Returns
    -------
    rotor_pts : list of tiziano points
        Tiziano points added to the drawing that form the rotor.
    rotor_lns : list of tiziano lines
        Tiziano lines added to the drawing that form the rotor.
    arc_lns1 + arc_lns2: list of tiziano lines
        Tiziano lines that close the rotor.
        Useful to easily impose physical tags for boundary conditions in next steps.
        Empty if Dre=0.
    """
    # in the case of NO consequent pole topologies and when an integer number of
    # pole pairs is passed, take all lines that form a pole and mirror them to have a pole pair
    # mirror action is needed to consider tilted magnets, that lead to pole asymmetry
    if typ_rot == 'norm' and int(p) == p:
        alpha_1 = alpha_p
        pts_torot, lns_torot, _, _ = draw.mirror(input_lns, 0, 0,
                                                 np.cos((alpha_p/2+ang_pos)*np.pi/180),
                                                 np.sin((alpha_p/2+ang_pos)*np.pi/180),
                                                 entity=1, coherence=coherence)
    # in the case of consequent pole topology or when a float (int +- 0.5) pole pairs is passed,
    # elements to be rotate lists are empy for the moment
    elif typ_rot == 'cons' or (2*p)%2 != 0:
        pts_torot = []
        lns_torot = []
        alpha_1 = 2*alpha_p
    else:
        print('Input error!')
        return [], [], []
    # complete lists of elements to be rotate by adding input lines and points
    # now we have a poles pair!
    pts_torot = pts_torot + input_pts
    lns_torot = lns_torot + input_lns

    alpha_sim = alpha_p*2*p*np.pi/180  # rotor slice to be simulated [rad]
    # if a float (int +- 0.5) number of pole pairs is passed, use p as odd number
    # of poles for the next copy rotation action and half pole angle
    if (2*p)%2 != 0:
        p = int(2*p)
        alpha_p = alpha_p/2
        alpha_1 = 2*alpha_p
    # take all lines that form a pole pairs
    # and copy rotate them with their endpoins to get the SPM rotor
    pts, lns, _, _ = draw.copy_rotate(lns_torot, 2*alpha_p, deg=True, N=int(p-1),
                                      entity=1, coherence=coherence)

    # if an external diameter is provided,
    if Dre != 0:
        if alpha_sim == 2*np.pi:  # no periodicity ==> entire machine
            # close the rotor with two 180 deg arcs
            p1 = draw.add_point(0, Dre/2, coherence=coherence)
            p2 = draw.add_point(0, -Dre/2, coherence=coherence)
            arc_pts1, arc_lns1 = draw.add_arc(p1, p2, np.pi, 6*np.pi/180, coherence=coherence)
            arc_pts2, arc_lns2 = draw.add_arc(p1, p2, -np.pi, 6*np.pi/180, coherence=coherence)
            # do not consider point already added
            arc_pts2 = arc_pts2[1:-1]
        else:                    # machine periodicity exploited
            # draw external rotor profile with only one arc
            start_ang = (ang_pos-alpha_1/2)*np.pi/180            # [rad]
            end_ang = (ang_pos-alpha_1/2)*np.pi/180 + alpha_sim  # [rad]
            p1 = draw.add_point(Dre/2*np.cos(start_ang), Dre/2*np.sin(start_ang), coherence=coherence)
            p2 = draw.add_point(Dre/2*np.cos(end_ang), Dre/2*np.sin(end_ang), coherence=coherence)
            arc_pts1, arc_lns1 = draw.add_arc(p1, p2, alpha_sim, 6*np.pi/180, coherence=coherence)
            arc_pts2 = []
            arc_lns2 = []
    else:  # else, no closing arcs are added
        arc_pts1 = []
        arc_lns1 = []
        arc_pts2 = []
        arc_lns2 = []

    # join lists to get rotor points and lines output lists
    rotor_pts = pts_torot + pts + arc_pts1 + arc_pts2
    rotor_lns = lns_torot + lns + arc_lns1 + arc_lns2

    return rotor_pts, rotor_lns, arc_lns1 + arc_lns2

# ==============================================================================
# Function to assign labels and physical tags to a stator of a radial flux electric machine
# using dolomites-tiziano for GetDP magnetostatic formulation
def ph_stator(draw, r, Az0_lns, Dg, Dse, Qs, alpha_s, wso, hso, hwed, hs, msh_areas, ang_pos=0):
    """Assign labels and physical tags to stator of a radial-flux Electric Machine.

    Labels and physical tags are imposed according to GetDP magnetostatic formulation.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to assign stator labels
        an physical tags.
    r : int
         1 --> airgap inside
               (outer stator)
        -1 --> airgap outside
               (inner stator)
    Az_lns : list of tiziano lines
        Stator boundary lines that need a physical tag for Dirichlet BC.
    Dg : int or float
        Stator diameter toward airgap [Dg > 0].
    Dse : int or float
        Stator external diameter [Dse >= 0].
        "External" means not toward airgap.
        External diameter for outer stators [Dse > Dg].
        Shaft and/or passing clamp diameter for inner stators [0 < Dse < Dg].
        It can be also set to zero if one does not want to close the stator (outer stators)
        or when shaft/passing clamp is not necessary in the drawing (inner stators).
    Qs : int
        Number of stator slots in the drawing [Qs > 0].
    alpha_s : int or float
        Slot mechanical angle in deg [alpha_s > 0].
    wso : int or float
        Slot opening width [0 <= wso].
    hso : int or float
        Slot opening height [0 <= hso <= hs].
    hwed : int or float
        Slot wedge height [-hso <= hwed <= hs-hso].
        -hso    --> slots closed at opening (fully filled slot for bar conductors)
          0     --> slots closed just after opening, (half) open slots
        hs-hso  --> empty slots, no conductors inside slots
    hs : int or float
        Slot total height.
    msh_areas : list of ints or floats
        Maximum areas of mesh triangles in stator iron, slots and air regions.
        [stat_iron_area, stat_slot_area, stat_air_area]
    ang_pos : int or float, default 0
        Mechanical angular position with respect to origin in deg which in the
        first slot was drawn.

    Returns
    -------
    [lab_s] + lab_slot + air_slot + hol : list of tiziano labels and holes
        Tiziano labels and holes added to the drawing for stator meshing.
        lab_s    -> stator iron, lamination (ph=1000)
        lab_slot -> stator slots            (ph=1000+1,..., 1000+Qs)
        air_slot -> stator slot openings    (ph=999)
        hol      -> No-Mesh central region, if needed for inner stators
    """
    # stator radius toward airgap
    Rg = Dg/2
    # mechanical slot angle in radians
    alpha_s = alpha_s*np.pi/180
    # first slot mechanical angular position in radians
    ang_pos = ang_pos*np.pi/180
    # compute stator back iron height
    hsbi = (abs(Dse - Dg) - 2*hs)/2
    # add label in stator iron (stator lamination)
    # compute label coordinates
    if r == 1:     # outer stator, inner rotor
        R_l = Rg + hs + hsbi/2
    elif r == -1:  # inner stator, outer rotor
        R_l = Dse/2 + hsbi/2
    else:
        print('Input error!')
        return []
    x_l = R_l*np.cos(ang_pos)
    y_l = R_l*np.sin(ang_pos)
    lab_s = draw.add_label(x_l, y_l, 1000, msh_areas[0])

    # add in-slots labels
    lab_slot = []  # labels in conducting regions inside slots
    air_slot = []  # labels in insulating/air regions in slots (i.e. openings)
    # labels in conducting regions are added only if slots are not empty
    if hwed != hs-hso:
        R_slot_cond = Rg + r*(hso/2+hwed/2+hs/2)
        for qq in range(0, Qs):
            lab_slot.append(draw.add_label(R_slot_cond*np.cos(ang_pos+qq*alpha_s),
                                           R_slot_cond*np.sin(ang_pos+qq*alpha_s),
                                           1000+qq+1, msh_areas[1]))
    # labels in slot openings are added only if slot are (half) open
    if wso != 0 and hso+hwed != 0:
        R_slot_air = Rg + 3*r*(hso+hwed)/4
        for qq in range(0, Qs):
            air_slot.append(draw.add_label(R_slot_air*np.cos(ang_pos+qq*alpha_s),
                                           R_slot_air*np.sin(ang_pos+qq*alpha_s),
                                           999, msh_areas[-1]))

    # impose physical tags on stator boundary lines
    for ll in Az0_lns:
        ll.ph_tag = 10

    # assign No-Mesh region for inner stators with shaft/passing clamp diameter
    hol = []
    if r == -1 and Dse != 0:
        hol.append(draw.add_hole(0, 0))
    # output label list
    return [lab_s] + lab_slot + air_slot + hol

# ==============================================================================
# Function to assign labels and physical tags to a squirrel cage rotor of a radial flux induction machine
# using dolomites-tiziano for GetDP formulation
def ph_IMrotor(draw, r, Az0_lns, Dg, Dre, Qr, alpha_r, wro, hro, hwed, hr, msh_areas, ang_pos=0):
    """Assign labels and physical tags to squirrel cage rotor of a radial-flux Induction Machine.

    Labels and physical tags are imposed according to GetDP formulation.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to assign stator labels
        an physical tags.
    r : int
         1 --> airgap outside
               (inner rotor)
        -1 --> airgap inside
               (outer rotor)
    Az_lns : list of tiziano lines
        Rotor boundary lines that need a physical tag for Dirichlet BC.
    Dg : int or float
        Rotor diameter toward airgap [Dg > 0].
    Dre : int or float
        Rotor external diameter [Dre >= 0].
        "External" means not toward airgap.
        External diameter for outer rotors [Dre > Dg].
        Shaft and/or passing clamp diameter for inner rotors [0 < Dre < Dg].
        It can be also set to zero if one does not want to close the rotor (outer rotors)
        or when shaft/passing clamp is not necessary in the drawing (inner rotors).
    Qr : int
        Number of rotor slots in the drawing [Qr > 0].
    alpha_r : int or float
        Slot mechanical angle in deg [alpha_r > 0].
    wro : int or float
        Slot opening width [0 <= wro].
    hro : int or float
        Slot opening height [0 <= hro <= hr].
    hwed : int or float
        Slot wedge height [-hro <= hwed <= hr-hro].
        -hso    --> slots closed at opening (fully filled slot for bar conductors)
          0     --> slots closed just after opening, (half) open slots
        hr-hro  --> empty slots, no conductors inside slots
    hr : int or float
        Slot total height.
    msh_areas : list of ints or floats
        Maximum areas of mesh triangles in rotor iron, slots and air regions.
        [rot_iron_area, rot_slot_area, rot_air_area]
    ang_pos : int or float, default 0
        Mechanical angular position with respect to origin in deg which in the
        first slot was drawn.

    Returns
    -------
    [lab_r] + lab_slot + air_slot + hol : list of tiziano labels and holes
        Tiziano labels and holes added to the drawing for stator meshing.
        lab_r    -> rotor iron, lamination (ph=100)
        lab_slot -> rotor slots            (ph=501,..., 500+Qr)
        air_slot -> rotor slot openings    (ph=99)
        hol      -> No-Mesh central region, if needed for inner rotors
    """
    # rotor radius toward airgap
    Rg = Dg/2
    # mechanical slot angle in radians
    alpha_r = alpha_r*np.pi/180
    # first slot mechanical angular position in radians
    ang_pos = ang_pos*np.pi/180
    # compute rotor back iron height
    hrbi = (abs(Dre - Dg) - 2*hr)/2
    # add label in rotor iron (rotor lamination)
    # compute label coordinates
    if r == 1:     # outer stator, inner rotor
        R_l = Dre/2 + hrbi/2
    elif r == -1:  # inner stator, outer rotor
        R_l = Rg + hr + hrbi/2
    else:
        print('Input error!')
        return []
    x_l = R_l*np.cos(ang_pos)
    y_l = R_l*np.sin(ang_pos)
    lab_r = draw.add_label(x_l, y_l, 100, msh_areas[0])

    # add in-slots labels
    lab_slot = []  # labels in conducting regions inside slots
    air_slot = []  # labels in insulating/air regions in slots (i.e. openings)
    # labels in conducting regions are added only if slots are not empty
    if hwed != hr-hro:
        R_slot_cond = Rg - r*(hro/2+hwed/2+hr/2)
        for qq in range(0, Qr):
            lab_slot.append(draw.add_label(R_slot_cond*np.cos(ang_pos+qq*alpha_r),
                                           R_slot_cond*np.sin(ang_pos+qq*alpha_r),
                                           500+qq+1, msh_areas[1]))
    # labels in slot openings are added only if slot are (half) open
    if wro != 0 and hro+hwed != 0:
        R_slot_air = Rg - r*(hro+hwed)/2
        for qq in range(0, Qr):
            air_slot.append(draw.add_label(R_slot_air*np.cos(ang_pos+qq*alpha_r),
                                           R_slot_air*np.sin(ang_pos+qq*alpha_r),
                                           99, msh_areas[-1]))

    # impose physical tags on stator boundary lines
    for ll in Az0_lns:
        ll.ph_tag = 10

    # assign No-Mesh region for inner stators with shaft/passing clamp diameter
    hol = []
    if r == 1 and Dre != 0:
        hol.append(draw.add_hole(0, 0))
    # output label list
    return [lab_r] + lab_slot + air_slot + hol

# ==============================================================================
# Function to assign labels and physical tags to a rotor of a radial flux SPM machine
# using dolomites-tiziano for GetDP magnetostatic formulation
def ph_SPMrotor(draw, r, Az0_lns, Dg, Dre, p, alpha_p, alpha_m, tm, msh_areas, alpha_tilt=0, tin=0, ang_pos=0, typ_pm='radial', typ_rot='norm'):
    """Assign labels and physical tags to rotor of a radial-flux SPM machine.

    Labels and physical tags are imposed according to GetDP magnetostatic formulation.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to draw the pole profile.
    r : int
         1 --> airgap outside (inner rotor --> inrunner SPM machine)
        -1 --> airgap inside (outer rotor --> outrunner SPM machine)
    Dg : int or float
        Rotor diameter toward airgap [Dg > 0].
    Dre : int or float
        Rotor external diameter [Dre >= 0].
        "External" means not toward airgap.
        External diameter for outer rotors.
        Shaft and/or passing clamp diameter for inner rotors.
        It can be also set to zero if one does not want to close the rotor (outer rotors)
        or when shaft/passing clamp is not necessary in the drawing (inner rotors).
    p : int
        Number of pole pairs in the drawing.
    alpha_p: int or float
        Pole mechnical angle in deg [0 < alpha_p < 180].
    alpha_m : int or float
        Mechanical angle in deg covered by magnet at the airgap. [0 < alpha_m < alpha_p].
    tm : int or float
        Magnet thickness [tm >= 0].
    msh_areas : list of ints or floats
        Maximum areas of mesh triangles in rotor iron, PMs and air regions.
        [rot_iron_area, rot_pms_area, rot_air_area]
    alpha_tilt : int or float, default 0
        Magnet tilting angle in deg [0 <= alpha_tilt <= (alpha_p - alpha_m)/2].
    tin : int or float, default 0
        Magnet thickness buried inside iron [0 <= tin <= tm].
        0  --> magnet all in the airgap (default)
        tm --> magnet totally buried inside rotor iron
    ang_pos : int or float, default 0
        Mechanical angular position with respect to origin in deg which in the
        first pole was drawn.
    typ_pm : str, default 'radial'
        Magnet block type, 'radial' (deafult) or 'square'.
    typ_rot : str, default 'norm'
        typ_rot = 'norm' --> normal or tilted asymmetric topologies
        typ = 'cons' --> consequent pole topologies

    Returns
    -------
    [lab_r] + lab_pm + air_rot + hol : list of tiziano labels and holes
        Tiziano labels and holes added to the drawing for SPM rotor meshing.
        lab_r    -> rotor iron, lamination   (ph=100)
        lab_pm   -> rotor magnets            (ph=101,..., 100+2*p)
        air_slot -> rotor air near pm        (ph=99)
        hol      -> No-Mesh central region, if needed for inner rotors
    """
    # rotor radius toward airgap
    Rg = Dg/2
    # mechanical pole angle in radians
    alpha_p = alpha_p*np.pi/180
    # mechanical angle covered by magnet in radians
    alpha_m = alpha_m*np.pi/180
    # magnet tilt angle in radians
    alpha_tilt = alpha_tilt*np.pi/180
    # first pole mechanical angular position in radians
    ang_pos = ang_pos*np.pi/180
    # compute magnet tickness in the air gap
    tout = tm - tin
    # add label in rotor iron (rotor lamination)
    # compute label coordinates
    if r == 1:     # outer stator, inner rotor
        R_l = Dre/2 + (Rg - Dre/2 - tm)/2
    elif r == -1:  # inner stator, outer rotor
        R_l = Rg + tm + (Dre/2 - Rg - tm)/2
    else:
        print('Input error!')
        return []
    x_l = R_l*np.cos(ang_pos)
    y_l = R_l*np.sin(ang_pos)
    lab_r = draw.add_label(x_l, y_l, 100, msh_areas[0])
    # add Permanent Magnets labels
    lab_pm = []   # labels in surface mounted permanent magnets
    air_rot = []  # labels in insulating/air regions near magnet inside airgap
    if tm != 0:
        if typ_pm == 'radial':
            rr = 0
        elif typ_pm == 'square':
            rr = Rg*(1-np.cos(alpha_m/2))
        else:
            print('Input error!')
            return []
        R_pm = Rg - r*tm/2 - rr
        if typ_rot == 'norm':
            step_rot = 1
            inv_pm = 1
            start = 0
        elif typ_rot == 'cons':
            step_rot = 2
            inv_pm = -1
            start = 1
        else:
            print('Input error!')
            return []
        # loop to assign pm labels
        t = 1
        for pp in range(0, int(2*p), step_rot):
            lab_pm.append(draw.add_label(R_pm*np.cos(ang_pos+pp*alpha_p+t*alpha_tilt),
                                         R_pm*np.sin(ang_pos+pp*alpha_p+t*alpha_tilt),
                                         100+pp+1, msh_areas[1]))
            t = -inv_pm*t

    if tout != 0:  # add labels in air regions near magnets
        if typ_rot == 'cons' or alpha_p > alpha_m:
            R_air = Rg - r*tout/2
            alpha_a = (alpha_p - alpha_m)/2 - alpha_tilt
            if typ_rot == 'cons' and alpha_p <= alpha_m:
                start_ang = alpha_tilt
                inv = 0
            elif -1e-10 < alpha_a < 1e-10:
                start_ang = -(alpha_p+alpha_m)/4 + alpha_tilt/2
                inv = -1
            else:
                start_ang = (alpha_p+alpha_m)/4 + alpha_tilt/2
                inv = 1
            tt = 1
            for pp in range(start, int(2*p), step_rot):
                air_rot.append(draw.add_label(R_air*np.cos(ang_pos+tt*start_ang+pp*alpha_p),
                                              R_air*np.sin(ang_pos+tt*start_ang+pp*alpha_p),
                                              99, msh_areas[-1]))
                tt = inv*tt
            # add another air rotor label --> needed in case of periodicity
            if alpha_p - alpha_m - alpha_a > 1e-10:
                air_rot.append(draw.add_label(R_air*np.cos(ang_pos-(alpha_p+alpha_m)/4 + alpha_tilt/2),
                                              R_air*np.sin(ang_pos-(alpha_p+alpha_m)/4 + alpha_tilt/2),
                                              99, msh_areas[-1]))

            # add others air regions in fron of pm in case of inner rotor with square magnets
            if typ_pm == 'square' and r == 1:
                R_air = Rg - rr/2
                t = 1
                for pp in range(0, int(2*p), step_rot):
                    air_rot.append(draw.add_label(R_air*np.cos(ang_pos+pp*alpha_p+t*alpha_tilt),
                                                  R_air*np.sin(ang_pos+pp*alpha_p+t*alpha_tilt),
                                                  99, msh_areas[-1]))
                    t = -inv_pm*t

    # impose physical tags on rotor boundary lines
    for ll in Az0_lns:
        ll.ph_tag = 10

    # assign No-Mesh region for inner rotors with shaft diameter
    hol = []
    if r == 1 and Dre != 0:
        hol.append(draw.add_hole(0, 0))
    # output label list
    return [lab_r] + lab_pm + air_rot + hol

def airgap(draw, g, Dint, da=1.5, per=1, ang_pos=0, ng=5, coherence=False):
    """Divide airgap into 5 or 3 equal layers according to GetDP synchronous machines formulation.

    5 or 3 equal layers in airgap are considered with Moving Band in the middle.
    This function does not impose any label nor physical tag on lines and points.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to draw the airgap layers.
    g : int or float
         Airgap thickness.
    Dint : int or float
        Machine internal diameter toward gap (inner rotor or inner stator) [Dint > 0].
    da : int or float, default 1.5
        Arcs discretization angle in deg [da > 0].
    per : int, default 1
        Machine periodicity in the drawing.
        Default: no periodicity, entire machine
    ang_pos: int or float, default 0
        Angular position of airgap arcs starting points in deg.
    ng : int, default 5
        Number of equal layers to divide airgap. Only 5 or 3 allowed!
    coherence : bool, default False
        True  --> check coherence when tiziano drawing methods are called
        False --> do not check coherence when tiziano drawing methods are called

    Returns
    -------
    pts : list of tiziano points
        Airgap arcs discretization points.
    lns : list of tiziano lines
        Airgap arcs discretization lines.
    lns_int : list of tiziano lines.
        Internal Moving Band arc lines.
    lns_out : list of tiziano lines
        External Moving Band arc lines.
    """
    # check total number of airgap layers
    if ng not in [3, 5]:
        print('machine_drawing ==> Invalid input!')
        return [], [], [], []
    # compute machine internal radius
    # (inner stator or inner rotor radius toward airgap)
    Rint = Dint/2
    # arcs discretization angle in radians
    da = da*np.pi/180
    # angular position of airgap arcs startin points in radians
    ang_pos = ang_pos*np.pi/180
    # arcs spanning angle in radians
    alpha_sim = 2*np.pi/per
    # compute equal airgap layers thickness
    t = g/ng
    pts = []
    lns = []
    if per == 1:   # no periodicity ==> entire machine
        for tt in range(1, ng):
            pts.append(draw.add_point(Rint + tt*t, 0, coherence=coherence))
            pts.append(draw.add_point(-(Rint + tt*t), 0, coherence=coherence))
        for ii in range(0, 2*(ng-1), 2):
            arc_pts1, arc_lns1 = draw.add_arc(pts[ii], pts[ii+1], np.pi, da, coherence=coherence)
            arc_pts2, arc_lns2 = draw.add_arc(pts[ii], pts[ii+1], -np.pi, da, coherence=coherence)
            pts.extend(arc_pts1[1:-1] + arc_pts2[1:-1])
            lns.extend(arc_lns1 + arc_lns2)
    else:          # machine periodicity exploited
        for tt in range(1, ng):
            pts.append(draw.add_point((Rint + tt*t)*np.cos(ang_pos),
                                      (Rint + tt*t)*np.sin(ang_pos), coherence=coherence))
            pts.append(draw.add_point((Rint + tt*t)*np.cos(ang_pos+alpha_sim),
                                      (Rint + tt*t)*np.sin(ang_pos+alpha_sim), coherence=coherence))
        for ii in range(0, 2*(ng-1), 2):
            arc_pts1, arc_lns1 = draw.add_arc(pts[ii], pts[ii+1], alpha_sim, da, coherence=coherence)
            pts.extend(arc_pts1[1:-1])
            lns.extend(arc_lns1)

    # take internal arc lines for next Moving Band physical tags imposition
    # valid for n=5 or n=3 airgap layers only!!!
    if ng == 5:
        # internal MB side lines
        lns_int = lns[int(len(lns)/4):int(len(lns)/4*2)]
        # outer MB side lines
        lns_out = lns[int(len(lns)/4*2):int(len(lns)/4*3)]
    elif ng == 3:
        lns_int = lns[:int(len(lns)/2)]
        lns_out = lns[int(len(lns)/2):]

    return pts, lns, lns_int, lns_out

def airgap_oneside(draw, r, g, Dg, side='stat', da=1.5, per=1, ang_pos=0, ng=5, coherence=False):
    """Discretize stator or rotor airgap side into 2 or 1 equal layer(s).

    5 or 3 equal layers in airgap are considered with Moving Band in the middle.
    This function does not impose any label nor physical tag on lines and points.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to draw the airgap layers.
    r : int
         1 --> inner rotor, outer stator --> inrunner machine
        -1 --> outer rotor, inner stator --> outrunner machine
    g : int or float
         Airgap thickness.
    Dg : int or float
         Stator or rotor diameter toward gap [Dg > 0].
    side : str, default 'stat'
        'stat'  --> draw arcs of stator side airgap layers
        'rot'   --> draw arcs of rotor side airgap layers
    da : int or float, default 1.5
        Arcs discretization angle in deg [da > 0].
    per : int, default 1
        Machine periodicity in the drawing.
        Deafult: no periodicity, entire machine
    ang_pos: int or float, default 0
        Angular position of airgap arcs starting points in deg.
    ng : int, default 5
        Number of equal layers to divide airgap. Only 5 or 3 allowed!
    coherence : bool, default False
        True  --> check coherence when tiziano drawing methods are called
        False --> do not check coherence when tiziano drawing methods are called

    Returns
    -------
    pts : list of tiziano points
        Airgap arcs discretization points.
    lns : list of tiziano lines
        Airgap arcs discretization lines.
    lns_mb : list of tiziano lines.
        Moving Band arc lines.
    """
    # check total number of airgap layers
    if ng not in [3, 5]:
        print('machine_drawing ==> Invalid input!')
        return [], [], []
    # set airgap layers orientation with respect to rotor or stator
    if side == 'stat':
        rr = -r
    elif side == 'rot':
        rr = r
    else:
        print('machine_drawing ==> Invalid input!')
        return [], [], []
    # number of airgap layers on stator or rotor side
    n = int((ng-1)/2)  # ng=5 ==> n=2 ;  ng=3 ==> n=1 ;
    # compute stator or rotor radius toward airgap
    Rg = Dg/2
    # arcs discretization angle in radians
    da = da*np.pi/180
    # angular position of airgap arcs startin points in radians
    ang_pos = ang_pos*np.pi/180
    # arcs spanning angle in radians
    alpha_sim = 2*np.pi/per
    # compute equal airgap layers thickness
    t = g/ng
    pts = []
    lns = []
    if per == 1:   # no periodicity ==> entire machine
        for tt in range(1, n+1):
            pts.append(draw.add_point(Rg + rr*tt*t, 0, coherence=coherence))
            pts.append(draw.add_point(-(Rg + rr*tt*t), 0, coherence=coherence))
        for ii in range(0, 2*n-1, 2):
            arc_pts1, arc_lns1 = draw.add_arc(pts[ii], pts[ii+1], np.pi, da, coherence=coherence)
            arc_pts2, arc_lns2 = draw.add_arc(pts[ii], pts[ii+1], -np.pi, da, coherence=coherence)
            pts.extend(arc_pts1[1:-1] + arc_pts2[1:-1])
            lns.extend(arc_lns1 + arc_lns2)
    else:          # machine periodicity exploited
        for tt in range(1, n+1):
            pts.append(draw.add_point((Rg + rr*tt*t)*np.cos(ang_pos),
                                      (Rg + rr*tt*t)*np.sin(ang_pos), coherence=coherence))
            pts.append(draw.add_point((Rg + rr*tt*t)*np.cos(ang_pos+alpha_sim),
                                      (Rg + rr*tt*t)*np.sin(ang_pos+alpha_sim), coherence=coherence))
        for ii in range(0, n+1, 2):
            arc_pts1, arc_lns1 = draw.add_arc(pts[ii], pts[ii+1], alpha_sim, da, coherence=coherence)
            pts.extend(arc_pts1[1:-1])
            lns.extend(arc_lns1)

    # take arc lines on boundary for next Moving Band physical tags imposition
    # valid for n=2 or n=1 airgap layers on stator/rotor side only!!!
    # MB side lines
    if n == 2:  # ng = 5; n =2
        lns_mb = lns[int(len(lns)/2):]
    elif n == 1:  # ng = 3; n = 1
        lns_mb = lns

    return pts, lns, lns_mb

def ph_airgap(draw, lns_int, lns_out, r, g, Dint, msh_area, ang_pos=0, ng=5):
    """Assign labels and pysical tags to airgap layers and lines according to GetDP formulation.

    5 or 3 equal layers in airgap are considered with Moving Band in the middle.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to put airgap labels and ph tags.
    lns_int : list of tiziano lines
        Lines that discretize internal Moving Band arc.
    lns_out : list of tiziano lines
        Lines that discretize outer Moving Band arc.
    r : int
         1 --> inner rotor, outer stator --> inrunner machine
        -1 --> outer rotor, inner stator --> outrunner machine
    g : int or float
         Airgap thickness.
    Dint : int or float
         Machine internal diameter toward gap [Dint > 0].
         Rotor diamter toward airgap for inrunner machines.
         Stator diamter toward airgap for outrunner machines.
    msh_area : int or float
        Maximum area allowed for mesh triangles in airgap layers.
    ang_pos: int or float, default 0
        Angular position which in to place airgap layers labels in deg.
        Set this parameter carefully if you want to exploit machine periodicity.
    ng : int, default 5
        Number of equal layers to divide airgap. Only 5 or 3 allowed!

    Returns
    -------
    lab : list of tiziano labels and holes
        Airgap labels and the hole in the central Moving Band layer.
    """
    # check total number of airgap layers
    if ng not in [3, 5]:
        print('machine_drawing ==> Invalid input!')
        return []
    # compute equal airgap layers thickness
    t = g/ng
    # compute machine internal radius
    # (inner stator or inner rotor radius toward airgap)
    Rint = Dint/2 + t/2
    # airgap labels angular position with respect to origin in radians
    ang_pos = ang_pos*np.pi/180
    cos = np.cos(ang_pos)
    sin = np.sin(ang_pos)
    if r == 1:  # inner rotor, outer stator
        ph_lns = [30, 20]
        if ng == 3:
            ph_lab = [98, 998]
        if ng == 5:
            ph_lab = [98, 97, 997, 998]
    elif r == -1:  # outer rotor, inner stator
        ph_lns = [20, 30]
        if ng == 3:
            ph_lab = [998, 98]
        if ng == 5:
            ph_lab = [998, 997, 97, 98]
    else:
        print('Input error.')
        return []
    lab = []  # list for labels and holes
    if ng == 3:
        lab.append(draw.add_label(Rint*cos, Rint*sin, ph_lab[0], msh_area))
        lab.append(draw.add_label((Rint+2*t)*cos, (Rint+2*t)*sin, ph_lab[1], msh_area))
        lab.append(draw.add_hole((Rint+t)*cos, (Rint+t)*sin))
    elif ng == 5:
        lab.append(draw.add_label(Rint*cos, Rint*sin, ph_lab[0], msh_area))
        lab.append(draw.add_label((Rint+t)*cos, (Rint+t)*sin, ph_lab[1], msh_area))
        lab.append(draw.add_label((Rint+3*t)*cos, (Rint+3*t)*sin, ph_lab[2], msh_area))
        lab.append(draw.add_label((Rint+4*t)*cos, (Rint+4*t)*sin, ph_lab[3], msh_area))
        lab.append(draw.add_hole((Rint+2*t)*cos, (Rint+2*t)*sin))
     # assign Moving Band lines physical tags
    for ll in lns_int:
        ll.ph_tag = ph_lns[0]
    for ll in lns_out:
        ll.ph_tag = ph_lns[1]

    return lab

def ph_airgap_oneside(draw, lns_mb, r, g, Dg, msh_area, side='stat', ang_pos=0, ng=5):
    """Assign labels and pysical tags to airgap layers and lines on stator or rotor side according to GetDP formulation.

    2 or 1 equal layers in airgap are considered on both rotor and stator side with Moving Band in the middle.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in to put airgap labels and ph tags.
    lns_mb : list of tiziano lines
        Lines that discretize a Moving Band arc.
    r : int
         1 --> inner rotor, outer stator --> inrunner machine
        -1 --> outer rotor, inner stator --> outrunner machine
    g : int or float
         Airgap thickness.
    Dg : int or float
         Rotor or stator diameter toward gap [Dg > 0].
    msh_area : int or float
        Maximum area allowed for mesh triangles in airgap layers.
    side : str, default 'stat'
        'stat'  --> assign labels and ph tags to stator side airgap layers
        'rot'   --> assign labels and ph tags to rotor side airgap layers
    ang_pos: int or float, default 0
        Angular position which in to place airgap layers labels in deg.
        Set this parameter carefully if you want to exploit machine periodicity.
    ng : int, default 5
        Number of equal layers to divide airgap. Only 5 or 3 allowed!

    Returns
    -------
    lab : list of tiziano labels
        Stator or rotor side airgap layers labels.
    """
    # check total number of airgap layers
    if ng not in [3, 5]:
        print('machine_drawing ==> Invalid input!')
        return []
    # set airgap layers orientation with respect to rotor or stator
    if side == 'stat':
        rr = -r
    elif side == 'rot':
        rr = r
    else:
        print('machine_drawing ==> Invalid input!')
        return []
    # number of airgap layers on stator/rotor side
    n = int((ng-1)/2)  # ng=5 ==> n=2 ;  ng=3 ==> n=1 ;
    # compute equal airgap layers thickness
    t = g/ng
    # compute stator or rotor radius toward airgap
    Rg = Dg/2
    # airgap labels angular position with respect to origin in radians
    ang_pos = ang_pos*np.pi/180
    cos = np.cos(ang_pos)
    sin = np.sin(ang_pos)
    if side == 'stat':    # stator side airgap layers labels physical tags
        ph_lab = [998, 997]
        ph_lns = 20
    elif side == 'rot':   # rotor side airgap layers labels physical tags
        ph_lab = [98, 97]
        ph_lns = 30
    else:
        print('Input error.')
        return []
    lab = []  # list for labels and holes
    lab.append(draw.add_label((Rg+rr*t/2)*cos, (Rg+rr*t/2)*sin, ph_lab[0], msh_area))
    if ng == 5:
        lab.append(draw.add_label((Rg+rr*(t/2+t))*cos, (Rg+rr*(t/2+t))*sin, ph_lab[1], msh_area))
    # assign Moving Band lines physical tags
    for ll in lns_mb:
        ll.ph_tag = ph_lns

    return lab

def close_per(draw, Dg, Dext, g, ph_tags, N, per, tout=0, N_air=None, ang_pos=0, ng=5, coherence=False):
    """Close stator or rotor to set machine periodicity boundary conditions.

    This function draws and discretizes two lines (i.e. master and slave lines
    in GetDP formulation) to close the given machine slice and assign to them
    the given physical tags.

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in the function will take effect.
    Dg : int or float
        Stator or rotor diamter toward airgap [Dg > 0].
    Dext : int or float
        Stator or rotor external diameter [Dext > 0].
    g : int or float
        Airgap thickness.
        Set g=0 if you do not want to close airgap layers.
    ph_tags : list of int
        [ph_tag for 'rigth' closing line, ph_tag for 'left' closing line]
        Physical tags for master and slave lines boundary lines.
    N : int
        Number of discretization points for lamination closing lines [N >= 0].
        (Not in airgap!!! Only stator or rotor lamination boundaries are discretized).
        Set N=0 if do not want to discretize lamination closing lines.
    per : int
        Machine periodicity in the drawing.
    tout : int or float, default 0
        Air thickness outside lamination (SPM rotor case - air regions near magnets).
        Set it to zero if you do not have air regions between laminai and airgap
        layers to close (e.g. when closing a stator).
    N_air : int, default None
        Number of discretization points for air regions closing lines [N_air >= 0].
        (E.g. for SPM rotors air regions; not airgap layers !!!).
        Set N_air=0 if do not want to discretize air regions closing lines.
    ang_pos : int or float, default 0
        Angular position which in arc toward airgap start in deg.
    ng : int, default 5
        Number of airgap equal discretization layers. Must be and odd integer.
    coherence : bool, default False
        True  --> check coherence when tiziano drawing methods are called
        False --> do not check coherence when tiziano drawing methods are called

    Returns
    -------
    lns _right : list fo tiziano lines
        Right closing lines.
    lns _left : list fo tiziano lines
        Left closing lines.
    """
    # radius toward airgap
    Rg = Dg/2
    # external radius of lamination
    Rext_lam = Dext/2
    if Rg > Rext_lam:  # inner rotor or inner stator
        rr = -1
    else:              # outer rotor or outer stator
        rr = 1
    # radius toward airgap of lamination
    Rg_lam = Rg + rr*tout
    # radius toward airgap (plus all airgap layers!!!)
    Rg_plus = Rg - rr*((ng-1)/2)*g/ng
    # radius between stator or rotor side airgap layers
    Rg_plus_1 = Rg - rr*((ng-1)/2-1)*g/ng
    # starting angular position in radians
    ang_pos = ang_pos*np.pi/180
    cos = np.cos(ang_pos)
    sin = np.sin(ang_pos)
    # mahine angular portion in radians
    alpha_sim = 2*np.pi/per
    # add points to draw 'right' closing line
    p1 = draw.add_point(Rext_lam*cos, Rext_lam*sin, coherence=coherence)
    p2 = draw.add_point(Rg_lam*cos, Rg_lam*sin, coherence=coherence)
    p3 = draw.add_point(Rg*cos, Rg*sin)
    p4 = draw.add_point(Rg_plus_1*cos, Rg_plus_1*sin, coherence=coherence)
    p5 = draw.add_point(Rg_plus*cos, Rg_plus*sin, coherence=coherence)
    # add right line to close lamination
    l_right_lam = draw.add_line(p1, p2, ph=ph_tags[0], coherence=coherence)
    # discretize lamination right closing line with N points, if needed
    # coherence always checked here!!!
    dd = abs(Rg_lam - Rext_lam)/(N+1)
    lns_lam = []
    if N != 0:
        for nn in range(1, N+1):
            _, new_lns = draw.add_point(p1.x-rr*nn*dd*cos, p1.y-rr*nn*dd*sin, get_lns=True)
            if nn != N:
                lns_lam.append(new_lns[0])
            else:
                lns_lam.extend(new_lns)
    else:
        lns_lam = l_right_lam
    # discretize air region right closing line with N_air points, if needed
    # coherence always checked here!!!
    lns_air = []
    if tout != 0:
        # add right line to close air regions (if any)
        l_right_air = draw.add_line(p2, p3, ph=ph_tags[0], coherence=coherence)
        if N_air is not None:
            dd = abs(Rg_lam - Rg)/(N_air+1)
            for nn in range(1, N_air+1):
                _, new_lns = draw.add_point(p2.x-rr*nn*dd*cos, p2.y-rr*nn*dd*sin, get_lns=True)
                if nn != N_air:
                    lns_air.append(new_lns[0])
                else:
                    lns_air.extend(new_lns)
        else:
            lns_air = l_right_air

    # add right line to close airgap layers (coherence always True here!!!)
    l_right_gap_1 = draw.add_line(p3, p4, ph=ph_tags[0])
    # no ph tag on the last gap layer closing line !!!
    l_right_gap   = draw.add_line(p4, p5, coherence=coherence)

    lns_right = lns_lam + lns_air + l_right_gap_1 + l_right_gap  # right closing lines
    # copy rotate right closing lines to get left closing lines
    _, lns_left, _, _ = draw.copy_rotate(lns_right, alpha_sim, entity=1, coherence=coherence)
    # assign new ph tag to left closing lines
    for ll in lns_left[:-1]:
        ll.ph_tag = ph_tags[1]

    return lns_right, lns_left

def cloe_MB_per(draw, per, mb_ph=30):
    """Add nodes and edges to triangle mesh to close properly periodic Moving Band for GetDP formulation.

    This function makes changes only on mesh entities' lists of input tiziano
    drawing instance.
    It does not return anything (no output).

    Parameters
    ----------
    draw : tiziano drawing object
        Drawing object from dolomites-tiziano which in the function will take effect.
    per : int
        Machine periodicity in the drawing [per>1].
    mb_ph : int, default 30
        Rotor (or stator) side Moving Band physical tag.
        Other MB tags will be: mb_ph+1, mb_ph+2, ..., mb_ph+per-1
        Default: 31, 32, 33, ..., 30+per-1
    """
    if len(draw.triangle_edges) == 0:
        print('Invalid input! Please call triangle first.')
        return
    alpha_sim = 2*np.pi/per  # rotation angle (simulation angle) in rad
    node_idx = []
    node_xy = np.array([[], []])
    for ee in draw.triangle_edges:
        if ee[1] == mb_ph:
            node1_idx = ee[0][0]
            node2_idx = ee[0][1]
            if node1_idx not in node_idx:
                node_idx.append(node1_idx)
                node_xy = np.hstack((node_xy,
                                     np.array([[draw.triangle_points[node1_idx][0][0]],
                                               [draw.triangle_points[node1_idx][0][1]]])))
            if node2_idx not in node_idx:
                node_idx.append(node2_idx)
                node_xy = np.hstack((node_xy,
                                     np.array([[draw.triangle_points[node2_idx][0][0]],
                                               [draw.triangle_points[node2_idx][0][1]]])))
    node_phi = np.arctan2(node_xy[1, :], node_xy[0, :])
    sorted_idx = np.argsort(node_phi)
    for nn in range(1, per):
        offset = len(draw.triangle_points)
        rot_mat = np.array([[np.cos(nn*alpha_sim), -np.sin(nn*alpha_sim)],
                            [np.sin(nn*alpha_sim),  np.cos(nn*alpha_sim)]])
        new_node_xy = rot_mat @ node_xy
        idx_to_loop = sorted_idx[1:]
        if nn == per-1:
            idx_to_loop = sorted_idx[1:-1]
        for ii in idx_to_loop:
            new_node = [(new_node_xy[0, ii], new_node_xy[1, ii]), 0]
            draw.triangle_points.append(new_node)
        if nn == 1:
            draw.triangle_edges.append(((node_idx[sorted_idx[-1]], offset), mb_ph+nn))
            offset = offset + 1
            idx_to_loop = idx_to_loop[:-1]
        if nn == per-1:
            draw.triangle_edges.append(((node_idx[sorted_idx[0]], offset+len(idx_to_loop)-1), mb_ph+nn))
        for ii, _ in enumerate(idx_to_loop):
            draw.triangle_edges.append(((ii+offset-1, ii+offset), mb_ph+nn))

def calc_w(D, Q, w_in):
    """Compute slot width or tooth width at a given diameter (center 0,0) with a given number of slots."""
    w_out = D*np.sin(np.pi/Q - np.arcsin(w_in/D))
    return w_out
