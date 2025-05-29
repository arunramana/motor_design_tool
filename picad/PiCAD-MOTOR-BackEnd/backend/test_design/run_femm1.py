import ezdxf
from random import random
from ezdxf.addons import r12writer
import numpy as np
from numpy import pi

import femm

import flask

from .trig import *

import pythoncom
import json

from dolomites import koil


class GenerateFEMM:

    def __init__(self, directory, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width, bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w):

        #file dxf generated
        self.directory = directory

        #origin
        self.xo = 0
        self.yo = 0

        self.no_slots = no_slots
        self.no_poles = no_poles

        self.stator_od = stator_od
        self.stator_id = stator_id
        self.rotor_od = rotor_od
        self.shaft_dia = shaft_dia
        self.slot_opening = slot_opening
        self.tooth_tip_angle = tooth_tip_angle
        self.tooth_width = tooth_width
        self.yoke_thickness = yoke_thickness
        self.pole_length = pole_length
        self.pole_width = pole_width

        self.bridge_thickness = bridge_thickness

        #duct values are in 3 arrays below and must have the same dimension
        self.duct_distances = duct_distances
        self.duct_radii = duct_radii
        self.duct_angles = duct_angles

        self.notch_depth = notch_depth
        self.notch_angle = notch_angle

        self.stack_depth = stack_depth

        #circuits
        self.no_turns = no_turns
        self.u = u
        self.v = v
        self.w = w
        #-----------------------------calculated variables ---------------------
        self.has_notch = notch_angle!=0 or notch_depth!=0


        self.stator_or = stator_od/2
        self.stator_ir = stator_id/2
        self.rotor_or = rotor_od/2
        self.shaft_rad = shaft_dia/2

        self.air_gap = stator_id - rotor_od

        self.theta_inc_stator = 360/no_slots

        self.stator_outter_circum = pi*stator_od
        self.stator_inner_circum = pi*stator_id

        self.stator_opening_total = self.stator_inner_circum - self.slot_opening*no_slots
        #s = r*theta
        self.tooth_arc_s = self.stator_opening_total/self.no_slots
        self.tooth_arc_theta = (self.tooth_arc_s*2)/self.stator_id

        self.tooth_arc_theta_deg = ((self.tooth_arc_theta*180)/pi)/2 #divide by 2 for start and end angle

        self.rotor_outter_circum = pi*self.rotor_od

        self.max_pole_length = self.rotor_outter_circum/self.no_poles

        #pole outter points distance from origin
        self.pole_dist = (self.rotor_od-2*self.bridge_thickness)/2
        self.pole_dist = self.rotor_or-self.bridge_thickness
        self.pole_angle = (self.pole_length/self.pole_dist)*(180/pi)

        self.pole_dist = (self.pole_dist**2 - (self.pole_length/2)**2)**0.5

        b = 360/self.no_poles

        a = 2*taninv(self.pole_length/(self.pole_dist*2))

        self.pole_arc = (a/b)*180




        #multiply by 2 since pole arc is 180 deg
        self.pole_angle = (self.pole_arc*2)/self.no_poles
        self.pole_gap = ((180-self.pole_arc)*2)/self.no_poles

        self.theta_inc_rotor = 360/self.no_poles



    def run_femm(self,verbose=True,saveplots=False,fast_run=False):

        pythoncom.CoInitialize()

        #open femm
        femm.openfemm()

        femm.main_minimize()
        doctype = 0 #for magnetic problem
        femm.newdocument(doctype)

        #define problem
        self.problem_definition()

        #add materials
        self.add_material()

        #add circuits u,v,w
        self.add_circuit()

        #import dxf
        self.read_dxf()

        #group rotor and stator
        self.group_stator_and_rotor()

        #add properties to various components
        self.add_properties()

        #add boundary
        self.add_boundary()

        #for line integral after analyze
        line_integral_coords = self.draw_line_integral()


        #for yoke & tooth B_n
        tooth_x,tooth_y,yoke_x,yoke_y = self.draw_yoke_and_tooth_points()


        #save file
        #fname = self.filename.replace(".dxf",".FEM")
        fname = self.directory+"/motor.FEM"
        femm.mi_saveas(fname)

        #turn off smart mesh to increase speed
        femm.mi_smartmesh(0)
        #mesh
        #femm.mi_createmesh()


        #run analysis
        '''
        The flag parameter controls whether the solver
        window is visible or minimized. For a visible window, specify 0. For a minimized window, flag
        should be set to 1.
        '''

        flag=1
        femm.mi_analyze(flag)


        #load solution
        femm.mi_loadsolution()



        if(saveplots and not fast_run):
            femm.main_maximize()

            #show density plot
            self.show_density_plot()

            femm.mo_zoom(-112,-56,112,56)
            f = fname.replace(".FEM","")
            f = f+"_density_plot.png"

            femm.mo_savebitmap(f)

            femm.main_minimize()

        #line integral to obtain Bavg



        self.select_line_integral(line_integral_coords)

        self.airgap_flux_density,self.airgap_flux = self.line_integral_results(verbose=verbose)

        #make plot

        if(not fast_run):
            if(saveplots):
                f = fname.replace(".FEM","")
                f = f+"_plot.txt"
                self.make_plot(Filename=f)
            else:
                self.make_plot()

        #clear contour
        femm.mo_clearcontour()


        #get pole line integral
        self.select_pole_line_integral()

        self.flux_leakage_factor = self.flux_leakage(self.airgap_flux,verbose=verbose)

        #clear contour
        femm.mo_clearcontour()

        #get b from tooth
        femm.mo_selectpoint(tooth_x,tooth_y)
        bx, by = femm.mo_getb(tooth_x,tooth_y)
        self.b_tooth = np.sqrt(bx**2+by**2)
        #clear contour
        femm.mo_clearcontour()


        #get b from yoke
        femm.mo_selectpoint(yoke_x,yoke_y)
        bx, by = femm.mo_getb(yoke_x,yoke_y)
        self.b_yoke = np.sqrt(bx**2+by**2)
        #clear contour
        femm.mo_clearcontour()

        self.b_yoke_to_b_tooth_ratio = self.b_yoke/self.b_tooth

        if(verbose):
            print("\nb_yoke : {}\nb_tooth: {}\nb_yoke/b_tooth: {}\n".format(self.b_yoke,self.b_tooth,self.b_yoke_to_b_tooth_ratio))



        femm.closefemm()

    def problem_definition(self):

        #problem definition

        freq=0
        units="millimeters" # ’inches’, ’millimeters’,’centimeters’, ’mils’, ’meters’, and ’micrometers’
        typ="planar" # ’planar’ 'axi'
        precision = '1E-8'
        depth = self.stack_depth
        minangle = 30
        acsolver = 0  # 0 for successive approximation, 1 for Newton.

        femm.mi_probdef(freq,units,typ,precision,depth,minangle,(acsolver))


    def add_material(self):

        #adding to material library

        femm.mi_getmaterial('Air')
        femm.mi_getmaterial('20 SWG')
        femm.mi_getmaterial('N42')

        femm.mi_getmaterial('M250-35A')

        propname="A=0"

        A0, A1, A2, Phi, Mu, Sig, c0, c1, BdryFormat, ia, oa = 0,0,0,0,0,0,0,0,0,0,0

        femm.mi_addboundprop(propname, A0, A1, A2, Phi, Mu, Sig, c0, c1, BdryFormat, ia, oa)
        #add steel
        #femm.mi_addmaterial('M250-35A', 1, 1, 0, 0, 0, 0.35, 0, 0.98, 0, 0, 0, 0, 0)


    def add_circuit(self):

        #adding properties circuit

        circuittype = 1 #0- parallel 1- series

        u=self.u#i current in amps
        v=self.v
        w=self.w

        femm.mi_addcircprop('U', self.u, circuittype)

        femm.mi_addcircprop('V', self.v, circuittype)

        femm.mi_addcircprop('W', self.w, circuittype)


    def read_dxf(self):

        #import dxf generated and add origin at 0,0
        self.dxf_filename = self.directory+"/motor.dxf"
        femm.mi_readdxf(self.dxf_filename)

        #add origin

        femm.mi_addnode(self.xo,self.yo)


    def group_stator_and_rotor(self):

        #select circles for groups

        #select all
        R = self.stator_or*1.25
        editmode = 4

        femm.mi_selectcircle(self.xo,self.yo,R,editmode)
        femm.mi_setgroup(1)


        #select rotor
        R = (self.rotor_or+self.stator_ir)/2

        femm.mi_selectcircle(self.xo,self.yo,R,editmode)
        femm.mi_setgroup(2)


    def add_properties(self):

        #add all properties

        self.add_air_property()
        self.add_steel_property()
        self.add_magnet_property()
        self.add_circuit_property()

    def add_air_property(self):

        coords = []

        x,y = -self.shaft_rad/2,0
        femm.mi_addblocklabel(x,y) #shaft

        coords.append([x,y])

        x,y =  -(self.rotor_or+self.stator_ir)/2,0
        femm.mi_addblocklabel(x,y) #airgap

        coords.append([x,y])

        x,y = -(self.stator_or+5),0
        femm.mi_addblocklabel(x,y) #air outside

        coords.append([x,y])

        #air for rotor ducts

        if(self.duct_distances!=[0,0,0] and self.duct_radii != [0,0,0]):

            reference_angle = 0

            for iter_no in range(0,self.no_poles):

                if(iter_no%2==0):
                    reference_angle += self.theta_inc_rotor
                    continue


                for i in range(len(self.duct_distances)):

                    duct_dist = self.duct_distances[i]
                    duct_angle = self.duct_angles[i]
                    duct_radius = self.duct_radii[i]


                    x = duct_dist*sin(reference_angle+duct_angle)
                    y = duct_dist*cos(reference_angle+duct_angle)

                    femm.mi_addblocklabel(x,y)
                    coords.append([x,y])


                reference_angle += self.theta_inc_rotor





        #select all points in coords and give air as property

        blockname = 'Air'
        automesh = 1
        meshsize=0
        incircuit='None'
        magdir=0
        group=0
        turns=1


        for crds in coords:
            x,y = crds
            femm.mi_selectlabel(x,y)  #to select the block

        femm.mi_setblockprop(blockname, automesh, meshsize, incircuit, magdir, group, turns) #set block properties

        femm.mi_clearselected() #important to clear selected


    def add_steel_property(self):

        #add steel property

        coords = []

        if(self.duct_distances != [0,0,0]):

            x,y = -(min(self.duct_distances)+self.shaft_rad)/2,0

        else:

            x,y = -(self.rotor_or+self.shaft_rad)/2.5,0

        femm.mi_addblocklabel(x,y) #rotor -> group 2

        coords.append([x,y,2])

        x,y =  -(self.stator_or-self.yoke_thickness/2.25),0
        femm.mi_addblocklabel(x,y) #stator -> group 1

        coords.append([x,y,1])


        #select all points in coords and give steel as property

        blockname = 'M250-35A'
        automesh = 1
        meshsize=0
        incircuit='None'
        magdir=0
        turns=1


        for crds in coords:
            x,y,group = crds
            femm.mi_selectlabel(x,y)

            femm.mi_setblockprop(blockname, automesh, meshsize, incircuit, magdir, group, turns) #set block properties

            femm.mi_clearselected() #important to clear selected


    def add_magnet_property(self):

        #add magnets

        blockname = 'N42'
        automesh = 1
        meshsize=0
        incircuit='None'
        group=2
        turns=1



        theta = self.pole_angle/2

        r1 = self.pole_length/(sin(theta)*2)

        r1 = r1-(self.pole_width/2.1) #2.1 due to some discrepancy


        reference_angle = 90

        for no in range(self.no_poles):

            x = r1*cos(reference_angle)
            y = r1*sin(reference_angle)


            femm.mi_addblocklabel(x,y)
            femm.mi_selectlabel(x,y)

            magdir = reference_angle

            if(no%2!=0):
                magdir -= 180

            femm.mi_setblockprop(blockname, automesh, meshsize, incircuit, magdir, group, turns) #set block properties

            femm.mi_clearselected() #important to clear selected

            reference_angle += self.theta_inc_rotor


    def add_circuit_property(self):

        #add circuit

        blockname = '20 SWG'
        automesh = 1
        meshsize=0
        magdir = 0
        incircuit='None'
        group=1
        turns=1


        reference_angle = 0

        i=0

        winding_pattern = GenerateFEMM.getWindingPattern(self.no_slots,self.no_poles)

        for no in range(self.no_slots):

            #theta = (reference_angle-self.theta_inc_stator)/2

            theta = (reference_angle-self.theta_inc_stator)/2

            r = (self.stator_ir+(self.stator_or-self.yoke_thickness))/2

            #approx
            theta1 = ((self.tooth_width/r)*(180/pi))/2

            #x2 = r*sin(theta+theta1)
            #y2 = r*cos(theta+theta1)

            x2 = r*sin(theta+2)
            y2 = r*cos(theta+2)


            #x1 = r*sin(theta-theta1)
            #y1 = r*cos(theta-theta1)

            x1 = r*sin(theta-2)
            y1 = r*cos(theta-2)


            #--------------------------------------------------------

            femm.mi_addblocklabel(x1,y1)
            femm.mi_selectlabel(x1,y1)


            winding = winding_pattern[i]
            incircuit,turns = GenerateFEMM.getCircuitAndTurns(winding,self.no_turns)

            femm.mi_setblockprop(blockname, automesh, meshsize, incircuit, magdir, group, turns) #set block properties
            femm.mi_clearselected() #important to clear selected
            i+=1

            #--------------------------------------------------------

            femm.mi_addblocklabel(x2,y2)
            femm.mi_selectlabel(x2,y2)

            winding = winding_pattern[i]
            incircuit,turns = GenerateFEMM.getCircuitAndTurns(winding,self.no_turns)

            femm.mi_setblockprop(blockname, automesh, meshsize, incircuit, magdir, group, turns) #set block properties
            femm.mi_clearselected() #important to clear selected
            i+=1

            reference_angle += self.theta_inc_stator*2


    def draw_yoke_and_tooth_points(self):

        tooth_x, tooth_y = 0,(self.stator_or+self.stator_ir)/2

        r_yoke = self.stator_or-self.yoke_thickness/2

        theta_yoke = (360/self.no_slots)/2

        yoke_x = r_yoke*sin(theta_yoke)
        yoke_y = r_yoke*cos(theta_yoke)

        femm.mi_addnode(tooth_x,tooth_y)
        femm.mi_addnode(yoke_x,yoke_y)

        return [tooth_x,tooth_y,yoke_x,yoke_y]

    def draw_line_integral(self):

        #line integral for Bavg - draw points

        r = (self.rotor_or+self.stator_ir)/2


        line_integral_coords = []
        no_points = 5

        arc = self.pole_angle

        theta = -self.pole_angle/2 #start point


        for i in range(no_points+1):

            x = r*sin(theta)
            y = r*cos(theta)

            femm.mi_addnode(x,y)

            theta += arc/no_points

            line_integral_coords.append([x,y])


        return line_integral_coords

    def select_line_integral(self,line_integral_coords):

        #select points for line integral after drawing points for line integral

        for coord in line_integral_coords:

            x,y = coord

            femm.mo_selectpoint(x,y)

    def line_integral_results(self,verbose=True):
        #call this after adding and selecting the points

        li_type = 0  #0 is for Bn

        line_integral_results = femm.mo_lineintegral(li_type)

        normal_flux, Bn_avg = line_integral_results

        if(verbose):
            print("Air gap:\nNormal Flux : {} Wb\nB.n avg. : {} T".format(normal_flux,Bn_avg))

        #return bavg for leakage flux factor
        return Bn_avg,normal_flux


    def select_pole_line_integral(self):

        reference_angle = 0

        theta = self.pole_angle/2

        r1 = self.pole_length/(sin(theta)*2)

        theta1 = reference_angle+theta
        theta2 = reference_angle-theta


        r_pole = r1*cos(theta)
        h = r_pole-self.pole_width

        theta_t = taninv(self.pole_length/(2*h))

        r2 = (self.pole_length)/(2*sin(theta_t))


        theta3 = reference_angle+theta_t
        theta4 = reference_angle-theta_t

        x1 = r1*sin(theta1)
        y1 = r1*cos(theta1)

        x2 = r1*sin(theta2)
        y2 = r1*cos(theta2)


        x3 = r2*sin(theta3)
        y3 = r2*cos(theta3)

        x4 = r2*sin(theta4)
        y4 = r2*cos(theta4)


        femm.mo_selectpoint(x3,y3)
        femm.mo_selectpoint(x4,y4)




    def flux_leakage(self,flux,verbose=True):
        #returns the leakgae flux
        #do this after selecting pole line integral

        li_type = 0  #0 is for Bn

        line_integral_results = femm.mo_lineintegral(li_type)

        normal_flux2, Bn_avg2 = line_integral_results

        if(verbose):
            print("Pole:\nNormal Flux : {} Wb\nB.n avg. : {} T".format(normal_flux2,Bn_avg2))

        flux_leakage_factor = np.abs(flux/normal_flux2)

        if(verbose):
            print("\nFlux leakage Factor: ",flux_leakage_factor)

        return flux_leakage_factor


    def add_boundary(self):
        #adding boundary

        n=7 #u1-u7
        R=self.stator_or+10 #radius

        bc=0 #0-dirichlet, 1-neumann

        #femm.mi_makeABC(n,R,self.xo,self.yo,bc)

        #for stator boundary only
        x=0
        y1 = self.stator_or+5
        y2 = -self.stator_or-5

        maxsegdeg=5
        propname = "A=0"
        group =1
        hide=0

        femm.mi_selectarcsegment(x,y1)
        femm.mi_selectarcsegment(x,y2)


        femm.mi_setarcsegmentprop(maxsegdeg, propname, hide, group)
        femm.mi_clearselected()



    def show_density_plot(self):

        #density plot

        legend= 1
        gscale=0
        upper_B=2
        lower_B=0.002
        typ = "bmag"

        femm.mo_showdensityplot(legend,gscale,upper_B,lower_B,typ)


    def make_plot(self,Filename=None,FileFormat=0):

        #call after line integral calc.
        PlotType = 1 # |B|

        NumPoints = 1000

        if(Filename==None):
            femm.mo_makeplot(PlotType,NumPoints)
        else:
            femm.mo_makeplot(PlotType,NumPoints,Filename,FileFormat)







    @staticmethod
    def getWindingPattern(slots,poles):

        w = koil.m_phase_winding()
        m = 3   # number of phases assumed
        Q = slots  # number of slots
        p = poles/2   # number of pole pairs


        # let ask koil to compute the symmetrical and balanced winding
        #no single layer assumed

        w.compute_winding(m,Q,p,single_layer=False)

        phase_labels = ["U","V","W"]

        winding_pattern = [[i,i+1] for i in range(Q)]
        arr = []

        for i in range(len(w.windings)):

            _w = w.windings[i]
            _p = phase_labels[i]

            for c in _w.coils:

                if(c.nc==-1):

                    winding_pattern[c.begin-1][1] = _p+"-"
                    winding_pattern[c.end-1][0] = _p


                else:
                    winding_pattern[c.begin-1][1] = _p
                    winding_pattern[c.end-1][0] = _p+"-"



        #flatten winding pattern 2d array
        winding_pattern = [item for row in winding_pattern for item in row]

        return winding_pattern

    @staticmethod
    def getCircuitAndTurns(winding,no_turns):

        circuit = winding[0]

        if(len(winding)==1):
            turns = no_turns
        else:
            turns = -no_turns

        return [circuit,turns]
