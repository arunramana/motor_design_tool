import ezdxf
from random import random
from ezdxf.addons import r12writer
import numpy as np
from numpy import pi

from trig import *
from dolomites import koil


#SPMSMDXF is a class that draws the dxf of a motor (SPMSM for now)
#uses ezdxf library
#The motor dimensions can be drawn based on given input variables
#You can draw a custom motor based on user input or
#Getting dimensions from motor wiz
#The generated dxf can then be imported into FEMM

class SPMSMDXF:

    def __init__(self, filename, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns):

        #file
        self.directory = filename

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

        self.no_turns = no_turns


        #-----------------------------calculated variables ---------------------
        self.has_notch = notch_angle!=0 or notch_depth!=0


        self.stator_or = stator_od/2
        self.stator_ir = stator_id/2
        
        self.rotor_or = rotor_od/2
        self.rotor_or_magnet = rotor_od/2
        self.rotor_or_lamination = self.rotor_or_magnet-2*pole_width
        
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



    def draw_motor(self):
        #draw the full motor
        if(self.stator_id <= self.rotor_od):
            print("Error stator ID <= rotor_od")



        self.dxf_file = self.directory+"/motor.dxf"

        with r12writer(self.dxf_file) as doc:

            #stator od
            doc.add_circle((self.xo,self.yo), radius=self.stator_or)


            #stator id
            #doc.add_circle((xo,yo), radius=stator_id)

            #rotor od

            if(not self.has_notch):
                doc.add_circle((self.xo,self.yo), radius=self.rotor_or_lamination)

            #shaft
            doc.add_circle((self.xo,self.yo), radius=self.rotor_or_lamination)
            doc.add_circle((self.xo,self.yo), radius=self.shaft_rad)

            # doc.add_circle((xo,yo), radius=stator_od-yoke_thickness,color = 1)

            #draw stator
            self.draw_stator(doc,True)

            #draw rotor
            self.draw_rotor(doc)


        #winding diagram

        winding_filename = self.directory+"/motor_winding.dxf"

        with r12writer(winding_filename) as doc:

            #stator od
            doc.add_circle((self.xo,self.yo), radius=self.stator_or)


            self.draw_stator(doc,True)

            self.draw_winding(doc,self.no_slots,self.no_poles,self.theta_inc_stator,self.stator_ir,self.stator_or,self.yoke_thickness,self.tooth_width)



    #---------------------------FUNCTIONS-------------------------------------------
    #------------------------- STATOR ----------------------------------------------
    def draw_stator(self,doc,enable_femm=False):
        #draw stator
        reference_angle = 0 # in degrees

        tooth_coord = []

        for no in range(0,self.no_slots):

            reference_angle += self.theta_inc_stator

            start_angle = reference_angle-self.tooth_arc_theta_deg
            end_angle = reference_angle+self.tooth_arc_theta_deg

            #draw tooth arc
            self.draw_tooth_arc(doc,self.xo,self.yo,self.stator_ir,start_angle,end_angle)

            if(enable_femm):
                self.connect_slot_opening(doc,self.xo,self.yo,self.no_slots,self.stator_ir,self.stator_or,self.theta_inc_stator,self.yoke_thickness,self.slot_opening)
                #doc.add_circle((xo,yo), radius=stator_ir)

            #draw tooth tip
            coord = self.draw_tooth_tip(doc,self.stator_ir,start_angle,end_angle,self.tooth_tip_angle,self.tooth_width,self.tooth_arc_s)

            #draw tooth base
            coord = self.draw_tooth_base(doc,self.stator_or,reference_angle,self.yoke_thickness,self.tooth_width,coord)

            #store the coordinates to connect the tooth
            tooth_coord.append(coord)


        #connect all the tooth

        r0,theta0,r1,theta1 = tooth_coord[0]


        for i in range(1,len(tooth_coord)):

            r2,theta2,r3,theta3 = tooth_coord[i]

            doc.add_arc((self.xo,self.yo),radius=r0, start=theta1+90, end=theta2+90)

            r1,theta1=r3,theta3



        doc.add_arc((self.xo,self.yo),radius=r0, start=theta1+90, end=theta0+90)


    def draw_tooth_arc(self,doc,xo,yo,radius,start_angle,end_angle):
        #draw arc of the stator tooth that faces the the rotor
        #radius is stator ir
        doc.add_arc((xo,yo),radius=radius, start=start_angle+90, end=end_angle+90)

    def draw_tooth_tip(self,doc,radius,start_angle,end_angle,tooth_tip_angle,tooth_width,tooth_arc_s):

        #draws  2 triangular tips of the stator tooth

        w = (tooth_arc_s -tooth_width)/2 # approximation

        r2 = w*tan(tooth_tip_angle)+radius

        theta2 = start_angle+(180*w)/(pi*radius)
        theta3 = end_angle-(180*w)/(pi*radius)


        x1 = radius*sin(start_angle)
        y1 = radius*cos(start_angle)

        x2 = r2*sin(theta2)
        y2 = r2*cos(theta2)


        x3 = radius*sin(end_angle)
        y3 = radius*cos(end_angle)

        x4 = r2*sin(theta3)
        y4 = r2*cos(theta3)

        doc.add_line((x1,y1),(x2,y2))
        doc.add_line((x3,y3),(x4,y4))

        #pass the below array as coordinates for draw_tooth_base fn
        return [x2,y2,x4,y4]

    def connect_slot_opening(self,doc,xo,yo,no_slots,stator_ir,stator_or,theta_inc_stator,yoke_thickness,slot_opening):

        #Connects stator tooth opening and draws a central line (half slots) for femm.

        reference_angle = 0

        str_or = stator_or-yoke_thickness

        theta = (slot_opening/stator_ir)*(180/pi)


        for no in range(0,no_slots):

            reference_angle += theta_inc_stator/2

            x1 = stator_ir*sin(reference_angle)
            y1 = stator_ir*cos(reference_angle)

            x2 = str_or*sin(reference_angle)
            y2 = str_or*cos(reference_angle)

            start_t = reference_angle + theta/2+90
            end_t = reference_angle - theta/2+90

            doc.add_line((x1,y1),(x2,y2), color=1)

            doc.add_arc((xo,yo),radius=stator_ir, start=end_t, end=start_t,color=1)

            reference_angle += theta_inc_stator/2

    def draw_tooth_base(self,doc,stator_or,reference_angle,yoke_thickness,tooth_width,coord):

        #draws the base of the stator tooth

        r = stator_or-yoke_thickness

        x1,y1,x2,y2 = coord

        w = tooth_width/2

        theta = (w/r)*(180/pi)

        theta1 = reference_angle-theta
        theta2 = reference_angle+theta

        x3 = r*sin(theta1)
        y3 = r*cos(theta1)


        x4 = r*sin(theta2)
        y4 = r*cos(theta2)

        doc.add_line((x1,y1),(x3,y3))
        doc.add_line((x2,y2),(x4,y4))

        #pass the below array to connect_tooth fn
        return r,theta1,r,theta2


    #------------------------- ROTOR------------------------------------------------
    def draw_rotor(self,doc):

        #draws the rotor

        reference_angle = 0 # in degrees

        notch_coord = []

        for iter_no in range(0,self.no_poles):

            reference_angle += self.theta_inc_rotor

            r2 = self.draw_pole(doc,self.rotor_or_magnet,self.rotor_or_lamination,self.pole_length,self.pole_width,self.pole_dist,self.pole_angle,reference_angle)



    def draw_pole(self,doc,rotor_or_magnet,rotor_or_lamination,pole_length,pole_width,pole_distance,pole_angle,reference_angle):
        

        #draws the pole(magnet) of the rotor
        theta = (pole_length/(2*rotor_or_lamination))*180/np.pi
        
        r1 = rotor_or_magnet
        r2 = rotor_or_lamination
        
        theta1 = reference_angle-theta
        theta2 = reference_angle+theta
        
        
        x0,y0=0,0
        
        #outer magnet points/arc
        x1 = r1*sin(theta1)
        y1 = r1*cos(theta1)
        
        x2 = r1*sin(theta2)
        y2 = r1*cos(theta2)
        
        #inner magnet points/arc
        x3 = r2*sin(theta1)
        y3 = r2*cos(theta1)
        
        x4 = r2*sin(theta2)
        y4 = r2*cos(theta2)
        
        
        #two half arcs
        doc.add_arc((x0,y0),radius=r1, start=theta1+90, end=reference_angle+90)
        doc.add_arc((x0,y0),radius=r1, start=reference_angle+90, end=theta2+90) 
        
        doc.add_line((x2,y2),(x4,y4))
        
        doc.add_line((x3,y3),(x1,y1))
        
        #return r2 for hole fn
        return r2

    def draw_rotor_ducts(self,doc,reference_angle,duct_dist,duct_angle,duct_radius,pole_dist,r2,iter_no):

        #draws ducts in the rotor (optional)
        #duct_angle +ve if right of reference_angle nad -ve if left
        #r2 comes from draw_pole fn

        #skip condition

        if(duct_dist==0 or duct_radius==0):
            #no duct
            return


        if(iter_no%2!=0):
            return

        if(r2==None):
            r2 = pole_dist

        if(pole_dist <= duct_dist or r2 <= duct_dist):
            print("Error: Change Duct Distance to prevent overlapping")

        x = duct_dist*sin(reference_angle+duct_angle)
        y = duct_dist*cos(reference_angle+duct_angle)


        doc.add_circle((x,y), radius=duct_radius)

    def draw_rotor_notch(self,doc,reference_angle,no_poles,has_notch,notch_depth,notch_angle,rotor_or):

        #draws notches on the rotor (optional)
        
        if(not has_notch):
            return
        
        if(notch_depth>=rotor_or):
            print("Invalid value for notch depth")
            return
        
        r1 =rotor_or-notch_depth
        theta1 = reference_angle-(360/no_poles)/2 #plus because only the left notch is drawn
        
        r2 = rotor_or
        theta2 = (theta1+notch_angle/2)
        theta3 = (theta1-notch_angle/2)
        
        x1 = r1*sin(theta1)
        y1 = r1*cos(theta1)
        
        x2 = r2*sin(theta2)
        y2 = r2*cos(theta2)
        
        x3 = r2*sin(theta3)
        y3 = r2*cos(theta3)
        
        
        doc.add_line((x1,y1),(x2,y2))#right
        doc.add_line((x1,y1),(x3,y3))#left
        
        
        #return coords to draw arcs
        return [r2,theta3,r2,theta2]


    def draw_winding(self,doc,no_slots,no_poles,theta_inc_stator,stator_ir,stator_or,yoke_thickness,tooth_width):

        #Generates Winding Diagram for reference

        reference_angle = 0

        i=0


        winding_pattern = SPMSMDXF.getWindingPattern(no_slots,no_poles)

        coords = []

        for no in range(no_slots):

            theta = (reference_angle-theta_inc_stator)/2

            r = (stator_ir+(stator_or-yoke_thickness))/2

            #approx
            theta1 = ((tooth_width/r)*(180/pi))/2

            x2 = r*sin(theta+theta1)
            y2 = r*cos(theta+theta1)


            x1 = r*sin(theta-theta1)
            y1 = r*cos(theta-theta1)

            #--------------------------------------------------------


            winding = winding_pattern[i]
            incircuit = SPMSMDXF.getCircuitLabels(winding)


            doc.add_point((x1,y1))
            doc.add_text(incircuit,(x1, y1-2),height=2,width=2)

            i+=1

            #--------------------------------------------------------

            winding = winding_pattern[i]
            incircuit  = SPMSMDXF.getCircuitLabels(winding)

            doc.add_point((x2,y2))
            doc.add_text(incircuit,(x2, y2-2),height=2,width=2)

            i+=1

            reference_angle += self.theta_inc_stator*2

            coords.append([x1,y1])
            coords.append([x2,y2])



        i=1

        while i < len(coords):


            if(i==len(coords)-1):
                x1,y1 = coords[i]
                x2,y2 = coords[0]

            else:
                x1,y1 = coords[i]
                x2,y2 = coords[i+1]


            #curve fitting


            if("U" in winding_pattern[i]):
                color=1
            elif("V" in winding_pattern[i]):
                color=4
            else:
                color=3


            doc.add_line((x1,y1),(x2,y2), color=color)

            i+=2



    @staticmethod
    def getWindingPattern(slots,poles):

        #generates winding pattern using Koil Library

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
    def getCircuitLabels(winding):

        #puts labels based on the winding pattern

        circuit = winding[0]

        if(len(winding)==1):
            circuit=circuit+"+"
        else:
            circuit=circuit+"-"

        return circuit


#Default values
#-------------------------------Input Variables --------------------------------------
no_slots = 12

no_poles = 10

xo,yo = 0,0

stator_od = 125
stator_id = 84
rotor_od = 83
shaft_dia = 20

slot_opening = 3

tooth_tip_angle = 20

tooth_width = 10

yoke_thickness = 5


pole_length = 18
pole_width = 2
pole_arc = 150

bridge_thickness = 1


#the three arrays below should be of same dimensions
duct_distances = [30,25,30]
duct_radii = [2,2,2]
duct_angles = [-18,0,18]
duct_distances = [0,0,0]

notch_depth = 0
notch_depth=0
notch_angle = 7


no_turns = 9

#-------------------------- Calculated Varaibles -------------------------------
has_notch = notch_angle!=0 or notch_depth!=0
has_notch = False
enable_femm = False

stator_or = stator_od/2
stator_ir = stator_id/2
rotor_or_magnet = rotor_od/2

rotor_or_lamination = rotor_or_magnet-2*pole_width
shaft_rad = shaft_dia/2

air_gap = stator_id - rotor_od

theta_inc_stator = 360/no_slots

stator_outter_circum = pi*stator_od
stator_inner_circum = pi*stator_id

stator_opening_total = stator_inner_circum - slot_opening*no_slots
#s = r*theta
tooth_arc_s = stator_opening_total/no_slots
tooth_arc_theta = (tooth_arc_s*2)/stator_id

tooth_arc_theta_deg = ((tooth_arc_theta*180)/pi)/2 #divide by 2 for start and end angle

rotor_outter_circum = pi*rotor_od

max_pole_length = rotor_outter_circum/no_poles

#pole outter points distance from origin
pole_dist = (rotor_od-2*bridge_thickness)/2
pole_angle = (pole_length/pole_dist)*(180/pi)


#multiply by 2 since pole arc is 180 deg
pole_angle = (pole_arc*2)/no_poles
pole_gap = ((180-pole_arc)*2)/no_poles



theta_inc_rotor = 360/no_poles

project_name = "data"

gdxf = SPMSMDXF(project_name, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width, bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns)

gdxf.draw_motor()


