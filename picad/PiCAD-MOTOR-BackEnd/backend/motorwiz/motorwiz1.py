import numpy as np
import pandas as pd

from .trig import *

from scipy.optimize import fsolve
import math


from .materials.load_materials import *

from .materials.thermal import *

from .topology.IPMSM.IPMSM_RADIAL import *

import json


#--------------------------- Read Me -------------------------------------------
#inputs to motor wiz
#1) Topology:
#   right now only supports IPMSM - Radial

#2) Peak Torque,Voltage Limit, Current Limit. Note current limit can be reduced

#3) mean_airgap_radius, TRV, airgap

#4) assumed winding temperature

#5) Number of slots, Number of Poles

#6) materials: (use config files unless specified as user input)

    #1) magnet -> supports N42
        #material type -> RE N42
        #i) b_r
        #ii)l_pm
        #iii)w_m
        #iv)k_pm
        #v)mu_ra
        #vi)density
        #vii) cost per g

    #2) wire -> 20 SWG copper

        #material type -> copper
        #1) wire dia ->user input
        #2) fill factor -> user input
        #3) tooth_depth factor -> user input
        #4) resistivity
        #5) temp. coeff of resistivity
        #6) density of conductor
        #7) cost_per_g

    #3) steel -> m350-50A

        #material type -> iron
        #1)B_sat
        #2)grade of steel
        #3)steel thickness
        #4)density
        #5)rs_per_g
        #6)steel loss

#7) Thermals:

    #1) HTC btw wire and steel -> in this case cu-fe
    #2) HTC with air -> based on cooling type -> user input
    #3) Fin area factor -> user input
    #4) housing/stack length -> user input
    #5) ambient temp -> user input

#output:

#1)motor geometry + cost
#2)Performance Analysis

#3) Next to motor design whose input will be the derived motor geomtry



#-------------------------------------------------------------------------------

#---------------------------Load Data-------------------------------------------

cooling_type_dict = {

    "natural_convection": 10,
    "forced_convection": 30,
}

#---------------------------------- class --------------------------------------
class MotorWiz:


    def __init__(self,peak_torque,peak_speed, current_limit_i_ph, voltage_limit, bridge_thickness, trv, ag, no_of_poles, no_of_slots, fill_factor, tooth_depth_factor, fin_area_factor, housing_by_stacklength_ratio,external_htc, ambient_temperature ,magnet,steel, wire, topology,motor_speed_rpm=0, motor_torque=0, assumed_slot_width_ratio=0.75, assumed_winding_temp = 158):


        #numerical values
        num_params = [peak_torque, current_limit_i_ph, voltage_limit, bridge_thickness, trv, ag, no_of_poles, no_of_slots, fill_factor, tooth_depth_factor, fin_area_factor, housing_by_stacklength_ratio, external_htc, ambient_temperature ,assumed_slot_width_ratio, assumed_winding_temp]

        i=0

        for p in num_params:

            if(type(p)!=int and type(p)!=float):
                raise ValueError("Entered None Numerical Values at param position {}".format(i))


        #magnet, steel, wire are tuples in the format:
        #("material_name", "material type")

        if(type(magnet)!=list):
            raise ValueError("Motor Wiz: Magnet not in list format")
        else:
            self.magnet_name = magnet[0]
            self.magnet_type = magnet[1]

        if(type(wire)!=list):
            raise ValueError("Motor Wiz: Wire not in list format")
        else:
            #make changes here
            #remove wire dia from params and calculate dia from wire_name
            #("20 SWG", "COPPER")
            self.wire_name = wire[0]
            self.wire_type = wire[1]


        if(type(steel)!=list):
            raise ValueError("Motor Wiz: Steel not in list format")
        else:
            #steel ("M350-50A", "IRON")
            self.steel_name = steel[0]
            self.steel_type = steel[1]


        if(type(topology)!=list):
            raise ValueError("Motor Wiz: topology not in list format")

        else:
            # topology = ("IPMSM","RADIAL")
            self.topology_name = topology[0]
            self.topology_type = topology[1]


        self.peak_torque = peak_torque
        self.peak_speed = peak_speed
        self.current_limit_i_ph = current_limit_i_ph
        self.voltage_limit = voltage_limit

        self.trv = trv
        self.ag = ag
        self.no_of_poles = no_of_poles
        self.no_of_slots = no_of_slots
        #self.pole_arc = pole_arc
        self.fill_factor = fill_factor

        self.tooth_depth_factor = tooth_depth_factor
        self.fin_area_factor = fin_area_factor
        self.housing_by_stacklength_ratio = housing_by_stacklength_ratio

        self.ambient_temperature = ambient_temperature
        self.assumed_slot_width_ratio = assumed_slot_width_ratio
        self.assumed_winding_temp = assumed_winding_temp
        self.external_htc = external_htc


        self.bridge_thickness = bridge_thickness

        #0.5 is rotordia/stacklength ratio
        self.r_ag = (((peak_torque*0.5)/(2*np.pi*trv))**(1/3))*10

        l_pm = 2*np.pi*self.r_ag/self.no_of_poles*0.80

        if(l_pm<15):
            l_pm=15



        self.magnet = Magnet(self.magnet_name)
        self.magnet.set_dimensions(l_pm,l_pm/5)

        self.steel = Steel(self.steel_name)

        self.wire = Wire(self.wire_type,self.fill_factor,self.tooth_depth_factor)
        self.wire_diameter = self.wire.wire_diameter

        self.thermal = Thermals(self.fin_area_factor,self.housing_by_stacklength_ratio,self.ambient_temperature,self.external_htc,self.wire_type,self.steel_type)



    def run(self,motor_speed_rpm=0,motor_torque=0,updated_bavg = False,verbose=True):

        #runs motor wiz
        #create motor geometry object in the end

        #re run below line after modifying current_limit and then store it
        #1)
        self.kt,self.ke,self.max_no_load_rpm = self.calculate_max_no_load_rpm(self.peak_torque,self.current_limit_i_ph,self.voltage_limit)


        #2) stack lenght will also be apart of motor geometry
        #self.stack_length = self.peak_torque/self.trv*1000/(np.pi*self.r_ag**2)

        #new calc for stack length 0.5 is rotor_dia/stack length

        self.stack_length = 2*self.r_ag/0.5

        self.σ = self.trv/2*(10**6)



        #3)make assumed slot width = derived slot width

        '''
        self.r = self.ag*self.magnet.w_m*self.no_of_poles/(2*np.pi*self.r_ag*self.magnet.l_pm)

        self.slot_pitch = 2*np.pi*self.r_ag/self.no_of_slots

        #to make w1 = w
        solution = fsolve(self.solve_slot_width, self.assumed_slot_width_ratio, args=(self.r,self.slot_pitch,self.ag,self.magnet.k_pm,self.magnet.b_r,self.magnet.mu_ra,self.σ))
        self.assumed_slot_width_ratio = solution[0]


        self.w,self.a,self.airgap_coeff,self.k_c,self.b_g,self.b_avg,self.K,self.slot_opening_width = self.calculate_slot_width_ratio(self.r,self.slot_pitch,self.assumed_slot_width_ratio,self.ag,self.magnet.k_pm,self.magnet.b_r,self.magnet.mu_ra,self.σ)
        '''

        self.w = self.assumed_slot_width_ratio

        if(not updated_bavg):
            #first time calling run -> assume b_avg
            self.b_avg = self.steel.b_sat*(1-self.assumed_slot_width_ratio)

        self.K = (4/np.pi)*(self.σ/self.b_avg)/1000

        self.slot_pitch = 2*np.pi*self.r_ag/self.no_of_slots

        self.i_total = round(self.K*2*np.pi*self.r_ag/1000)


        #make no. of strands whole number by changing i_ph

        solution = fsolve(self.solve_number_of_turns,self.current_limit_i_ph, args=(self.no_of_slots,self.i_total))
        self.current_limit_i_ph = solution[0]

        #update values after modifying a_rms as mention in 1)
        self.kt,self.ke,self.max_no_load_rpm = self.calculate_max_no_load_rpm(self.peak_torque,self.current_limit_i_ph,self.voltage_limit)

        self.number_of_turns = round((self.i_total*1000)/(2*self.no_of_slots*self.current_limit_i_ph))

        self.tooth_thickness = round(2*np.pi*self.r_ag*(1-self.w)/self.no_of_slots)

        self.yoke_thickness =round(1.2*(self.r_ag/self.no_of_poles))


        #4) Call IPMSM class and do all the necessary operation like, run, performance, cost and geometry

        # make changes here: this code base will only calculate for IPMSM_RADIAL
        # make it modular so other topologies can be calculated

        self.motor = IPMSM_RADIAL(self.magnet, self.steel, self.wire, self.thermal, self.tooth_depth_factor, self.r_ag,self.w,self.b_avg,self.current_limit_i_ph,self.σ,self.kt,self.ag,self.yoke_thickness,self.tooth_thickness,self.number_of_turns,self.no_of_slots,self.no_of_poles,self.slot_pitch,self.fill_factor,self.stack_length,self.assumed_winding_temp,self.ambient_temperature)

        #run motor
        self.motor.run(self.peak_speed,self.peak_torque)

        #update data in this case no strands and fill factor

        self.fill_factor = self.motor.fill_factor
        self.number_of_strands = self.motor.number_of_strands

        #display data

        if(verbose):
            self.display_all()

            self.magnet.display()
            self.steel.display()
            self.wire.display()

            self.motor.display_all(motor_speed_rpm,motor_torque)


    def calculate_max_no_load_rpm(self,peak_torque_τ,current_limit_i_ph,voltage_limit):

        #calculates maximum no load rpm
        #re run this when changing current_limit_i_ph

        kt = peak_torque_τ/current_limit_i_ph

        ke = kt/(3**0.5)*(1000*2*np.pi/60)

        max_no_load_rpm = voltage_limit/(2**0.5)/ke*1000

        return kt,ke,max_no_load_rpm


    def calculate_slot_width_ratio(self,r,slot_pitch,assumed_slot_width_ratio,ag,k_pm,b_r,μ_ra,σ):

        #calculates slot width ratio
        slot_opening_width = slot_pitch*assumed_slot_width_ratio

        a = slot_opening_width/(2*ag)

        airgap_coeff = 4/np.pi*(a*np.arctan(a)-np.log(np.sqrt(1+a**2)))

        k_c = slot_pitch/(slot_pitch-airgap_coeff*ag)

        b_g = k_pm*b_r/(1+μ_ra*k_c*r)

        b_avg = (2/np.pi)*b_g

        K = (4/np.pi)*(σ/b_avg)/1000

        w = 1-b_avg/self.steel.b_sat

        return [w,a,airgap_coeff,k_c,b_g,b_avg,K,slot_opening_width]


    def solve_slot_width(self,assumed_slot_width_ratio,r,slot_pitch,ag,k_pm,b_r,μ_ra,σ):

        #solves for slot width by changing assumed_slot_width_ratio
        w = self.calculate_slot_width_ratio(r,slot_pitch,assumed_slot_width_ratio,ag,k_pm,b_r,μ_ra,σ)[0]

        return w/assumed_slot_width_ratio-1

    def solve_number_of_turns(self,current_limit_i_ph,no_of_slots,i_total):

        #solves for  number of turns by adjusting phase current limit
        number_of_turns = (i_total*1000)/(2*no_of_slots*current_limit_i_ph)

        return np.rint(number_of_turns)-number_of_turns


    def display_inputs(self):

        #display input data
        print("\nMotorwiz Inputs")
        print("peak_torque:",self.peak_torque)
        print("current limit:", self.current_limit_i_ph)
        print("voltage_limit:", self.voltage_limit)
        print('mean airgap radius:',self.r_ag)
        print("TRV:",self.trv)

        print("air gap:", self.ag)
        print("no. poles:",self.no_of_poles)
        print("no. slots:", self.no_of_slots)
        print("assumed_slot_width_ratio", self.assumed_slot_width_ratio)

        print("fill factor:", self.fill_factor)
        print("wire dia:", self.wire_diameter)
        print("tooth_depth_factor:", self.tooth_depth_factor)

        print("External HTC",self.external_htc)
        print("fin_area_factor:", self.fin_area_factor)
        print("housing/stack length:", self.housing_by_stacklength_ratio)
        print("ambient_temperature:",self.ambient_temperature)
        print("assumed_winding_temp:", self.assumed_winding_temp)

        print("topology:",self.topology_name,self.topology_type)
        print("steel:",self.steel_name,self.steel_type)
        print("magnet:",self.magnet_name,self.magnet_type)
        print("wire:",self.wire_name,self.wire_type)


    def display_calculated_vars(self):

        #displays input variables
        print("\nCalculated Variables:")
        print('Kt (Nm/A_rms)',self.kt)
        print('Ke (V/Krpm)',self.ke)
        print('Max no-load rpm',self.max_no_load_rpm)
        print('Stack length (L)',self.stack_length)
        print('σ (N/m2)',self.σ)

        '''
        print('r = g*W_m*P/[2π*r_g*L_pm]',self.r)
        print('Slot pitch p =  2π*r_g/S',self.slot_pitch)
        print('Slot opg width W_s = p*w',self.slot_opening_width)
        print('a = W_s/(2g)',self.a)
        print('Airgap coeff γ',self.airgap_coeff)
        print('Carter Coeff K_c = p/(p-γg)',self.k_c)
        print('Peak Airgap flux B_g ',self.b_g)
        print('B_avg =(2/π)*B_g',self.b_avg)
        print('K (kA/m) = (4/π)*(σ/B_avg)',self.K)
        print('Actual slot-width ratio w',self.w)
        '''
        print('Actual slot-width ratio w',self.w)
        print('B_avg =(2/π)*B_g',self.b_avg)
        print('K (kA/m) = (4/π)*(σ/B_avg)',self.K)
        print('I_tot (KA) = K*(2πr_g)',self.i_total)
        print('Turns N = I_tot/(2S*I_ph)',self.number_of_turns)
        print('Tooth thickness t ',self.tooth_thickness)
        print('Yoke thickness y',self.yoke_thickness)



    def display_all(self):

        self.display_inputs()
        self.display_calculated_vars()
