import numpy as np
import pandas as pd

from .trig import *

from scipy.optimize import fsolve
import math


from .materials.load_materials1 import *

from .materials.thermal1 import *

from .topology.RADIAL.INNER_ROTOR import *

import json

import warnings
warnings.filterwarnings("ignore")


#--------------------------- Read Me -------------------------------------------
#inputs to motor wiz
#1) Topology:
#   right now only supports IPMSM - Radial and SPMSM

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
class MotorWiz_Radial_InnerRotor:


    def __init__(self, peak_torque,max_no_load_rpm, voltage_limit, trv, aspect_ratio, no_of_poles, no_of_slots, coil_span, g_airgap, L_pm, tooth_height, fill_factor,  ambient_temperature,htc_cond_housing_stator, htc_cond_slotpaper, htc_conv, housing_stack_length_diff, fin_area_factor, magnet, steel, wire, topology, assumed_slot_width_ratio=0.55, assumed_winding_temp = 60,overhang_looseness = 1.5,margin_harmonics= 1.3,t_tot_max = 100):


        #numerical values

        num_params = [peak_torque,max_no_load_rpm, voltage_limit, trv, aspect_ratio, no_of_poles, no_of_slots, coil_span, g_airgap, L_pm, tooth_height, fill_factor,  ambient_temperature, htc_cond_housing_stator, htc_cond_slotpaper, htc_conv, housing_stack_length_diff, fin_area_factor, assumed_slot_width_ratio, assumed_winding_temp,overhang_looseness,margin_harmonics,t_tot_max]
        i=0

        for p in num_params:

            if(type(p)!=int and type(p)!=float):
                raise ValueError("Entered None Numerical Values at param position {}".format(i))


        #magnet, steel, wire are tuples in the format:
        #("material_name", "material type")

        if(type(magnet)!=list):
            raise ValueError("Motor Wiz: Magnet not in list format")
        else:
            #("N42","RE")
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
            # topology = ("IPMSM","RADIAL", "INTERIOR"), ("IPMSM","RADIAL", "SURFACE")
            
            self.topology_name = topology[0]
            self.topology_type = topology[1]
            self.magnet_position = topology[2]
            
        #parameter
        self.τ_max = peak_torque
        self.max_no_load_rpm = max_no_load_rpm
        
        self.V_DC = voltage_limit
        
        #sizing
        self.TRV_max = trv
        self.aspect_ratio_k = aspect_ratio
        
        #finding B_av
        self.no_poles = no_of_poles
        self.no_slots = no_of_slots
        self.assumed_slot_width_ratio = assumed_slot_width_ratio
        self.coil_span = coil_span
        self.g_airgap = g_airgap
        self.L_pm = L_pm
        
        #R_ph and Km
        self.tooth_height = tooth_height
        self.fill_factor = fill_factor
        self.margin_harmonics = margin_harmonics
        self.overhang_looseness = overhang_looseness
        
        #thermal
        self.ambient_temperature = ambient_temperature
        self.assumed_winding_temp = assumed_winding_temp
        self.t_tot_max = t_tot_max 
        self.htc_cond_housing_stator = htc_cond_housing_stator
        self.htc_cond_slotpaper = htc_cond_slotpaper 
        
        self.htc_conv = htc_conv
        self.housing_stack_length_diff = housing_stack_length_diff 
        self.fin_area_factor = fin_area_factor
        

        #materials
        
        self.magnet = Magnet(self.magnet_name)
        #remember to update magnet dimensions after running calculate params fn

        self.steel = Steel(self.steel_name)

        self.wire = Wire(self.wire_name, self.wire_type,self.fill_factor)
        self.wire_diameter = self.wire.wire_diameter

        #thermal
        self.thermal = Thermals(self.fin_area_factor,self.ambient_temperature,self.htc_cond_housing_stator, self.htc_cond_slotpaper,self.htc_conv, self.t_tot_max, self.wire_type, self.steel_type)
    


    def run(self,b_g=0,updated_bavg = False,verbose=True):

        #runs motor wiz
        #create motor geometry object in the end

        #re run below line after modifying current_limit and then store it
        #1) calculate params
        self.calculate_params()
        

        #2) set magnet dimensions
        self.magnet.set_dimensions(self.L_pm, self.W_m)
        
        
        #3)make assumed slot width = derived slot width

        #to make w1 = w
        self.calculate_slot_width(self.assumed_slot_width_ratio)
        
        #with warnings.catch_warnings(action="ignore"):
        solution = fsolve(self.solve_slot_width, self.assumed_slot_width_ratio,xtol=0.1)
        
        self.assumed_slot_width_ratio = solution[0]

        self.w = self.assumed_slot_width_ratio
        
        self.calculate_slot_width(self.w)
        

        if(updated_bavg):
            #b_avg from femm
            
            self.B_g = b_g
            self.B_av =  2 / np.pi * self.B_g
            self.w = 1 - self.B_av / self.steel.b_sat


        #4)find number of turns and make it a whole number
        self.calculate_turns()
        
        turns = self.Turns
        
        self.Turns = round(self.Turns)
        
        if(self.Turns==0):
            self.Turns = 1
       
        I_max_required = (np.pi*self.r_g*self.k_rms)/(int(self.Turns)*self.no_slots)
        #print("-------------Start-------------")
        #print("before goalseek \nturns {}, turns_round {},\nmaxnl {},\nimax {}, imax_req {}".format(turns, self.Turns, self.max_no_load_rpm, self.I_max, I_max_required))

        #with warnings.catch_warnings(action="ignore"):
        solution = fsolve(self.solve_Imax, self.max_no_load_rpm, args=(I_max_required),xtol=0.01)
            
    
        self.max_no_load_rpm = solution[0]
        
        
        turns = self.Turns
     
        self.Turns = int(self.Turns)
        
        
        #print("---------------------------")
        #update params
        self.calculate_params()
        
        #print("after goalseek \nturns {}, turns_round {},\nmaxnl {},\nimax {}, imax_req {}".format(turns, self.Turns, self.max_no_load_rpm, self.I_max, I_max_required))

        #print(self.τ_max)
        
        #5)find number of strands and make it a whole number, also calculate Km and Rph
        
        self.L_ph_actv = 2 * self.NS / 3 * self.L
        self.mean_slot_radius = self.r_g * (1 + self.tooth_height / 2)
        self.overhang = 2 * np.pi * self.mean_slot_radius * self.coil_span * 2 * self.Turns / 3 * self.overhang_looseness
        self.L_ph = self.L_ph_actv + self.overhang
        self.tooth_depth_d = self.tooth_height * self.r_g
        self.tooth_thickness = 2 * np.pi * self.r_g * (1 - self.w) / self.no_slots * self.margin_harmonics
        self.actual_slot_width_ratio = 1 - self.tooth_thickness / self.sp
        
        
        self.fill_factor = self.wire.fill_factor
        self.calculate_strands(self.wire.fill_factor)
        
        self.strands = round(self.strands)
        
        if(self.strands==0):
            self.strands = 1
        
        #with warnings.catch_warnings(action="ignore"):
        solution = fsolve(self.solve_strands, self.fill_factor, args=(self.strands),xtol=0.01)
        self.fill_factor = solution[0]
        self.wire.fill_factor = solution[0]
            
        self.calculate_strands(self.wire.fill_factor)
        self.strands = int(self.strands)
        

        #6) run motor (geometry and thermal)
        self.motor = INNER_ROTOR(self.magnet_position, self.magnet, self.steel, self.wire, self.thermal, self.r_g, self.tooth_thickness, self.tooth_height, self.tooth_depth_d, self.mean_slot_radius, self.stack_length, self.actual_slot_width_ratio, self.no_poles, self.no_slots, self.ambient_temperature, self.L_ph, self.A_ph, self.R_ph, self.Kt, self.housing_stack_length_diff)
        
        self.motor.create()      

        if(verbose):
            self.display_all()

            self.magnet.display()
            self.steel.display()
            self.wire.display()

            self.motor.display_all()


    def calculate_params(self):

        self.Ke = (self.V_DC / np.sqrt(2)) / (self.max_no_load_rpm / 1000)
        self.Kt = np.sqrt(3) * self.Ke * (3 / (100 * np.pi))
        self.I_max = self.τ_max / self.Kt
        self.σ = (self.TRV_max / 2 )* 10**6
        self.RV = self.τ_max / self.TRV_max * 10**3
        self.r_g = (self.RV *self.aspect_ratio_k / np.pi)**(1/3)
        self.L = self.r_g / self.aspect_ratio_k
        
        self.stack_length = self.L
        
        #width of magnet
        self.W_m = 0.8*2 * (self.r_g - self.g_airgap / 2) * np.sin (np.pi / self.no_poles)
        

    
    def calculate_slot_width(self,assumed_slot_width_ratio):

        self.rel = (self.g_airgap * self.W_m * self.no_poles) / (2 * np.pi * self.r_g * self.L_pm)
        self.sp = 2 * np.pi * self.r_g / self.no_slots
        self.W_s = self.sp * assumed_slot_width_ratio
        self.a = self.W_s / (2 * self.g_airgap)
        self.γ = 4 / np.pi * (self.a * np.arctan(self.a) - np.log (np.sqrt(1 + self.a**2)))
        
        self.K_c = self.sp / (self.sp - self.γ * self.g_airgap)
        
        
        #new K_c calculation
        #print(1+(np.pi*self.W_s)/(4*self.g_airgap))
        #self.K_c = 1/(1-self.W_s/self.sp+(4*self.g_airgap)/(np.pi*self.sp)*np.log(1+(np.pi*self.W_s)/(4*self.g_airgap)))
        
        self.B_g = self.magnet.k_pm * self.magnet.b_r / (1 + self.magnet.mu_ra* self.K_c * self.rel)
        self.B_av =  2 / np.pi * self.B_g
        self.w = 1 - self.B_av / self.steel.b_sat
        self.slot_width = self.w
        
        #print("calc slot width", assumed_slot_width_ratio, self.rel, self.sp,self.W_s,self.a,self.γ,self.ert2,self.B_g, self.B_av,self.w,self.slot_width)
        

    def solve_slot_width(self,assumed_slot_width_ratio):
        
        self.calculate_slot_width(assumed_slot_width_ratio)
        
        slot_width = self.slot_width

        return slot_width-assumed_slot_width_ratio


    def calculate_turns(self):
        
        self.B_av = (2 / np.pi) * self.B_g
        self.k_rms = 2 * np.sqrt(2) / np.pi * self.σ / self.B_av/ 1000
        self.NS = np.pi * self.r_g * self.k_rms / self.I_max
        self.Turns = self.NS /self.no_slots
        
    def solve_Imax(self,max_no_load_rpm,I_max_required):
       
        self.max_no_load_rpm = max_no_load_rpm
        self.calculate_params()
        i_max = self.I_max
        return i_max-I_max_required
    

    def calculate_strands(self,fill_factor):
    
        
        self.A_ph = (2 * np.pi * self.mean_slot_radius - self.no_slots * self.tooth_thickness) * self.tooth_depth_d * fill_factor / (2 * self.Turns * self.no_slots)
        
        self.R_ph = self.wire.resistivity * self.L_ph / self.A_ph * 10 ** 6
        self.j_max = self.I_max / self.A_ph
        self.Km = self.Kt / np.sqrt(3 * self.R_ph * 10 ** (-3))
        self.strands = self.A_ph / (np.pi * self.wire.wire_diameter ** 2 / 4)
    

    def solve_strands(self,fill_factor,strands_required):
    
        self.calculate_strands(fill_factor)
        strands = self.strands
        return strands-strands_required

    
        

    def display_inputs(self):

        #display input data
        print("\nMotorwiz Inputs")
        print('peak torque',self.τ_max)
        print('max no load rpm',self.max_no_load_rpm)
        print('voltage limit',self.V_DC)
        print('trv max',self.TRV_max)
        print('aspect ratio',self.aspect_ratio_k)
        print('no of poles',self.no_poles)
        print('no of slots',self.no_slots)
        print('coil span',self.coil_span)
        print('g airgap',self.g_airgap)
        print('L pm',self.L_pm)
        print('tooth height factor',self.tooth_height)
        print('fill factor',self.wire.fill_factor)
        print('ambient temperature',self.ambient_temperature)
        print('htc cond slotpaper',self.htc_cond_slotpaper)
        print('htc conv',self.htc_conv)
        print('assumed slot width ratio',self.assumed_slot_width_ratio)
        print('assumed winding temp',self.assumed_winding_temp)
        print('overhang looseness',self.overhang_looseness)
        print('margin harmonics',self.margin_harmonics)
        print('t tot max',self.t_tot_max)

        print("topology:",self.topology_name,self.topology_type, self.magnet_position)
        print("steel:",self.steel_name,self.steel_type)
        print("magnet:",self.magnet_name,self.magnet_type)
        print("wire:",self.wire_name,self.wire_type)


    def display_calculated_vars(self):

        #displays input variables
        print("\nCalculated Variables:")
        print('Kt (Nm/A_rms)',self.Kt)
        print('Ke (V/Krpm)',self.Ke)
        print('Max no-load rpm',self.max_no_load_rpm)
        print('Stack length (L)',self.stack_length)
        print('σ (N/m2)',self.σ)
        print('I max (A RMS)',self.I_max)
        print('Rotor Volume RV (mm3)',self.RV)
        print('R_g (mm)',self.r_g)
        
    

        
        print('r = g*W_m*P/[2π*r_g*L_pm]',self.rel)
        print('Slot pitch p =  2π*r_g/S',self.sp)
        print('Slot opg width W_s = p*w',self.W_s)
        print('a = W_s/(2g)',self.a)
        print('Airgap coeff γ',self.γ)
        print('Carter Coeff K_c = p/(p-γg)',self.K_c)
        print('Peak Airgap flux B_g ',self.B_g)
        print('B_avg =(2/π)*B_g',self.B_av)
        print('Actual slot-width ratio w',self.w)
        
        
        print('K_RMS ', self.k_rms)
        print('NS',self.NS)
        print('Turns ',self.Turns)
    
        print("L_ph_act", self.L_ph_actv)
        print("Mean Slot Radius",self.mean_slot_radius)
        print("Overhang", self.overhang)
        print("tooth depth",self.tooth_depth_d)
        print("tooth thickness", self.tooth_thickness)
        print("L_ph", self.L_ph)
        print('A PH',self.A_ph)
        print('R PH',self.R_ph)
        print('J MAX',self.j_max)
        print('KM',self.Km)
        print("Strands",self.strands)
        




    def display_all(self):

        self.display_inputs()
        self.display_calculated_vars()


