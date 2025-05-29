import numpy as np

from scipy.optimize import fsolve
import math

from backend.motorwiz.materials.load_materials import *

import json


class IPMSM_RADIAL:

    def __init__(self, magnet, steel, wire, thermal, tooth_depth_factor, r_ag,w,b_avg,current_limit_i_ph,σ,kt,ag,yoke_thickness,tooth_thickness,number_of_turns,no_of_slots,no_of_poles,slot_pitch,fill_factor,stack_length,assumed_winding_temp,ambient_temperature):

        #materials
        self.magnet = magnet
        self.steel = steel
        self.wire = wire

        self.thermal = thermal

        #variables
        self.tooth_depth_factor = tooth_depth_factor
        self.r_ag = r_ag

        self.w = w
        self.b_avg = b_avg
        self.current_limit_i_ph = current_limit_i_ph
        self.σ = σ
        self.kt = kt
        self.wire_diameter = self.wire.wire_diameter
        self.fill_factor = fill_factor
        self.yoke_thickness = yoke_thickness
        self.tooth_thickness = tooth_thickness
        self.number_of_turns = number_of_turns

        self.no_of_slots = no_of_slots
        self.no_of_poles = no_of_poles

        self.slot_pitch = slot_pitch
        self.stack_length = stack_length

        self.assumed_winding_temp = assumed_winding_temp
        self.ambient_temperature = ambient_temperature


        #1) calculations
        self.tooth_depth = self.tooth_depth_factor*self.r_ag

        #first scaler z for radial
        self.z = (self.tooth_depth+2*self.r_ag)/(2*self.r_ag)
        self.ag = ag



        #make number of strands whole number by changing fillfactor -> change this in motorwiz as well fill factor and number of strands


        solution = fsolve(self.solve_number_of_strands,self.fill_factor,args=(self.w,self.tooth_depth,self.b_avg,self.current_limit_i_ph,self.z,self.σ,self.wire_diameter))
        self.fill_factor = solution[0]

        self.number_of_strands = round((self.fill_factor*self.w*self.tooth_depth*self.b_avg*self.current_limit_i_ph*self.z/(self.σ*(self.wire_diameter**2))*1000)/100)



        #2) calculate IPMSM geometry

        self.slot_depth = self.tooth_depth

        self.rotor_od = 2*self.r_ag-self.ag

        self.slot_od = self.rotor_od+2*self.ag+2*self.tooth_depth

        self.stator_od = self.slot_od+self.yoke_thickness*2

        self.rotor_id = self.rotor_od-2*self.yoke_thickness

        self.mean_slot_dia = (self.rotor_od+2*self.ag+self.slot_od)/2


        #3)calculate scaler2 z1:

        #calc for radial

        self.z2 = self.mean_slot_dia/(2*self.r_ag)

        self.l_ph = 2*self.number_of_turns*self.no_of_slots/3*(self.stack_length+self.slot_pitch)*self.z2

        self.a_ph = 2*np.pi*(self.fill_factor/100)*self.w*self.tooth_depth*self.r_ag/(2*self.number_of_turns*self.no_of_slots)*self.z

        self.r_ph = self.wire.resistivity*(1+self.wire.temp_coeff_of_resistivity/100*(self.assumed_winding_temp-25))*self.l_ph/self.a_ph*(10**6)*1.33

        self.performance_dict = {}


    def run(self,motor_speed_rpm,motor_torque):

        self.calculate_thermal_resistance()

        #calculate this before performance analysis
        self.calculate_cost_of_materials()

        torq_div = 5
        speed_div = 10

        torq_inc = motor_torque/torq_div
        speed_inc = motor_speed_rpm/speed_div

        torq_arr = []
        speed_arr = []
        eff_arr = []
        w_temp_arr = []
        e_freq_arr = []

        for t in range(1,torq_div+1):

            torq = t*torq_inc
            torq_arr.append(torq)

        for s in range(1,speed_div+1):

            speed = s*speed_inc
            speed_arr.append(speed)


        for torq in torq_arr:

            for speed in speed_arr:

                solution = fsolve(self.solve_winding_temperature,self.assumed_winding_temp,args=(speed,torq))
                self.assumed_winding_temp = solution[0]

                if(self.assumed_winding_temp!=-1):
                    self.performance_analysis(self.assumed_winding_temp,speed,torq)
                else:
                    self.assumed_winding_temp = 158 #default
                    self.winding_temperature = -1

                if(type(self.efficiency)==np.ndarray or type(self.efficiency)==list):
                    self.efficiency=self.efficiency[0]

                if(type(self.winding_temperature)==np.ndarray or type(self.winding_temperature)==list):
                    self.winding_temperature=self.winding_temperature[0]


                if(type(self.electrical_frequency)==np.ndarray or type(self.electrical_frequency)==list):
                    self.electrical_frequency=self.electrical_frequency[0]


                eff_arr.append(self.efficiency)
                w_temp_arr.append(self.winding_temperature)
                e_freq_arr.append(self.electrical_frequency)





        self.performance_dict = {"torque":torq_arr,"speed": speed_arr, "efficiency": eff_arr, "winding_temperature": w_temp_arr, "electrical_frequency": e_freq_arr}



    def update_rph(self,assumed_winding_temp):

        if(assumed_winding_temp==-1 or self.assumed_winding_temp==-1):
            assumed_winding_temp = 158 #default

        self.assumed_winding_temp = assumed_winding_temp
        self.r_ph = self.wire.resistivity*(1+self.wire.temp_coeff_of_resistivity/100*(assumed_winding_temp-25))*self.l_ph/self.a_ph*(10**6)*1.33



    def solve_number_of_strands(self,fill_factor,w,tooth_depth,b_avg,current_limit_i_ph,z,σ,wire_dia):

        number_of_strands = (fill_factor*w*tooth_depth*b_avg*current_limit_i_ph*z/(σ*(wire_dia**2))*1000)/100

        return np.rint(number_of_strands)-number_of_strands



    def calculate_cost_of_materials(self):

        #cost calculation radial

        self.stator_steel_mass = (np.pi*(self.stator_od**2-self.slot_od**2)/4+self.tooth_thickness*self.slot_depth*self.no_of_slots)*self.stack_length/1000*self.steel.density

        self.stator_steel_cost = self.stator_steel_mass*self.steel.rs_per_g


        self.magnet_mass = self.magnet.w_m*self.magnet.l_pm*self.no_of_poles*self.stack_length/1000*self.magnet.density

        self.magnet_cost = self.magnet_mass*self.magnet.rs_per_g


        self.rotor_steel_mass = (np.pi*(self.rotor_od**2-self.rotor_id**2)/4*self.stack_length-self.magnet_mass/self.magnet.density)*self.steel.density/1000

        self.rotor_steel_cost = self.rotor_steel_mass*self.steel.rs_per_g


        self.conductor_mass = 3*self.l_ph*self.a_ph/1000*self.wire.density

        self.conductor_cost = self.conductor_mass*self.wire.rs_per_g


        self.total_mass = self.stator_steel_mass+self.rotor_steel_mass+self.magnet_mass+self.conductor_mass
        self.total_cost = self.stator_steel_cost+self.rotor_steel_cost+self.magnet_cost+self.conductor_cost



    def calculate_thermal_resistance(self):

        #calc radial thermal resistance (tr)

        #celsius/W
        self.slot_perimeter = self.slot_pitch-self.tooth_thickness+2*self.slot_depth

        self.conduction_area = self.slot_perimeter*self.stack_length*self.no_of_slots*(10**-6)

        self.R_conduction = 1/self.thermal.wire_steel_htc/self.conduction_area

        self.convection_area = np.pi*self.stator_od*(self.stack_length*self.thermal.fin_area_factor*self.thermal.housing_by_stacklength_ratio+self.stator_od/2)*(10**-6)

        self.R_convection = 1/(self.convection_area*self.thermal.htc)





    def performance_analysis(self,assumed_winding_temp,motor_speed_rpm,motor_torque):
        #peformance_analysis
        #Losses

        #need to solve for r_ph
        #calc
        self.electrical_frequency = motor_speed_rpm/60*self.no_of_poles/2

        self.a_rms = motor_torque/self.kt

        self.cogging_frequency = self.electrical_frequency*self.no_of_slots

        self.calculate_thermal_resistance()

        #update r_ph using winding temp
        self.update_rph(assumed_winding_temp)


        self.copper_loss = 3*self.a_rms**2*self.r_ph/1000

        #self.steel_loss = self.stee #given for m350-50a

        #calculate steel loss based on rpm and stator weight

        self.steel.calculate_steel_loss(motor_speed_rpm,self.stator_steel_mass,self.no_of_poles)


        self.total_loss = self.copper_loss+self.steel.steel_loss

        #power kW

        self.output_power = motor_speed_rpm*motor_torque*2*np.pi/60/1000

        self.input_power = self.output_power+self.total_loss/1000

        self.efficiency = self.output_power/self.input_power


        #thermal steady state

        self.heat_flux = self.total_loss/self.convection_area*10**-4
        self.body_temperature = self.total_loss*self.R_convection+self.ambient_temperature

        if(self.heat_flux<=0.3):

            self.winding_temperature = self.copper_loss*self.R_conduction+self.body_temperature

        else:
            #unstable thermal condition
            self.winding_temperature = -1
            self.assumed_winding_temp = -1





    def solve_winding_temperature(self,assumed_winding_temp,motor_speed_rpm,motor_torque):

        self.assumed_winding_temp = assumed_winding_temp
        self.update_rph(self.assumed_winding_temp)
        self.performance_analysis(self.assumed_winding_temp,motor_speed_rpm,motor_torque)

        return self.assumed_winding_temp-self.winding_temperature



    def display_geometry(self):

        print("\nIPMSM RADIAL GEOMETRY")
        print("rotor_od:",self.rotor_od)
        print("tooth_thickness:",self.tooth_thickness)
        print("slot_od:",self.slot_od)
        print("slot_depth:",self.slot_depth)
        print("yoke_thickness:",self.yoke_thickness)
        print("stator_od:",self.stator_od)
        print("rotor_id:",self.rotor_id)
        print("mean_slot_dia:",self.mean_slot_dia)
        print("stack_length:",self.stack_length)


    def display_cost(self):

        self.calculate_cost_of_materials()

        print("\nCOST")
        print("Stator Steel: {} g, {} INR".format(self.stator_steel_mass,self.stator_steel_cost))
        print("Rotor Steel:  {} g, {} INR".format(self.rotor_steel_mass, self.rotor_steel_cost))
        print("Conductor:    {} g, {} INR".format(self.conductor_mass, self.conductor_cost))
        print("Magnet:       {} g, {} INR".format(self.magnet_mass, self.magnet_cost))
        print("TOTAL:        {} g, {} INR".format(self.total_mass, self.total_cost))



    def display_performance(self,motor_speed_rpm,motor_torque):

        self.performance_analysis(self.assumed_winding_temp,motor_speed_rpm,motor_torque)

        print("\nPERFORMANCE ANALYSIS")
        print("Inputs")
        print("Motor Speed (rpm):",motor_speed_rpm)
        print("Motor Torque (Nm):", motor_torque)
        print("Electrical Frequency (Hz):", self.electrical_frequency)
        print("Phase current (A_rms):", self.a_rms)
        print("Losses")
        print("Copper Loss (W):", self.copper_loss)
        print("Steel Loss (W):", self.steel.steel_loss)
        print("Total Loss (W):", self.total_loss)
        print("Power")
        print("Output power (kW):", self.output_power)
        print("Input power (kW):", self.input_power)
        print("Efficiency", self.efficiency)
        print("Thermal (steady-state):")
        print("Heat Flux (W/cm2):", self.heat_flux)
        print("Body temp (deg-C):", self.body_temperature)
        print("Winding temp (deg-C) (T):", self.winding_temperature)
        print("Assumed winding temp (T'):", self.assumed_winding_temp)
        print("Cogging Freq = f * S:", self.cogging_frequency)


    def display_thermal_resistance(self):

        self.calculate_thermal_resistance()

        print("\nTHERMAL RESISTANCE")
        print("Slot perimeter:", self.slot_perimeter)
        print("Conduction area (m2):", self.conduction_area)
        print("R_conduction:", self.R_conduction)
        print("Convection area (m2)",self.convection_area)
        print("R_convection", self.R_convection)



    def display_scalers(self):

        print("\nSCALERS")
        print("Scaler for strands (z):",self.z)
        print("Strands:",self.number_of_strands)
        print("Tooth depth d",self.tooth_depth)
        print()
        print("Scaler for radial L_ph (z2):",self.z2)
        print("L_ph (mm)", self.l_ph)
        print("A_ph (mm2)", self.a_ph)
        print("R_ph (mOhm)", self.r_ph)


    def display_all(self,motor_speed_rpm,motor_torque):

        self.display_scalers()
        self.display_geometry()
        self.display_cost()
        self.display_thermal_resistance()
        self.display_performance(motor_speed_rpm,motor_torque)


    def geometry_to_dict(self):


        geo_dict = {
            'rotor_od': self.rotor_od,
            'tooth_thickness': self.tooth_thickness,
            'slot_od': self.slot_od,
            'slot_depth': self.slot_depth,
            'yoke_thickness': self.yoke_thickness,
            'stator_od': self.stator_od,
            'rotor_id': self.rotor_id,
            'mean_slot_dia': self.mean_slot_dia,
            'stack_length': self.stack_length,
        }

        return geo_dict


    def cost_of_materials_to_dict(self):

        #calculate cost before calling this fn
        self.calculate_cost_of_materials()

        #mass in g, Cost in INR
        material_dict = {
            "stator_steel": [self.stator_steel_mass,self.stator_steel_cost],
            "rotor_steel": [self.rotor_steel_mass, self.rotor_steel_cost],
            "conductor": [self.conductor_mass, self.conductor_cost],
            "magnet": [self.magnet_mass, self.magnet_cost],
            "total": [self.total_mass, self.total_cost],
        }

        return material_dict



    def performance_analysis_to_dict(self,motor_speed_rpm,motor_torque):

        #ensure value of rph is fixed before running this -> motor.run
        self.performance_analysis(self.winding_temperature,motor_speed_rpm,motor_torque)

        performance_dict = {
            "motor_speed_rpm": float(motor_speed_rpm),
            "motor_torque": float(motor_torque),
            "electrical_frequency": float(self.electrical_frequency),
            "a_rms": float(self.a_rms),
            "copper_loss": float(self.copper_loss),
            "steel_loss": float(self.steel.steel_loss),
            "total_loss": float(self.total_loss),
            "output_power": float(self.output_power),
            "input_power": float(self.input_power),
            "efficiency":  float(self.efficiency),
            "heat_flux": float(self.heat_flux),
            "body_temperature": float(self.body_temperature),
            "winding_temperature": float(self.winding_temperature),
            "cogging_frequency": float(self.cogging_frequency),
        }

        return performance_dict
