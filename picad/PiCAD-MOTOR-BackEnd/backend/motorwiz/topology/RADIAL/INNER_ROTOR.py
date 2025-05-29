import numpy as np
import pandas as pd

from scipy.optimize import fsolve
import scipy
import math

import seaborn as sn 

from backend.motorwiz.materials.load_materials1 import *
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon

import warnings
warnings.filterwarnings("ignore")

import backend.motorwiz.topology.plots.motorOperatingParamsController as mo
import backend.motorwiz.topology.plots.motorFns as mtr

from ..plots.plot_curves import curves_to_dict

import json

import base64


class INNER_ROTOR:

    def __init__(self, magnet_position, magnet, steel, wire, thermal, r_g, tooth_thickness, tooth_depth_d, tooth_height, mean_slot_radius, stack_length, actual_slot_width_ratio, no_poles, no_slots, ambient_temperature, L_ph, A_ph, R_ph, Kt, housing_stack_length_diff):

        
        #materials
        self.magnet_position = magnet_position
        self.magnet = magnet
        self.steel = steel
        self.wire = wire

        self.thermal = thermal

        #variables
        self.r_g = r_g
        self.tooth_thickness = tooth_thickness
        self.tooth_depth_d = tooth_depth_d
        self.tooth_height = tooth_height
        #self.slot_outer_rad = slot_outer_rad
        self.mean_slot_radius = mean_slot_radius
        self.stack_length = stack_length
        self.actual_slot_width_ratio = actual_slot_width_ratio
        self.no_poles = no_poles
        self.no_slots = no_slots
        
        self.L_ph = L_ph
        self.A_ph = A_ph
        self.R_ph = R_ph
        self.ambient_temperature = ambient_temperature
        
        self.Kt = Kt
        
        self.housing_stack_length_diff = housing_stack_length_diff
        
        
        
       


    def create(self):
        
        #1) generate geometry params
        self.calculate_geometry()
        
        #2) calculate mass according to magnet position
        
        if(self.magnet_position=="INTERIOR"):
            self.calculate_mass_of_materials_interior()
        else:
            self.calculate_mass_of_materials_surface()
            

        #3) calculate cost 
        
        self.calculate_cost_of_materials()
        
        
        
    def run(self,project_name,i_max,v_dc,ld,lq,psi_m, max_no_load_rpm,mat_step = 30):
        
        
        curves_dct = curves_to_dict(i_max, v_dc, ld, lq, psi_m, self.no_poles/2,project_name, self.magnet_position, max_no_load_rpm)
        
        #get operating points
        operating_o = curves_dct["TORQUE SPEED"]["operating torque speed"][0]
        operating_t = curves_dct["TORQUE SPEED"]["operating torque speed"][1]
        
        
        τ_max = curves_dct["PEAK TORQUE"]
        max_no_load_rpm = curves_dct["MAX NL RPM"]
        
        #close polygon
        p = [[0,0]]
        for i in range(len(operating_o)):
        
            p.append([operating_o[i], operating_t[i]])
            
        poly = Polygon(p)
        
        
        #go through torque speed for efficiecy and winding temp map
        
        efficiency_matrix = np.ones((mat_step,mat_step))*50
        temperature_matrix = np.ones((mat_step,mat_step))*self.ambient_temperature
        
        
        tau_lim = operating_t[0]
        omega_lim = operating_o[-1]
        
        tau_step = tau_lim/mat_step
        omega_step = omega_lim/mat_step
        
        i_arr = np.zeros(mat_step)
        j_arr = np.zeros(mat_step)
        
        i=0
        i_mat = 0
        
        while i <= tau_lim+tau_step:
        
            if(i_mat>=mat_step):
                break
        
            j=0
            j_mat = 0
            
            flag = True
        
            
            while j <= omega_lim+omega_step:
        
                if(j_mat >= mat_step):
                    break
        
                
                if(poly.contains_point([j,i])):
                    
                    
                    self.assumed_winding_temp = 55
                    #solve winding temp
                    
                    rated_torque = i
                    rated_speed = j
                    
                    #with warnings.catch_warnings(action="ignore"):
                    solution = fsolve(self.solve_winding_temp, self.assumed_winding_temp, args=(rated_torque,rated_speed,τ_max))
                    self.winding_temperature = solution[0]
                        
                   
                    #re-calculate thermal
                    self.calculate_thermal(self.winding_temperature,rated_torque,rated_speed,τ_max)
                    
                    
                    if(self.efficiency*100 < 50):
                        self.efficiency = 50
                    elif(self.efficiency*100 > 100):
                        self.efficiency = 50
                    else:
                        efficiency_matrix[i_mat][j_mat] = self.efficiency*100
                    
                    
                    if self.winding_temperature > 150:
                        temperature_matrix[i_mat][j_mat] = 150
                        
                    elif self.winding_temperature < self.ambient_temperature:
                        temperature_matrix[i_mat][j_mat] = self.ambient_temperature
                    else:
                        temperature_matrix[i_mat][j_mat] = self.winding_temperature
                    
                    
                    i_arr[i_mat] = i
                    j_arr[j_mat] = j
                else:
                    
                    efficiency_matrix[i_mat][j_mat] = 0
                    temperature_matrix[i_mat][j_mat] = 0
                    
        
                j+=omega_step
                j_mat+=1
            
            i+=tau_step
            i_mat+=1
        
                
        #print(omega_lim,tau_lim)
        '''
        thermal_file = project_name+"/thermal_map.png"
        eff_file = project_name+"/efficiency_map.png"
        
        plt.figure(figsize=(10,10))
        #plt.contourf(j_arr,i_arr,efficiency_matrix,50,cmap=sn.color_palette("blend:#FFFFFF,#FF0000,#00FF00", as_cmap=True),vmin=0, corner_mask=False)
        plt.contourf(j_arr,i_arr,efficiency_matrix,50,cmap="jet",corner_mask=False)
        plt.colorbar()
        #plt.show()
        ax = plt.gca()
        ax.set_xlim([0, omega_lim])
        ax.set_ylim([0, tau_lim])
        plt.xlabel("Speed in RPM")
        plt.ylabel("Torque in Nm")
        
        plt.savefig(eff_file, bbox_inches='tight')
        plt.close()
        
        plt.figure(figsize=(10,10))
        #plt.contourf(j_arr,i_arr,temperature_matrix,50,cmap=sn.color_palette("blend:#FFFFFF,#0000FF,#FF0000", as_cmap=True),vmin=0, corner_mask=False)
        plt.contourf(j_arr,i_arr,temperature_matrix,50,cmap="jet",corner_mask=False)
        plt.colorbar()
        ax = plt.gca()
        #plt.show()
        ax.set_xlim([0, omega_lim])
        ax.set_ylim([0, tau_lim])
        plt.xlabel("Speed in RPM")
        plt.ylabel("Torque in Nm")
        
        plt.savefig(thermal_file, bbox_inches='tight')
        plt.close()
        '''
        
        dct = {}
        
        efficiency_matrix = scipy.ndimage.gaussian_filter(efficiency_matrix, sigma=0.6, mode="nearest")
        temperature_matrix = scipy.ndimage.gaussian_filter(temperature_matrix, sigma=0.6, mode="nearest")
        
        dct["data"] = {}
        dct["data"]["efficiency"] = {
            "x": j_arr.tolist(),
            "y": i_arr.tolist(),
            "z": efficiency_matrix.tolist(),
        }
        
        dct["data"]["thermal"] = {
            "x": j_arr.tolist(),
            "y": i_arr.tolist(),
            "z": temperature_matrix.tolist(),
        }
    
        '''
        with open(thermal_file, "rb") as image_file:
            t_string = base64.b64encode(image_file.read()).decode('utf8')
            
        with open(eff_file, "rb") as image_file:
            e_string = base64.b64encode(image_file.read()).decode('utf8')
            
        dct["images"] = {
            "efficiency": e_string,
            "thermal": t_string,
        }
        '''
        
        
        curves_dct["performance"] = dct
        
        return curves_dct
                
    
    def calc_nearest_point(self,omega_arr, tau_arr, x, y):
        
        dist = np.sqrt((x-omega_arr[0])**2+(y-tau_arr[0])**2)
        
        i1 = 0
        j1 = 0
        
        for i in range(len(tau_arr)):
            
            for j in range(len(omega_arr)):
                
                dist1 = np.sqrt((x-omega_arr[j])**2+(y-tau_arr[i])**2)
                
                if(dist1<dist):
                    dist = dist1
                    i1 = i
                    j1 = j
                    
        return i1, j1
                
                

    def calculate_geometry(self):
        
        self.slot_outer_rad  = self.r_g *(1+ self.tooth_depth_d)

        self.rotor_od = 2 * self.r_g
        
        #changes her for spmsm
        if(self.magnet_position=="SURFACE"):
            self.rotor_or_magnet = self.rotor_od/2
            self.rotor_or_lamination = self.rotor_or_magnet-self.magnet.w_m
        
        self.tooth_thickness = self.tooth_thickness
        self.slot_depth = self.tooth_height
        self.slot_od = 2 * self.slot_outer_rad
        self.mean_slot_dia = 2 * self.mean_slot_radius
        #self.stack_length = self.L
        self.yoke_thickness = (1 - self.actual_slot_width_ratio ) * np.pi * self.r_g / self.no_poles
        self.stator_od =2 *(self.slot_outer_rad + self.yoke_thickness)
        self.rotor_id=2*(self.r_g - self.yoke_thickness)
        
    

    def calculate_mass_of_materials_interior(self):
        
        self.stator_steel_mass = (np.pi *( self.stator_od - self.yoke_thickness )* self.yoke_thickness +self.tooth_thickness * self.slot_depth * self.no_slots )*self.stack_length *10**(-3) *self.steel.density

        self.rotor_steel_mass = (np.pi *( self.rotor_od ** 2 - self.rotor_id**2)/4-self.magnet.w_m * self.magnet.l_pm * self.no_poles)* self.stack_length *(10**(-3))*self.steel.density
        
        self.magnet_mass = self.magnet.w_m *self.magnet.l_pm *self.no_poles *self.stack_length *10**(-3)*self.magnet.density
        
        self.conductor_mass =3*self.L_ph*self.A_ph /1000*self.wire.density
        
        self.total_mass = self.stator_steel_mass + self.rotor_steel_mass + self.magnet_mass +self.conductor_mass


    def calculate_mass_of_materials_surface(self):
        
        self.stator_steel_mass = (np.pi *( self.stator_od - self.yoke_thickness )* self.yoke_thickness +self.tooth_thickness * self.slot_depth * self.no_slots )*self.stack_length *10**(-3) *self.steel.density

        #changes here when compared to interior permanent magnet
        self.rotor_steel_mass = (np.pi *( (self.rotor_or_lamination*2) ** 2 - self.rotor_id**2)/4)* self.stack_length *(10**(-3))*self.steel.density
        
        #pole span for spmsm
        pole_span = self.magnet.l_pm/self.rotor_or_lamination
        
        area_one_pole = pole_span*(self.rotor_or_magnet**2 - self.rotor_or_lamination**2)/2
        
        self.magnet_mass = area_one_pole *self.no_poles *self.stack_length *10**(-3)*self.magnet.density
        
        self.conductor_mass =3*self.L_ph*self.A_ph /1000*self.wire.density
        
        self.total_mass = self.stator_steel_mass + self.rotor_steel_mass + self.magnet_mass +self.conductor_mass
        



    def calculate_cost_of_materials(self):

        #cost calculation radial


        self.stator_steel_cost = self.stator_steel_mass*self.steel.rs_per_g



        self.magnet_cost = self.magnet_mass*self.magnet.rs_per_g


        self.rotor_steel_cost = self.rotor_steel_mass*self.steel.rs_per_g


        self.conductor_cost = self.conductor_mass*self.wire.rs_per_g


        self.total_cost = self.stator_steel_cost+self.rotor_steel_cost+self.magnet_cost+self.conductor_cost

    
    
    def calculate_thermal(self,assumed_winding_temp,rated_torque,rated_speed,τ_max):
        
        #calculate thermal after calculating mass
        
        self.op_current = rated_torque/self.Kt
        self.electrical_freq = rated_speed/ 60* self.no_poles /2
        self.j_rtd = self.op_current / self.A_ph
        self.p_out = rated_torque * rated_speed * 2 * np.pi /60
        
        
        self.cond_loss = 3 * (self.op_current**2) * (self.R_ph*1e-3)*(1+(assumed_winding_temp - 25)* self.wire.temp_coeff_of_resistivity/100)
        
                             
        self.steel.calculate_steel_loss(self.stator_steel_mass/1000, self.electrical_freq)
        
        self.steel.steelLoss = self.steel.steelLoss*self.steel.loss_factor*0.5*(1+rated_torque/τ_max)
        
        self.total_loss = self.cond_loss+self.steel.steelLoss
        self.p_in = self.p_out + self.total_loss
        self.efficiency = self.p_out/self.p_in
    
        self.slot_outer_rad  = self.r_g *(1+ self.tooth_depth_d)
        
        #after calculating geometry
        self.stator_or = self.stator_od/2
       
        self.A_cond = (2*np.pi*self.slot_outer_rad+(2*self.tooth_height - self.tooth_thickness)*self.no_slots)*self.stack_length*10**(-6)
        self.R_th_cond = 1/(self.thermal.htc_cond_slotpaper * self.A_cond)
        
        self.T_cond1 = self.total_loss*self.R_th_cond
        
        self.A_cond2 = 2*np.pi*self.stator_or*self.stack_length
        
        self.T_cond2 = self.cond_loss/(self.thermal.htc_cond_housing_stator * self.A_cond2)
        
        
        self.T_cond = self.T_cond1+self.T_cond2
        
        conv_area_factor = self.thermal.fin_area_factor * (self.housing_stack_length_diff+self.stack_length)/self.stack_length
        
        self.A_conv = 2* np.pi *self.stator_or*(self.stator_or+self.stack_length)*conv_area_factor*10**(-6)
        
        self.R_th_conv = 1/(self.thermal.htc_conv *self.A_conv)
        
        self.T_conv = self.total_loss * self.R_th_conv
        
        self.T_tot = self.T_cond + self.T_conv
        
        
        self.winding_temperature = self.T_tot+self.ambient_temperature
        

        
        
    def solve_winding_temp(self,assumed_winding_temp,rated_torque,rated_speed,τ_max):
        
        self.calculate_thermal(assumed_winding_temp, rated_torque, rated_speed, τ_max)
        
        return self.winding_temperature - assumed_winding_temp


    def display_IPMSM_geometry(self):

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
        
    def display_SPMSM_geometry(self):
        
        
        print("\nSPMSM RADIAL GEOMETRY")
        print("rotor_magnet_od:",self.rotor_or_magnet*2)
        print("rotor_lamination_od:",self.rotor_or_lamination*2)
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


    def display_thermal(self):
        #run this after solving winding temperature

        print("\nTHERMAL")
        print("outer rad", self.slot_outer_rad)
        print("A cond",self.A_cond)
        print("R th cond",self.R_th_cond)
        print("T cond",self.T_cond)
        print("A conv",self.A_conv)
        print("R th conv",self.R_th_conv)
        print("T conv",self.T_conv)
        print("T tot",self.T_tot)
        print("ambient",self.ambient_temperature)
        print("assumed_winding_temp",self.assumed_winding_temp)
        print("winding temperature",self.winding_temperature)
        
        
    def display_loss_and_efficiency(self):
        #run this after solving winding temperature
        print("\nLOSS AND EFFICIENCY")
        print("op current",self.op_current)
        print("electrical freq",self.electrical_freq)
        print("j rtd",self.j_rtd)
        print("cond loss",self.cond_loss)
        print("Eddy Loss",self.steel.eddyLoss)
        print("Hysterisis Loss",self.steel.hyst_loss)
        print("Extra Loss",self.steel.extraLoss)
        print("Steel Loss",self.steel.steelLoss)
        print("total loss",self.total_loss)
        print("p out",self.p_out)
        print("p in",self.p_in)
        print("efficiency",self.efficiency)

    
    def display_all(self):
        
        if(self.magnet_position=="INTERIOR"):
            self.display_IPMSM_geometry()
        else:
            self.display_SPMSM_geometry()
            
        
        self.display_cost()
        #self.display_loss_and_efficiency()
        #self.display_thermal()



    def geometry_to_dict(self):


        geo_dict = {
            'rotor_od': self.rotor_od,
            'tooth_thickness': self.tooth_thickness,
            'slot_od': self.slot_od,
            'slot_depth': self.slot_depth,
            'yoke_thickness': self.yoke_thickness,
            'stator_od': self.stator_od,
            'shaft_dia': self.rotor_od/4,
            'mean_slot_dia': self.mean_slot_dia,
            'stack_length': self.stack_length,
        }

        return geo_dict


    def cost_of_materials_to_dict(self):

        #calculate cost before calling this fn
        #self.calculate_cost_of_materials()

        #mass in g, Cost in INR
        material_dict = {
            "stator_steel": [self.stator_steel_mass,self.stator_steel_cost],
            "rotor_steel": [self.rotor_steel_mass, self.rotor_steel_cost],
            "conductor": [self.conductor_mass, self.conductor_cost],
            "magnet": [self.magnet_mass, self.magnet_cost],
            "total": [self.total_mass, self.total_cost],
        }

        return material_dict



    def loss_and_efficiency_to_dict(self):


        performance_dict = {
        	"op current": self.op_current,
        	"electrical freq": self.electrical_freq,
        	"j rtd": self.j_rtd,
        	"cond loss": self.cond_loss,
        	"Eddy Loss": self.steel.eddyLoss,
        	"Hysterisis Loss": self.steel.hyst_loss,
        	"Extra Loss": self.steel.extraLoss,
        	"Steel Loss": self.steel.steelLoss,
        	"total loss": self.total_loss,
        	"p out": self.p_out,
        	"p in": self.p_in,
        	"efficiency": self.efficiency,
        }

        return performance_dict
