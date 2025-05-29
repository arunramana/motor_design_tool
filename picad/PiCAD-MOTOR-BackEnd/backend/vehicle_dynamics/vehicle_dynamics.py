import numpy as np
import pandas as pd

import json
from .trig import *
from .vehicle_dynamics_methods import *

from scipy.optimize import fsolve
from scipy.optimize import minimize

import math
import warnings
warnings.filterwarnings("ignore")



class VehicleDynamics:

    def __init__(self,gvw,gear_ratio,wheel_radius,frontal_area,cd,rolling_resistance_coeff,gear_efficiency,rated_speed_kmph, continuous_gradient, max_speed, gradeability, time_to_cross_grade_from_rest, length_of_grade, acceleration_from_rest_to_speed, time_to_accelerate, v_dc):



        #Vehicle Parameters
        self.gvw = float(gvw)
        self.gear_ratio = float(gear_ratio)
        self.wheel_radius = float(wheel_radius)
        self.frontal_area = float(frontal_area)
        self.cd = float(cd)
        self.rolling_resistance_coeff = float(rolling_resistance_coeff)
        self.gear_efficiency = float(gear_efficiency)

        #Performance Specs
        self.rated_speed_kmph = float(rated_speed_kmph)
        self.continuous_gradient = float(continuous_gradient)
        self.max_speed = float(max_speed)
        self.gradeability = float(gradeability)
        self.time_to_cross_grade_from_rest = float(time_to_cross_grade_from_rest)
        self.length_of_grade = float(length_of_grade)
        self.acceleration_from_rest_to_speed = float(acceleration_from_rest_to_speed)
        self.time_to_accelerate = float(time_to_accelerate)


        #Voltage Input
        self.v_dc = float(v_dc)


        #important physical constants
        self.g = 9.8 #acc due to gravity m/s2
        self.rho_air = 1.2 #denisty of air kg/m3



    def run(self,verbose=True, factor=4, time_arr=[], speed_arr=[]):

        #display inputs
        if(verbose):
            self.display_inputs()

        #run output
        #self.motor_spec_output()
        self.motor_spec_output_modified(factor=factor)
        
        if(time_arr!=[] and speed_arr!=[]):
            
            self.watt_hour_per_km, total_distance, total_time = self.calculate_drive_cycle(time_arr, speed_arr)
            
            self.torque_equivalent, self.speed_equivalent, self.slope_equivalent =  self.calculate_equivalent_operating_points(self.watt_hour_per_km, total_distance*1000, total_time)
            
            self.speed_rpm_equivalent = (self.speed_equivalent/3.6)*(60/(2*np.pi))*self.gear_ratio/self.wheel_radius #rpm
            
            #print(self.speed_rpm_equivalent)
            #print(self.torque_equivalent, self.speed_equivalent, self.slope_equivalent)
            
        else:
            
            self.watt_hour_per_km,total_distance, total_time = 0,0,0


        #display output
        if(verbose):
            self.display_outputs()

        #display graphs



    def motor_spec_output(self,efficiency=80):
        
        

        self.motor_spec_rated_condition(slope=self.continuous_gradient)
        self.motor_spec_peak_speed()

       
        
        #at start for continuous -> after getting rated_condition()
        self.motor_speed_rpm_continuous_start = 0
        self.motor_torque_continuous_start = self.motor_torque_rated_speed
        self.motor_power_continuous_start = 0


        #no load speed for continuous
        self.motor_speed_rpm_no_load = self.motor_torque_peak_speed*(self.motor_speed_rpm_peak_speed-self.motor_speed_rpm_rated_speed)/(self.motor_torque_rated_speed-self.motor_torque_peak_speed)+self.motor_speed_rpm_peak_speed
        self.motor_torque_no_load = 0
        self.motor_power_no_load = 0


        self.motor_spec_peak_torque()
        
        #self.motor_spec_peak_power()
        
        self.motor_spec_peak_power_modified()
        
        

        #at start for peak -> after getting peak torque
        self.motor_speed_rpm_peak_start = 0
        self.motor_torque_peak_start = self.motor_torque_peak_torque
        self.motor_power_peak_start = 0


        #voltage specs
        

        self.a_dc_max = self.motor_power_peak_power*1000/(efficiency/100)/self.v_dc

        self.ke =self.v_dc/(2**0.5)/(self.motor_speed_rpm_no_load/1000)

        self.ke_rad = self.ke/(1000*2*np.pi/60)

        self.kt = (3**0.5)*self.ke_rad

        self.max_A_rms = self.motor_torque_peak_start/self.kt
        
        
    def motor_spec_output_modified(self,efficiency=80, factor=4):
        
        
        self.motor_spec_rated_condition(slope=self.continuous_gradient)
        self.motor_spec_peak_speed()

       
        
        #at start for continuous -> after getting rated_condition()
        self.motor_speed_rpm_continuous_start = 0
        self.motor_torque_continuous_start = self.motor_torque_rated_speed
        self.motor_power_continuous_start = 0


        #no load speed for continuous
        #self.motor_speed_rpm_no_load = self.motor_torque_peak_speed*(self.motor_speed_rpm_peak_speed-self.motor_speed_rpm_rated_speed)/(self.motor_torque_rated_speed-self.motor_torque_peak_speed)+self.motor_speed_rpm_peak_speed
        self.motor_torque_no_load = 0
        self.motor_power_no_load = 0


        self.motor_spec_peak_torque()
        
        #self.motor_spec_peak_power()
        
        self.motor_spec_peak_power_modified()
        
        

        #at start for peak -> after getting peak torque
        self.motor_speed_rpm_peak_start = 0
        self.motor_torque_peak_start = self.motor_torque_peak_torque
        self.motor_power_peak_start = 0

        
        
        
        self.A_speed = 0
        self.A_torque = self.motor_torque_peak_torque
        self.A_power = 0
        
        
        self.B_speed = self.motor_speed_rpm_peak_torque
        self.B_torque = self.motor_torque_peak_torque
        self.B_power = self.motor_power_peak_torque
        
        self.C_speed = factor*self.motor_speed_rpm_peak_torque
        self.C_torque = self.motor_torque_peak_torque
        self.C_power = self.C_speed*self.C_torque/9.55/1000
        
        self.D_speed = self.motor_speed_rpm_peak_speed
        self.D_torque = self.motor_torque_peak_speed
        self.D_power = self.D_speed*self.D_torque/9.55/1000
        
        m = (self.C_torque-self.D_torque)/(self.C_speed-self.D_speed)
        
        
        x = self.D_speed-self.D_torque/m
        
        self.E_speed = abs(x)
        self.E_torque = 0
        self.E_power = 0
        
        self.no_load_speed = self.E_speed
        self.motor_speed_rpm_no_load = self.E_speed
        
        
        #rated
        self.F_speed = 0
        self.F_torque = self.motor_torque_rated_speed
        self.F_power = 0
        
        self.G_speed = self.motor_speed_rpm_rated_speed
        self.G_torque = self.motor_torque_rated_speed
        self.G_power = self.G_speed * self.G_torque/9.55/1000
        
        
        c = self.D_torque-m*self.D_speed
        
        x = (self.G_torque-c)/m
        
        self.H_speed = abs(x)
        self.H_torque = self.G_torque
        self.H_power = self.H_speed * self.H_torque/9.55/1000
        
        
        
        #voltage specs
        
        self.motor_power_peak_power = self.C_power
        self.motor_torque_peak_start = self.A_torque
        
        self.a_dc_max = self.motor_power_peak_power*1000/(efficiency/100)/self.v_dc

        self.ke =self.v_dc/(2**0.5)/(self.motor_speed_rpm_no_load/1000)
        

        self.ke_rad = self.ke/(1000*2*np.pi/60)

        self.kt = (3**0.5)*self.ke_rad

        self.max_A_rms = self.motor_torque_peak_start/self.kt



    

    def motor_spec_rated_condition(self,slope=5):

        #check rated speed

        A = calculate_parameter_A(self.rho_air, self.cd, self.frontal_area, self.gvw, self.rolling_resistance_coeff)

        motor_torque = calc_motor_torque_rated_speed(self.rated_speed_kmph,A,self.g,self.rolling_resistance_coeff,self.gvw,self.wheel_radius,self.gear_ratio,self.gear_efficiency,slope)

        wheel_torque = calculate_wheel_torque(motor_torque,self.gear_ratio,self.gear_efficiency)

        B = calculate_parameter_B(wheel_torque,slope, self.gvw, self.wheel_radius, self.g,self.rolling_resistance_coeff)


        #find time taken for steady state = rated speed
        #with warnings.catch_warnings(action="ignore"):
        solution = fsolve(time_for_rated_speed, 0, args=(A,B,self.rated_speed_kmph),xtol=0.1)
        t = solution[0]


        speed = self.rated_speed_kmph/3.6


        motor_rpm = calculate_motor_rpm(speed,self.wheel_radius,self.gear_ratio)


        power = calculate_power(motor_torque,motor_rpm)


        self.motor_speed_rpm_rated_speed = round(motor_rpm)
        self.motor_torque_rated_speed = round(motor_torque)
        self.motor_power_rated_speed = round(power,2)



    def motor_spec_peak_speed(self,slope=0.001):

        A = calculate_parameter_A(self.rho_air, self.cd, self.frontal_area, self.gvw, self.rolling_resistance_coeff)

        motor_torque = calc_motor_torque_peak_speed(self.max_speed,A,self.g,self.rolling_resistance_coeff,self.gvw,self.wheel_radius,self.gear_ratio,self.gear_efficiency,slope)

        wheel_torque = calculate_wheel_torque(motor_torque,self.gear_ratio,self.gear_efficiency)

        B = calculate_parameter_B(wheel_torque,slope, self.gvw, self.wheel_radius, self.g,self.rolling_resistance_coeff)

        #find time taken for steady state = peak speed
        #with warnings.catch_warnings(action="ignore"):
        solution = fsolve(time_for_peak_speed, 0, args=(A,B,self.max_speed),xtol=0.1)
        t = solution[0]


        speed = self.max_speed/3.6


        motor_rpm = calculate_motor_rpm(speed,self.wheel_radius,self.gear_ratio)


        power = calculate_power(motor_torque,motor_rpm)


        self.motor_speed_rpm_peak_speed = round(motor_rpm)
        self.motor_torque_peak_speed = round(motor_torque)
        self.motor_power_peak_speed = round(power,2)


    def motor_spec_peak_torque(self):

    
        A = calculate_parameter_A(self.rho_air, self.cd, self.frontal_area, self.gvw, self.rolling_resistance_coeff)

        # Solve for B using fsolve
       
        #with warnings.catch_warnings(action="ignore"):
        solution = fsolve(distance_equation, 0.5, args=(A, self.time_to_cross_grade_from_rest, self.length_of_grade))
        
        B = solution[0]
        
        

        motor_torque = calc_motor_torque_from_B(B,self.g,self.rolling_resistance_coeff,self.gvw,self.wheel_radius,self.gear_ratio,self.gear_efficiency,self.gradeability)


        wheel_torque = calculate_wheel_torque(motor_torque,self.gear_ratio,self.gear_efficiency)

        #speed = self.length_of_grade/self.time_to_cross_grade_from_rest

        speed1 = calculate_speed(A,B,self.time_to_cross_grade_from_rest)

        motor_rpm = calculate_motor_rpm(speed1,self.wheel_radius,self.gear_ratio)

        power = calculate_power(motor_torque,motor_rpm)

        self.motor_speed_rpm_peak_torque = round(motor_rpm)
        self.motor_torque_peak_torque = round(motor_torque)
        self.motor_power_peak_torque = round(power,2)


    def motor_spec_peak_power(self):

        A = calculate_parameter_A(self.rho_air, self.cd, self.frontal_area, self.gvw, self.rolling_resistance_coeff)

        # Solve for B using fsolve
        #with warnings.catch_warnings(action="ignore"):
        solution = fsolve(velocity_equation, 1, args=(A, self.time_to_accelerate, self.acceleration_from_rest_to_speed))
        B = solution[0]
        


        motor_torque = calc_motor_torque_from_B(B,self.g,self.rolling_resistance_coeff,self.gvw,self.wheel_radius,self.gear_ratio,self.gear_efficiency,0)

        wheel_torque = calculate_wheel_torque(motor_torque,self.gear_ratio,self.gear_efficiency)

        speed = self.acceleration_from_rest_to_speed/3.6

        motor_rpm = calculate_motor_rpm(speed,self.wheel_radius,self.gear_ratio)

        power = calculate_power(motor_torque,motor_rpm)

        self.motor_speed_rpm_peak_power = round(motor_rpm)
        self.motor_torque_peak_power = round(motor_torque)
        self.motor_power_peak_power = round(power,2)
        
        
    
    def motor_spec_peak_power_modified(self):
        
        #assume peak power = power @ pk speed * 1.5
        
        factor = 4
        power = self.motor_power_peak_torque*factor
        
        
        
        
        #calculate time at which we reach peak power assuming peak torque
        A = calculate_parameter_A(self.rho_air, self.cd, self.frontal_area, self.gvw, self.rolling_resistance_coeff)
        
        wheel_torque = calculate_wheel_torque(self.motor_torque_peak_torque, self.gear_ratio,self.gear_efficiency)

        B = calculate_parameter_B(wheel_torque,0, self.gvw, self.wheel_radius, self.g,self.rolling_resistance_coeff)
        
        
        solution = fsolve(time_for_power, 1, args=(A, B, self.motor_torque_peak_torque, self.wheel_radius, self.gear_ratio, power),xtol=0.1)
        
        t1 = solution[0] 
        
        v1 = calculate_speed(A,B,t1)
        
        v2 = 100/3.6 #100kmph in m/s
        
        a1 = v1/t1
        
        c = a1*v1
        
        t2 = (v2**2-v1**2)/(2*c)
        
        t = t1+t2
        
        motor_rpm = calculate_motor_rpm(v1, self.wheel_radius, self.gear_ratio)
        
        self.motor_speed_rpm_peak_power = round(motor_rpm)
        self.motor_torque_peak_power = round(self.motor_torque_peak_torque)
        self.motor_power_peak_power = round(power,2)
        
        
    
    
    def calculate_drive_cycle(self,time_arr, speed_arr, regen=0):
        
        
        speed_prev = speed_arr[1]
        time_prev = time_arr[1]
        
        wh = 0
        distance_total = 0
        
        
        
        for i in range(2,len(time_arr)):
            
            speed = speed_arr[i]
            time = time_arr[i]
            
            duration = time-time_prev
            
            
            
            #print()
            if(speed_prev==0 and speed>0):
                
                #calculate 0 to speed
                wh_km, distance = self.calculate_wh_km_0_to_v(duration, speed)
                
                wh += wh_km*distance/1000
                distance_total += distance/1000
                #print(1, speed/3.6, duration, wh_km, distance)
                
            
            elif(speed_prev>0 and speed==0):
                
                #calculate speed to 0
                wh_km, distance = self.calculate_wh_km_0_to_v(duration, speed_prev)
                
                wh -= regen*wh_km*distance/1000
                distance_total += distance/1000
                #print(2, speed/3.6, duration, -wh_km, distance)
                
            
            elif(speed_prev < speed):
                
                wh_km, distance = self.calculate_wh_km_u_to_v(duration, speed_prev, speed)
                
                wh += wh_km*distance/1000
                distance_total += distance/1000
                
                #print(3, speed/3.6, duration, wh_km, distance)
                
                
            elif(speed < speed_prev):
                
                wh_km, distance = self.calculate_wh_km_u_to_v(duration, speed, speed_prev)
                
                wh -= regen*wh_km*distance/1000
                distance_total += distance/1000
                #print(4, speed/3.6, duration, -wh_km, distance)
                
            elif(speed == speed_prev):
                
                wh_km, distance = self.calculate_wh_km_constant(speed, duration)
                
                wh += wh_km*distance/1000
                distance_total += distance/1000
                #print(5, speed/3.6, duration, wh_km, distance)
                
            time_prev = time
            speed_prev= speed
            
                
        return wh/distance_total, distance_total, time_arr[-1]
                    
    
    def calculate_equivalent_operating_points(self, wh_km, total_distance, total_time):
        
        ss_speed = (total_distance/total_time)*3.6
        
        
        
        solution = fsolve(solve_wh_per_km_for_torque, 10, args=(self.gear_ratio, self.wheel_radius, wh_km))
        
        torque = solution[0]
        
        
        solution = fsolve(solve_ss_speed_for_slope, 0.1, args=(torque, self.rho_air, self.cd, self.frontal_area, self.rolling_resistance_coeff, self.gvw, self.wheel_radius, self.g, self.gear_ratio, self.gear_efficiency, ss_speed))

        slope = solution[0]
        
        #print((ss_speed/3.6)*(60/(2*np.pi))*self.gear_ratio/self.wheel_radius)
        
        return torque, ss_speed, slope
        
    def calculate_wh_km_0_to_v(self, time, speed):
        
        
        solution = fsolve(calculate_torque_zero_to_v, 10, args=(self.rho_air, self.cd, self.frontal_area, self.rolling_resistance_coeff, self.gvw, self.wheel_radius, self.g, self.gear_ratio, self.gear_efficiency, time, speed))
        
        torque = solution[0]
        
        wh_km = torque*self.gear_ratio/(3.6*self.wheel_radius)
        
        distance = ((speed/3.6)/2)*time
        
        return wh_km, distance
    

    def calculate_wh_km_u_to_v(self, time, u, v):
        
        u = u/3.6
        v = v/3.6
        
        
        solution = minimize(calculate_torque_u_to_v, 10, args=(self.rho_air, self.cd, self.frontal_area, self.rolling_resistance_coeff, self.gvw, self.wheel_radius, self.g, self.gear_ratio, self.gear_efficiency, time, u, v), method="Nelder-Mead")
        
        torque = solution.x[0]
        
        
        wh_km = torque*self.gear_ratio/(3.6*self.wheel_radius)
        
        distance = ((v+u)/2)*time
        
        
        return wh_km, distance
        
    
    def calculate_wh_km_constant(self, speed, time):
        
        speed = speed/3.6
        
        #calculate_torque_constant_speed(15.78,self.rho_air, self.cd, self.frontal_area, self.rolling_resistance_coeff, self.gvw, self.wheel_radius, self.g, self.gear_ratio, self.gear_efficiency, speed)
        
        solution = fsolve(calculate_torque_constant_speed, 10, args=(self.rho_air, self.cd, self.frontal_area, self.rolling_resistance_coeff, self.gvw, self.wheel_radius, self.g, self.gear_ratio, self.gear_efficiency, speed))
        
        torque = solution[0]
        
        
        wh_km = torque*self.gear_ratio/(3.6*self.wheel_radius)
        
        distance = speed*time
        
        return wh_km, distance
        
        

    def display_inputs(self):

        print("\nINPUTS:")
        print("\nVehicle Parameters (assumed)")
        print("gvw :",self.gvw)
        print("gear ratio :",self.gear_ratio)
        print("wheel radius :",self.wheel_radius)
        print("frontal area :",self.frontal_area)
        print("cd :",self.cd)
        print("rolling resistance coeff :",self.rolling_resistance_coeff)
        print("gear efficiency :",self.gear_efficiency)

        print("\nPerformance Specs")
        print("rated speed kmph :",self.rated_speed_kmph)
        print("rated speed kmph :",self.continuous_gradient)
        print("max speed :",self.max_speed)
        print("gradeability :",self.gradeability)
        print("time to cross grade from rest :",self.time_to_cross_grade_from_rest)
        print("length of grade :",self.length_of_grade)
        print("acceleration from rest to speed :",self.acceleration_from_rest_to_speed)
        print("time to accelerate :",self.time_to_accelerate)

        print("\nVoltage Input")
        print("v dc :",self.v_dc)

        print("\nPhysical Constants:")
        print("g:",self.g)
        print("rho_air:",self.rho_air)


    def display_outputs(self):

        #run motor spec outputs() before displaying outputs
        print("\nVehicle Dynamics Output:")


        print("\nContinuous Mode:\n")
        print("\t\trpm\tNm\tkW")
        print("At Start:\t{}\t{}\t{}".format(self.motor_speed_rpm_continuous_start,self.motor_torque_continuous_start, self.motor_power_continuous_start))
        print("Rated Speed:\t{}\t{}\t{}".format(self.motor_speed_rpm_rated_speed,self.motor_torque_rated_speed, self.motor_power_rated_speed))
        print("Peak Speed:\t{}\t{}\t{}".format(self.motor_speed_rpm_peak_speed,self.motor_torque_peak_speed, self.motor_power_peak_speed))
        print("No Load Speed:\t{}\t{}\t{}".format(self.motor_speed_rpm_no_load ,self.motor_torque_no_load,self.motor_power_peak_start))

        print("\nPeak Mode:\n")
        print("\t\trpm\tNm\tkW")
        print("At Start:\t{}\t{}\t{}".format(self.motor_speed_rpm_peak_start,self.motor_torque_peak_start, self.motor_power_peak_start))
        print("Peak Torque:\t{}\t{}\t{}".format(self.motor_speed_rpm_peak_torque,self.motor_torque_peak_torque, self.motor_power_peak_torque))
        print("Peak Power:\t{}\t{}\t{}".format(self.motor_speed_rpm_peak_power,self.motor_torque_peak_power, self.motor_power_peak_power))

        print("\nVoltage Specs:")
        print("A DC MAX:",self.a_dc_max)
        print("Ke:",self.ke)
        print("Ke in rad:", self.ke_rad)
        print("Kt:", self.kt)
        print("Max A RMS:", self.max_A_rms)
        
        print("\nDrive Cycle\n")
        print(self.watt_hour_per_km,"WH/km")
        


    def toJson(self):

        response_dict = {}

        response_dict["continuous"]={}

        #speed(rpm),  torque Nm, Power Kw
        response_dict["continuous"]["start"] = [self.motor_speed_rpm_continuous_start,self.motor_torque_continuous_start, self.motor_power_continuous_start]
        response_dict["continuous"]["rated_speed"] = [self.motor_speed_rpm_rated_speed,self.motor_torque_rated_speed, self.motor_power_rated_speed]
        response_dict["continuous"]["peak_speed"] = [self.motor_speed_rpm_peak_speed,self.motor_torque_peak_speed, self.motor_power_peak_speed]
        response_dict["continuous"]["no_load_speed"] = [self.motor_speed_rpm_no_load ,self.motor_torque_no_load,self.motor_power_peak_start]

        #speed(rpm),  torque Nm, Power Kw
        response_dict["peak"] = {}

        response_dict["peak"]["start"] = [self.motor_speed_rpm_peak_start,self.motor_torque_peak_start, self.motor_power_peak_start]
        response_dict["peak"]["peak_torque"] = [self.motor_speed_rpm_peak_torque,self.motor_torque_peak_torque, self.motor_power_peak_torque]
        response_dict["peak"]["peak_power"] = [self.motor_speed_rpm_peak_power,self.motor_torque_peak_power, self.motor_power_peak_power]

        response_dict["voltage_specs"] = {}

        response_dict["voltage_specs"]["v_dc"] = self.v_dc
        response_dict["voltage_specs"]["a_dc_max"] = self.a_dc_max
        response_dict["voltage_specs"]["ke"] = self.ke
        response_dict["voltage_specs"]["ke_rad"] = self.ke_rad
        response_dict["voltage_specs"]["kt"] = self.kt
        response_dict["voltage_specs"]["max_a_rms"] = self.max_A_rms

        #remove plot from here

        p1 = response_dict["continuous"]
        p2 = response_dict["peak"]

        #self.vehicle_plot(p1,"Continuous Mode","continuous_plot")

        #self.vehicle_plot(p2,"Peak Mode","peak_plot")

        response_json = json.dumps(response_dict, indent = 4)

        return response_json
    

    def toJson_modified(self):

        response_dict = {}

        response_dict["table"]={}
        
        response_dict["table"]["A"] = [self.A_speed, self.A_torque, self.A_power]
        
        response_dict["table"]["B"] = [self.B_speed, self.B_torque, self.B_power]
        
        response_dict["table"]["C"] = [self.C_speed, self.C_torque, self.C_power]
        
        response_dict["table"]["D"] = [self.D_speed, self.D_torque, self.D_power]
        
        response_dict["table"]["E"] = [self.E_speed, self.E_torque, self.E_power]
        
        
        
        
        response_dict["rated"]={}
        
                
        response_dict["rated"]["F"] = [self.F_speed, self.F_torque, self.F_power]
        
        response_dict["rated"]["G"] = [self.G_speed, self.G_torque, self.G_power]
        
        response_dict["rated"]["H"] = [self.H_speed, self.H_torque, self.H_power]
        
        
        


        response_dict["voltage_specs"] = {}

        response_dict["voltage_specs"]["v_dc"] = self.v_dc
        response_dict["voltage_specs"]["a_dc_max"] = self.a_dc_max
        response_dict["voltage_specs"]["ke"] = self.ke
        response_dict["voltage_specs"]["ke_rad"] = self.ke_rad
        response_dict["voltage_specs"]["kt"] = self.kt
        response_dict["voltage_specs"]["max_a_rms"] = self.max_A_rms
        
        response_dict["wh_km"] = self.watt_hour_per_km
        response_dict["torque_equivalent"] = self.torque_equivalent
        response_dict["speed_equivalent"] = self.speed_equivalent
        response_dict["speed_rpm_equivalent"] = self.speed_rpm_equivalent
        response_dict["slope_equivalent"] = self.slope_equivalent



        response_json = json.dumps(response_dict, indent = 4)

        return response_json



    def toDict_modified(self):

        response_dict = {}

        response_dict["table"]={}
        
        response_dict["table"]["A"] = [self.A_speed, self.A_torque, self.A_power]
        
        response_dict["table"]["B"] = [self.B_speed, self.B_torque, self.B_power]
        
        response_dict["table"]["C"] = [self.C_speed, self.C_torque, self.C_power]
        
        response_dict["table"]["D"] = [self.D_speed, self.D_torque, self.D_power]
        
        response_dict["table"]["E"] = [self.E_speed, self.E_torque, self.E_power]
        
        
        
        
        response_dict["rated"]={}
        
                
        response_dict["rated"]["F"] = [self.F_speed, self.F_torque, self.F_power]
        
        response_dict["rated"]["G"] = [self.G_speed, self.G_torque, self.G_power]
        
        response_dict["rated"]["H"] = [self.H_speed, self.H_torque, self.H_power]
        
        

        response_dict["voltage_specs"] = {}

        response_dict["voltage_specs"]["v_dc"] = self.v_dc
        response_dict["voltage_specs"]["a_dc_max"] = self.a_dc_max
        response_dict["voltage_specs"]["ke"] = self.ke
        response_dict["voltage_specs"]["ke_rad"] = self.ke_rad
        response_dict["voltage_specs"]["kt"] = self.kt
        response_dict["voltage_specs"]["max_a_rms"] = self.max_A_rms
        
        response_dict["wh_km"] = self.watt_hour_per_km
        
        response_dict["torque_equivalent"] = self.torque_equivalent
        response_dict["speed_equivalent"] = self.speed_equivalent
        response_dict["speed_rpm_equivalent"] = self.speed_rpm_equivalent
        response_dict["slope_equivalent"] = self.slope_equivalent


        return response_dict

    def toDict(self):

        response_dict = {}

        response_dict["continuous"]={}

        #speed(rpm),  torque Nm, Power Kw
        response_dict["continuous"]["start"] = [self.motor_speed_rpm_continuous_start,self.motor_torque_continuous_start, self.motor_power_continuous_start]
        response_dict["continuous"]["rated_speed"] = [self.motor_speed_rpm_rated_speed,self.motor_torque_rated_speed, self.motor_power_rated_speed]
        response_dict["continuous"]["peak_speed"] = [self.motor_speed_rpm_peak_speed,self.motor_torque_peak_speed, self.motor_power_peak_speed]
        response_dict["continuous"]["no_load_speed"] = [self.motor_speed_rpm_no_load ,self.motor_torque_no_load,self.motor_power_peak_start]

        #speed(rpm),  torque Nm, Power Kw
        response_dict["peak"] = {}

        response_dict["peak"]["start"] = [self.motor_speed_rpm_peak_start,self.motor_torque_peak_start, self.motor_power_peak_start]
        response_dict["peak"]["peak_torque"] = [self.motor_speed_rpm_peak_torque,self.motor_torque_peak_torque, self.motor_power_peak_torque]
        response_dict["peak"]["peak_power"] = [self.motor_speed_rpm_peak_power,self.motor_torque_peak_power, self.motor_power_peak_power]

        response_dict["voltage_specs"] = {}

        response_dict["voltage_specs"]["v_dc"] = self.v_dc
        response_dict["voltage_specs"]["a_dc_max"] = self.a_dc_max
        response_dict["voltage_specs"]["ke"] = self.ke
        response_dict["voltage_specs"]["ke_rad"] = self.ke_rad
        response_dict["voltage_specs"]["kt"] = self.kt
        response_dict["voltage_specs"]["max_a_rms"] = self.max_A_rms

        #remove plot from here

        p1 = response_dict["continuous"]
        p2 = response_dict["peak"]

        #self.vehicle_plot(p1,"Continuous Mode","continuous_plot")

        #self.vehicle_plot(p2,"Peak Mode","peak_plot")


        return response_dict

    '''
    def vehicle_plot(self,p1,title,fname):

        x,y1,y2 = [],[],[]

        for p in p1:

            arr = p1[p]

            x.append(arr[0])
            y1.append(arr[1])
            y2.append(arr[2])


        fig = plt.figure(figsize=(10,10))
        plt.plot(x,y1,"b",label="torque (Nm)",marker="o")
        plt.plot(x,y2,"orange", label="power (kw)", marker="o")

        plt.legend()


        plt.title(title)
        plt.xlabel("Motor Speed RPM")
        plt.ylabel("Power/Torque")

        plt.savefig("data/{}.jpg".format(fname),bbox_inches='tight')

    '''
