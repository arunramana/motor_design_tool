import numpy as np
import pandas as pd

from .trig import *

from scipy.optimize import fsolve
import math



def calculate_parameter_A(rho_air, cd, frontal_area, gvw, rolling_resistance_coeff):

    return rho_air*cd*frontal_area/(2*gvw)


def calculate_parameter_B(wheel_torque,slope_deg, gvw, wheel_radius, g,rolling_resistance_coeff):

    slope_radians = slope_deg*np.pi/180

    return wheel_torque/(gvw*wheel_radius)-g*(rolling_resistance_coeff*cos(slope_deg)+sin(slope_deg))


def calculate_wheel_torque(motor_torque,gear_ratio,gear_efficiency):

    return (motor_torque*gear_ratio*gear_efficiency)/100



def calculate_steady_state_speed(A,B):
    #kmph
    return np.sqrt(B/A)*3.6

def calculate_speed(A,B,t):

    return np.sqrt(B/A)*(1-2/(np.exp(2*np.sqrt(A*B)*t)+1))

def calculate_motor_rpm(speed,wheel_radius,gear_ratio):

    return speed/(2*np.pi*wheel_radius)*gear_ratio*60

def calculate_power(motor_torque,motor_rpm):

    #in kilo watts

    return motor_torque*motor_rpm/9.55/1000

def calculate_distance(A,B,t):

    # in distance

    return (np.log(np.exp(2*np.sqrt(A*B)*t)+1)/A-np.sqrt(B/A)*t-np.log(2)/A)/1000


def calculate_friction_torque(rolling_resistance_coeff,gvw,g,slope_deg,wheel_radius):

    return rolling_resistance_coeff*gvw*g*cos(slope_deg)*wheel_radius


def calculate_gradient_torque(gvw,g,slope_deg,wheel_radius):

    return gvw*g*sin(slope_deg)*wheel_radius


def calculate_air_drag_torque(rho_air,cd,frontal_area,wheel_radius,A,B,t):

    speed = calculate_speed(A,B,t)
    return 1/2*rho_air*cd*frontal_area*(speed**2)*wheel_radius


def calculate_acceleration_torque(wheel_torque,friction_torque,gradient_torque,air_drag_torque):

    wheel_torque-(friction_torque+gradient_torque+air_drag_torque)


def calc_motor_torque_rated_speed(rated_speed,A,g,rolling_resistance_coeff,gvw,wheel_radius,gear_ratio,gear_efficiency,slope=5):

    #rated speed and slope are set slope = 5

    motor_torque = ((((rated_speed**2)*A/12.96)+g*(rolling_resistance_coeff*cos(slope)+sin(slope)))*gvw*wheel_radius*100)/(gear_ratio*gear_efficiency)
    return motor_torque

def calc_motor_torque_peak_speed(peak_speed,A,g,rolling_resistance_coeff,gvw,wheel_radius,gear_ratio,gear_efficiency,slope=0):

    #rated speed and slope are set slope = 5

    term1 = ((peak_speed**2)*A)/12.96

    term2 = g*(rolling_resistance_coeff*cos(slope) + sin(slope))

    term3 = (100*gvw*wheel_radius)/(gear_ratio*gear_efficiency)

    motor_torque = (term1+term2)*term3
    return motor_torque


def distance_equation(B, A, t, S):
    
    
    return (1/A)*np.log(np.exp(2*np.sqrt(A*B)*t)+1) - np.sqrt(B/A)*t - (np.log(2)/A) - S


def calc_motor_torque_from_B(B,g,rolling_resistance_coeff,gvw,wheel_radius,gear_ratio,gear_efficiency,slope):

    motor_torque = (B+g*(rolling_resistance_coeff*cos(slope)+sin(slope)))*((gvw*wheel_radius*100)/(gear_ratio*gear_efficiency))

    return motor_torque


def velocity_equation(B, A, t, V):
    
    return np.sqrt(B/A)*(1-2/(np.exp(2*np.sqrt(A*B)*t)+1))-V/3.6


def time_for_rated_speed(t,A,B,V):

    s=calculate_speed(A,B,t)*3.6

    return s-V


def time_for_peak_speed(t,A,B,V):

    s=calculate_speed(A,B,t)*3.6

    return s-V

def time_for_power(t,A,B,motor_torque,wheel_radius,gear_ratio,p):
    
    speed=calculate_speed(A,B,t) #m/s
    
    speed_rpm = calculate_motor_rpm(speed,wheel_radius,gear_ratio)
    
    power = calculate_power(motor_torque,speed_rpm)
    
    return power-p
    
    
def calculate_torque_zero_to_v(torque, rho_air, cd, frontal_area, rolling_resistance_coeff, gvw, wheel_radius, g, gear_ratio, gear_efficiency, t, speed):
    
    #pass speed in kmph

    slope_deg = 0
    
    
    
    wheel_torque = torque * gear_ratio * gear_efficiency/100
    
    A = calculate_parameter_A(rho_air, cd, frontal_area, gvw, rolling_resistance_coeff)
    
    B = calculate_parameter_B(wheel_torque,slope_deg, gvw, wheel_radius, g,rolling_resistance_coeff)
    
    
    speed_calc = calculate_speed(A,B,t)
    
    speed_ms = speed / 3.6
    
    
    return speed_calc - speed_ms


def calculate_torque_u_to_v(torque,rho_air, cd, frontal_area, rolling_resistance_coeff, gvw, wheel_radius, g, gear_ratio, gear_efficiency, t, u, v):
    
    slope_deg = 0
    
    
    
    wheel_torque = torque * gear_ratio * gear_efficiency/100
    
    A = calculate_parameter_A(rho_air, cd, frontal_area, gvw, rolling_resistance_coeff)
    
    B = calculate_parameter_B(wheel_torque,slope_deg, gvw, wheel_radius, g,rolling_resistance_coeff)
    
    if(A/B < 0):
        return 10
    
    y1 = u*np.sqrt(A/B)
    
    y2 = v*np.sqrt(A/B)
    
    z1 = np.log((1+y1)/(1-y1))
    
    z2 = np.log((1+y2)/(1-y2))
    
    t1 = z1/(2*np.sqrt(A*B))
    
    t2 = z2/(2*np.sqrt(A*B))
    
    
    return abs((t2-t1)-t)
    

def calculate_torque_constant_speed(torque,rho_air, cd, frontal_area, rolling_resistance_coeff, gvw, wheel_radius, g, gear_ratio, gear_efficiency, speed):
    
    slope_deg = 0
    
    
    
    wheel_torque = torque * gear_ratio * gear_efficiency/100
    
    A = calculate_parameter_A(rho_air, cd, frontal_area, gvw, rolling_resistance_coeff)
    
    B = calculate_parameter_B(wheel_torque,slope_deg, gvw, wheel_radius, g,rolling_resistance_coeff)
    
    if(B/A < 0):
        return 1000
    
    ss = np.sqrt(B/A)
    
    return ss-speed

def calculate_torque_constant_speed(torque,rho_air, cd, frontal_area, rolling_resistance_coeff, gvw, wheel_radius, g, gear_ratio, gear_efficiency, speed):
    
    slope_deg = 0
    
    
    
    wheel_torque = torque * gear_ratio * gear_efficiency/100
    
    A = calculate_parameter_A(rho_air, cd, frontal_area, gvw, rolling_resistance_coeff)
    
    B = calculate_parameter_B(wheel_torque,slope_deg, gvw, wheel_radius, g,rolling_resistance_coeff)
    
    if(B/A < 0):
        return 1000
    
    ss = np.sqrt(B/A)
    
    return ss-speed
    
    
def solve_wh_per_km_for_torque(motor_torque, gear_ratio, wheel_radius, wh_km):
    
    wh_km_calc = motor_torque*gear_ratio/(3.6*wheel_radius)
    
    return wh_km_calc-wh_km

def solve_ss_speed_for_slope(slope_deg, torque,rho_air, cd, frontal_area, rolling_resistance_coeff, gvw, wheel_radius, g, gear_ratio, gear_efficiency, ss_speed):
    
    wheel_torque = torque * gear_ratio * gear_efficiency/100
    
    A = calculate_parameter_A(rho_air, cd, frontal_area, gvw, rolling_resistance_coeff)
    
    B = calculate_parameter_B(wheel_torque,slope_deg, gvw, wheel_radius, g,rolling_resistance_coeff)
    
    if(B/A < 0):
        return -100
    
    ss_speed_calc = calculate_steady_state_speed(A,B)
    
    return ss_speed - ss_speed_calc
    