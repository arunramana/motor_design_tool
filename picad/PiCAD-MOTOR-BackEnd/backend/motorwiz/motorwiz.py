from .motorwiz_radial_inner_rotor import *

class MotorWiz:
    
    def __init__(self, peak_torque,max_no_load_rpm, voltage_limit, trv, aspect_ratio, no_of_poles, no_of_slots, coil_span, g_airgap, L_pm, tooth_height, fill_factor,  ambient_temperature, htc_cond_housing_stator, htc_cond_slotpaper, htc_conv, housing_stack_length_diff, fin_area_factor, magnet, steel, wire, topology):
        
        if(topology == ["IPMSM","RADIAL", "INTERIOR"] or topology == ["SPMSM","RADIAL", "SURFACE"]):
            
            self.motorwiz_object =  MotorWiz_Radial_InnerRotor(peak_torque,max_no_load_rpm, voltage_limit, trv, aspect_ratio, no_of_poles, no_of_slots, coil_span, g_airgap, L_pm, tooth_height, fill_factor,  ambient_temperature, htc_cond_housing_stator, htc_cond_slotpaper, htc_conv, housing_stack_length_diff, fin_area_factor, magnet, steel, wire, topology)
        
        else:     
            raise ValueError("run_femm: Invalid Topology")