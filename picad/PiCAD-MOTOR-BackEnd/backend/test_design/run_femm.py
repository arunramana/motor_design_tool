
from .topology.femm.IPMSM_femm import *
from .topology.femm.SPMSM_femm import *

#this class assigns to the correct topology

class GenerateFEMM:
    
    def __init__(self, topology, directory, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width, bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle,no_turns, stack_depth,u, v, w):
        
        
        if(topology == ["IPMSM","RADIAL", "INTERIOR"]):
            
            self.run_femm_object = GenerateFEMM_IPMSM(directory, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width, bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w)
        
        elif(topology == ["SPMSM","RADIAL", "SURFACE"]):
            
            self.run_femm_object = GenerateFEMM_SPMSM(directory, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width, bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w)

        else:
            
            raise ValueError("run_femm: Invalid Topology")
            
