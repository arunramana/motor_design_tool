from .topology.dxf.IPMSM_dxf import *
from .topology.dxf.SPMSM_dxf import *



class GenerateDXF:
    
    def __init__(self, topology, filename, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns):
        
        if(topology == ["IPMSM","RADIAL", "INTERIOR"]):
            
            self.run_dxf_object =  GenerateDXF_IPMSM(filename, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns)
        
        elif(topology == ["SPMSM","RADIAL", "SURFACE"]):
            
            self.run_dxf_object = GenerateDXF_SPMSM(filename, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns)

        else:     
            raise ValueError("run_femm: Invalid Topology")