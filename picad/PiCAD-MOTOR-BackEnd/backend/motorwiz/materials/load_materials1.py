import configparser

import pandas as pd
import numpy as np

#print(config.get("M350-50A","B_sat"))

#when testing from this file change dir_name to ""

import os



class Magnet:

    def __init__(self,magnet_name,dir_name="backend/motorwiz/materials/magnet.ini"):
        
        
        
        material_dir = os.path.dirname(os.path.abspath(dir_name))
        material_file = os.path.join(material_dir, 'magnet.ini')
        
    
        magnet_config = configparser.ConfigParser()
        magnet_config.read(material_file)
        
        
        magnet_name = magnet_name.strip()
        
        try:
            self.b_r = float(magnet_config.get(magnet_name,"b_r"))
            self.k_pm = float(magnet_config.get(magnet_name,"k_pm"))
            self.mu_ra = float(magnet_config.get(magnet_name,"mu_ra"))
            self.density = float(magnet_config.get(magnet_name,"density"))
            self.rs_per_g = float(magnet_config.get(magnet_name,"rs_per_g"))

            self.l_pm = 0
            self.w_m = 0
        except Exception as e:
            
            raise ValueError("Magnet Name does not exist in Material Library")


    def set_dimensions(self,l_pm,w_m):

        self.l_pm = w_m
        self.w_m = l_pm


    def display(self):

        print("\nMagnet Properties:")
        print("Br:",self.b_r)
        print("L Pm:",self.l_pm)
        print("W m:", self.w_m)
        print("K Pm:",self.k_pm)
        print("Mu ra:",self.mu_ra)
        print("density:",self.density)
        print("rs per g:",self.rs_per_g)




class Steel:

    def __init__(self,steel_name,dir_name="backend/motorwiz/materials/steel.ini"):



        material_dir = os.path.dirname(os.path.abspath(dir_name))
        material_file = os.path.join(material_dir, 'steel.ini')
        

        steel_config = configparser.ConfigParser()
        steel_config.read(material_file)

        try:
            #config
            self.b_sat = float(steel_config.get(steel_name,"b_sat"))
            self.grade = float(steel_config.get(steel_name,"grade"))
            self.thickness = float(steel_config.get(steel_name,"thickness"))
            self.density = float(steel_config.get(steel_name,"density"))
            self.rs_per_g = float(steel_config.get(steel_name,"rs_per_g"))

            self.hyst_factor = float(steel_config.get(steel_name,"hyst_factor"))
            self.eddy_factor = float(steel_config.get(steel_name,"eddy_factor"))
            self.extra_factor = float(steel_config.get(steel_name,"extra_factor"))
            self.steel_loss = 0

        except:
            raise ValueError("Steel Name does not exist in Material Library")



        self.loss_factor = 2



    def interpolate(self,df,xval, ycol):
        return np.interp([xval], df["B"], df[ycol])[0]

    def calculate_steel_loss(self,stator_steel_mass, electrical_freq, steel_saturation=1.5, dir_name="backend/motorwiz/materials/steel_loss.xlsx"):
        
        #mass should be in kg
        
        material_dir = os.path.dirname(os.path.abspath(dir_name))
        material_file = os.path.join(material_dir, "steel_loss.xlsx")

        
        steel_loss_df = pd.read_excel(material_file)
        
        query = steel_loss_df[steel_loss_df["B"]==steel_saturation]
    
        if(len(query)==0):
            #interpolate steel loss
            coeffs = [self.interpolate(steel_loss_df,steel_saturation,"aj0"),self.interpolate(steel_loss_df,steel_saturation,"aj1"),self.interpolate(steel_loss_df,steel_saturation,"aj2")]
        else:
            coeffs = list(steel_loss_df[steel_loss_df["B"]==steel_saturation].iloc[0])
    
        hystCoeff = coeffs[1]
        eddyCoeff = coeffs[2]
        extraCoeff = coeffs[3]
    
        hystLoss = (hystCoeff * stator_steel_mass) * electrical_freq * self.grade / 360
        eddyLoss = eddyCoeff * stator_steel_mass * electrical_freq * electrical_freq * self.thickness / 0.5 * self.grade / 360
        extraLoss = extraCoeff * stator_steel_mass * (electrical_freq ** 1.5) * self.grade / 360
    
        steelLoss = hystLoss + eddyLoss + extraLoss
        
        del(steel_loss_df)
        
        self.hyst_loss = hystLoss
        self.eddyLoss = eddyLoss
        self.extraLoss = extraLoss
        self.steelLoss = steelLoss
    


    def display(self):

        print("\nSteel properties:")
        print("B sat:",self.b_sat)
        print("Grade:",self.grade)
        print("Thickness:",self.thickness)
        print("density:",self.density)
        print("RS Per Gram:",self.rs_per_g)
        
        print("\nSteel Loss Params:")
        print("Hysterisis factor:",self.hyst_factor)
        print("Eddy Factor:",self.eddy_factor)
        print("Extra Factor:",self.extra_factor)
      


class Wire:

        #1) fill factor -> user input
        #2) tooth_depth factor -> user input

    def __init__(self, gauge, wire_name,fill_factor,dir_name="backend/motorwiz/materials/wire.ini"):

        material_dir = os.path.dirname(os.path.abspath(dir_name))
        material_file = os.path.join(material_dir, 'wire.ini')
        swg_file = os.path.join(material_dir, 'SWG.csv')
        

        wire_config = configparser.ConfigParser()
        wire_config.read(material_file)
        
        try:
            #config
            
            if("SWG" in gauge):
                
                df = pd.read_csv(swg_file)
                
                val = int(gauge.split(" ")[0])
                dia = float(list(df[df["SWG"] == val]["Diameter"])[0])
                self.wire_diameter = dia
                del(df)
                
            else:
                raise ValueError("Invalid Gauge")
                
            
            self.resistivity = float(wire_config.get(wire_name,"resistivity"))
            self.temp_coeff_of_resistivity = float(wire_config.get(wire_name,"temp_coeff_of_resistivity"))
            self.density = float(wire_config.get(wire_name,"density"))
            self.rs_per_g = float(wire_config.get(wire_name,"rs_per_g"))
            
        except ValueError as e:
            
            raise ValueError("Load Material: {}".format(str(e)))

        except:
            raise ValueError("Wire Name does not exist in Material Library")


        if(type(fill_factor)== int or type(fill_factor)== float):
            self.fill_factor = fill_factor
        else:
            raise ValueError("Fill Factor should be a Number")






    def display(self):

        print("\nWire Properties")
        print("resitivity:",self.resistivity)
        print("temp coeff of res.:",self.temp_coeff_of_resistivity)
        print("density:",self.density)
        print("rs per g:",self.rs_per_g)
        print("wire dia:",self.wire_diameter)
        print("fill factor",self.fill_factor)




