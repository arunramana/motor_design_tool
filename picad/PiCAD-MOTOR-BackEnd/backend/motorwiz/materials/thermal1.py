import configparser
import os



#change dir name to "" while testing from here

class Thermals:

    def __init__(self,fin_area_factor,ambient_temperature, htc_cond_housing_stator, htc_cond_slotpaper,htc_conv, delta_t_tot_max, wire_type,steel_type, dir_name="backend/motorwiz/materials/htc.ini"):

        material_dir = os.path.dirname(os.path.abspath(dir_name))
        material_file = os.path.join(material_dir, 'htc.ini')
        

        thermal_config = configparser.ConfigParser()
        thermal_config.read(material_file)
        


        if(type(fin_area_factor)==int or type(fin_area_factor)==float):
            self.fin_area_factor = fin_area_factor
        else:
            raise ValueError("Fin Area Factor should be a number")

        if(type(htc_conv)==int or type(htc_conv)==float):
            self.htc_conv = htc_conv
        else:
            raise ValueError("HTC conv should be a number")


        if(type(ambient_temperature)==int or type(ambient_temperature)==float):
            self.ambient_temperature = ambient_temperature
        else:
            raise ValueError("Ambient should be a number")
            
        
        

        if(type(htc_cond_housing_stator)==int or type(htc_cond_housing_stator)==float):
            self.htc_cond_housing_stator = htc_cond_housing_stator
        else:
            raise ValueError("HTC cond housing stator should be a number")


        if(type(htc_cond_slotpaper)==int or type(htc_cond_slotpaper)==float):
            self.htc_cond_slotpaper = htc_cond_slotpaper
        else:
            raise ValueError("HTC cond slotpaper should be a number")
            
        
        if(type(delta_t_tot_max)==int or type(delta_t_tot_max)==float):
            self.delta_t_tot_max = delta_t_tot_max
        else:
            raise ValueError("delta_t_tot_max should be a number")


        try:
            self.wire_steel_htc = float(thermal_config.get(str(wire_type).upper()+"-"+str(steel_type).upper(),"htc"))
        except:
            raise ValueError("HTC for given wire & steel not defined")



    def display(self):

        print('fin area factor:',self.fin_area_factor)
        print('ambient temperature:',self.ambient_temperature)
        print('htc cond slotpaper:',self.htc_cond_slotpaper)
        print('htc conv:',self.htc_conv)
        print('delta t tot max:',self.delta_t_tot_max)
        print('wire-steel htc',self.wire_steel_htc)

