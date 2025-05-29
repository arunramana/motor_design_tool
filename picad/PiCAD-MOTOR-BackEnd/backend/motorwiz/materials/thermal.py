import configparser



#change dir name to "" while testing from here

class Thermals:

    def __init__(self,fin_area_factor,housing_by_stacklength_ratio,ambient_temperature,external_htc,wire_type,steel_type,dir_name="backend/motorwiz/materials"):

        thermal_config = configparser.ConfigParser()
        thermal_config.read("{}/htc.ini".format(dir_name))

        if(type(fin_area_factor)==int or type(fin_area_factor)==float):
            self.fin_area_factor = fin_area_factor
        else:
            raise ValueError("Fin Area Factor should be a number")


        if(type(housing_by_stacklength_ratio)==int or type(housing_by_stacklength_ratio)==float):
            self.housing_by_stacklength_ratio = housing_by_stacklength_ratio
        else:
            raise ValueError("Housing/StackLength ratio should be a number")

        if(type(ambient_temperature)==int or type(ambient_temperature)==float):
            self.ambient_temperature = ambient_temperature
        else:
            raise ValueError("Ambient should be a number")


        if(type(ambient_temperature)==int or type(ambient_temperature)==float):
            self.htc = external_htc
        else:
            raise ValueError("External HTC should be a number")


        try:
            self.wire_steel_htc = float(thermal_config.get(str(wire_type).upper()+"-"+str(steel_type).upper(),"htc"))
        except:
            raise ValueError("HTC for given wire & steel not defined")



    def display(self):

        print(self.fin_area_factor,self.housing_by_stacklength_ratio,self.ambient_temperature,self.external_htc,self.wire_steel_htc)
