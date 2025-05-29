import configparser


#print(config.get("M350-50A","B_sat"))

#when testing from this file change dir_name to ""


class Magnet:

    def __init__(self,magnet_name,dir_name="backend/motorwiz/materials"):

        magnet_config = configparser.ConfigParser()
        magnet_config.read("{}/magnet.ini".format(dir_name))

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

    def __init__(self,steel_name,dir_name="backend/motorwiz/materials"):



        steel_config = configparser.ConfigParser()
        steel_config.read("{}/steel.ini".format(dir_name))

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


    def calculate_steel_loss(self,speed,stator_weight,poles):

        loss_factor = 2

        stator_wt = stator_weight/1000

        pole_pairs = poles/2

        freq = (speed/60)*pole_pairs

        hyst_loss = self.hyst_factor*stator_wt*freq*(self.grade/360)

        eddy_loss = self.eddy_factor*stator_wt*(freq**2)*(self.thickness/0.5)*(self.grade/360)

        extra_loss = self.extra_factor*stator_wt*(freq**1.5)*(self.grade/360)

        steel_loss = (hyst_loss+eddy_loss+extra_loss)*loss_factor

        #print(speed, stator_wt, eddy_loss,hyst_loss,extra_loss,steel_loss)

        self.steel_loss = steel_loss



    def display(self):

        print("\nSteel properties:")
        print("B sat:",self.b_sat)
        print("Grade:",self.grade)
        print("Thickness:",self.thickness)
        print("density:",self.density)
        print("RS Per Gram:",self.rs_per_g)


class Wire:

        #1) fill factor -> user input
        #2) tooth_depth factor -> user input

    def __init__(self,wire_name,fill_factor,tooth_depth_factor,dir_name="backend/motorwiz/materials"):

        wire_config = configparser.ConfigParser()
        wire_config.read("{}/wire.ini".format(dir_name))

        try:
            #config
            self.wire_diameter = float(wire_config.get(wire_name,"wire_diameter"))
            self.resistivity = float(wire_config.get(wire_name,"resistivity"))
            self.temp_coeff_of_resistivity = float(wire_config.get(wire_name,"temp_coeff_of_resistivity"))
            self.density = float(wire_config.get(wire_name,"density"))
            self.rs_per_g = float(wire_config.get(wire_name,"rs_per_g"))

        except:
            raise ValueError("Wire Name does not exist in Material Library")


        if(type(fill_factor)== int or type(fill_factor)== float):
            self.fill_factor = fill_factor
        else:
            raise ValueError("Fill Factor should be a Number")


        if(type(tooth_depth_factor)== int or type(tooth_depth_factor)== float):
            self.tooth_depth_factor = tooth_depth_factor
        else:
            raise ValueError("Tooth Depth Factor should be a Number")





    def display(self):

        print("\nWire Properties")
        print("resitivity:",self.resistivity)
        print("temp coeff of res.:",self.temp_coeff_of_resistivity)
        print("density:",self.density)
        print("rs per g:",self.rs_per_g)
        print("wire dia:",self.wire_diameter)
        print("fill factor",self.fill_factor)
        print("tooth depth factor:",self.tooth_depth_factor)



#Magnet("N42").display()
