from backend.motorwiz.motorwiz import *
from backend.vehicle_dynamics.vehicle_dynamics import *

from backend.test_design.run_femm import *
from backend.test_design.dxf import *

import backend.database.database as db
import backend.database.otp_database as otp_db
import backend.database.otp_email as otp_email

import numpy as np
import pandas as pd
from datetime import datetime

from scipy.optimize import fsolve

from flask import Flask, redirect, url_for, request,send_from_directory
from flask import send_file


from flask_cors import CORS, cross_origin

import os
import io

from multiprocessing import Process

from waitress import serve 



import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import ezdxf
from ezdxf.addons.drawing import RenderContext, Frontend
from ezdxf.addons.drawing.properties import LayoutProperties
from ezdxf.addons.drawing.matplotlib import MatplotlibBackend
import re
import base64
import json

import shutil





class ServerException(BaseException):
    def __init__(self, m):
        self.message = m
    def __str__(self):
        return self.message


class DXF2IMG(object):

    default_img_format = '.png'
    default_img_res = 300
    def convert_dxf2img(self, names, img_format=default_img_format, img_res=default_img_res):
        for name in names:
            doc = ezdxf.readfile(name)
            msp = doc.modelspace()
            # Recommended: audit & repair DXF document before rendering
            auditor = doc.audit()
            # The auditor.errors attribute stores severe errors,
            # which *may* raise exceptions when rendering.
            if len(auditor.errors) != 0:
                raise ServerException("The DXF document is damaged and can't be converted!")
            else :
                fig = plt.figure()
                ax = fig.add_axes([0, 0, 1, 1])
                ctx = RenderContext(doc)
                msp_properties = LayoutProperties.from_layout(msp)
                msp_properties.set_colors("#ffffff")

                #ctx.current_layout_properties.background_color()
                out = MatplotlibBackend(ax)
                Frontend(ctx, out).draw_layout(msp, finalize=True, layout_properties=msp_properties)

                img_name = re.findall("(\S+)\.",name)  # select the image name that is the same as the dxf file name
                first_param = ''.join(img_name) + img_format  #concatenate list and string
                fig.savefig(first_param, dpi=img_res)
                plt.close()



def dxf_img(project_name):

    name = str(project_name)+'/motor_full.dxf'
    try:
        first =  DXF2IMG()
        b = first.convert_dxf2img([name],img_format='.png')
    except Exception as e:

        raise ServerException("dxf image: "+str(e))


def dxf_winding(project_name):

    name = str(project_name)+'/motor_winding.dxf'
    try:
        first =  DXF2IMG()
        b = first.convert_dxf2img([name],img_format='.png')
    except Exception as e:

        raise ServerException("winding dxf image: "+str(e))



def load_data_for_dxf(motor_wiz):


    no_slots = motor_wiz.no_slots

    no_poles = motor_wiz.no_poles



    stator_od = round(motor_wiz.motor.stator_od)

    rotor_od = round(motor_wiz.motor.rotor_od)

    stator_id = rotor_od + motor_wiz.g_airgap*2

    shaft_dia = round(motor_wiz.motor.rotor_od/4)

    slot_opening = 3 #motor_wiz.slot_opening_width

    tooth_tip_angle = 20 #ask sir -> advanced

    tooth_width = round(motor_wiz.motor.tooth_thickness)

    yoke_thickness = round(motor_wiz.motor.yoke_thickness)




    bridge_thickness = 1 #for now
    pole_length = round(motor_wiz.motor.magnet.l_pm)
    pole_width = round(motor_wiz.motor.magnet.w_m)

    #the three arrays below should be of same dimensions
    duct_distances = [0,0,0]
    duct_radii = [0,0,0]
    duct_angles = [0,0,0]

    notch_depth = 0
    notch_angle = 0


    no_turns = motor_wiz.Turns

    u = 0
    v = 0
    w = 0

    stack_depth = round(motor_wiz.stack_length)


    return no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width, bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w


def calibrate_motorwiz(topology,project_name,motor_wiz,motor_torque,motor_speed_rpm):

    try:
        #run motor wiz -> remember to remove verbose

        motor_wiz.run(verbose=False)

        no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w = load_data_for_dxf(motor_wiz)

        #run once to get bavg

        #generate dxf

        gdxfObj = GenerateDXF(topology, project_name, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width, bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns)
        
        gdxf = gdxfObj.run_dxf_object
        
        gdxf.draw_motor()

        #run femm

        rfObj = GenerateFEMM(topology, project_name, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w)
        
        rf = rfObj.run_femm_object

        
        rf.run_femm(verbose=False,saveplots=True)

        #iteration 1

        t_new = (1/rf.b_yoke_to_b_tooth_ratio)*motor_wiz.tooth_thickness

        w =1-(t_new*motor_wiz.no_slots)/(2*np.pi*motor_wiz.g_airgap)

        #update values


        motor_wiz.w = w
        motor_wiz.assumed_slot_width_ratio = w

        
        motor_wiz.B_g = rf.peak_airgap_flux

        #print("\n---------- second run -----------------\n")
        motor_wiz.run(b_g=motor_wiz.B_g ,updated_bavg=True,verbose=False)


        no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w = load_data_for_dxf(motor_wiz)

        #generate dxf 2
        gdxfObj = GenerateDXF(topology, project_name, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns)

        gdxf = gdxfObj.run_dxf_object
        
        gdxf.draw_motor()

        #run femm 2

        rfObj = GenerateFEMM(topology, project_name, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w)
        
        rf = rfObj.run_femm_object
        
        rf.run_femm(verbose=False,saveplots=True,prevSolution=True)
        
        

        return motor_wiz,rf

    except Exception as e:
        print("calibration"+str(e))
        raise(ServerException(str(e)))


def airgap_plot(filename):
    
    dct = {}
    dct1 = {}
    dct2 = {}
    
    df = pd.read_csv(filename+"/new_airgapflux.csv")
    length_arr = list(df["Length"])
    b_arr = list(df["B.n"])
    time_arr = list(df["time"])
    
    dct1 = {
       "time": time_arr,
       "b": b_arr,
       "length": length_arr,
    }
    
    df = pd.read_csv(filename+"/original_airgapflux.csv")
    length_arr = list(df["Length"])
    b_arr = list(df["B.n"])
    time_arr = list(df["time"])
    
    dct2 = {
       "time": time_arr,
       "b": b_arr,
       "length": length_arr,
    }
    
    dct = {
        "fundamental": dct1,
        "original": dct2,
    }
    
    
    del(df)

    return dct


def get_density_plot(project_name):
    #run converts density plot base64

    file_name = project_name+'/motor_density_plot.png'

    try:

        with open(file_name, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf8')

            #print("s='{}'".format(encoded_string))

        return encoded_string
    except Exception as e:
        raise(ServerException("dp plot:"+str(e)))


def get_airgap_plot(project_name):
    #run this after running first api

    try:
    
        agdict = airgap_plot(project_name)

        return agdict

    except Exception as e:
        raise(ServerException("ag plot:"+str(e)))


def get_dxf(project_name):
    #returns motor base64 image from dxf file

    try:
        dxf_img(project_name)

        file_name=str(project_name)+'/motor_full.png'

        with open(file_name, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf8')


        return encoded_string


    except Exception as e:

        raise(ServerException("get dxf:"+str(e)))



def get_winding(project_name):
    
    #run this after running first api


    try:
        dxf_winding(project_name)
        file_name=str(project_name)+'/motor_winding.png'

        with open(file_name, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf8')


        return encoded_string


    except Exception as e:

        raise(ServerException("get winding:"+str(e)))


        
def background_task(user_id, vehicle_name, motor_name, motor_wiz,topology,input_data, rated_speed=3000,rated_torque=20):
    
        
    try:
        
        dir_name = "data/{}".format(user_id.strip())
        
        if(not os.path.isdir(dir_name)):
            os.mkdir(dir_name)
            
        dir_name = dir_name+"/"+vehicle_name.strip()
        
        if(not os.path.isdir(dir_name)):
            os.mkdir(dir_name)
            
        
        dir_name = dir_name+"/"+motor_name.strip()
        
        if(not os.path.isdir(dir_name)):
            os.mkdir(dir_name)
            
        project_name = dir_name
        
        
        motor_wiz, rf = calibrate_motorwiz(topology,project_name,motor_wiz, rated_torque, rated_speed)
        

        
        
        currentPk = 100
        current_angle = -5
        
        i_max = motor_wiz.I_max
        i_max_pk = i_max*np.sqrt(2)
        v_dc = motor_wiz.V_DC/np.sqrt(2)
        τ_max = motor_wiz.τ_max
        
        max_no_load_rpm = motor_wiz.max_no_load_rpm*1.1
        
        #motor_wiz.display_all()
        
        pole_pairs = motor_wiz.motor.no_poles/2
        
        #print(pole_pairs)
        
        #imax is RMS, 
        #for Ld,LQ because we want psi M in pk, we give imax also in pk
        if(i_max/2 > 100):
            ld,lq,psi_m = rf.calculate_LDLQ(i_max/2, current_angle, pole_pairs)
            psi_m = abs(psi_m) 
        else:
            ld,lq,psi_m = rf.calculate_LDLQ(i_max/2, current_angle, pole_pairs)
            psi_m = abs(psi_m) #this is also RMS
        
        #print("HERE!!!!!!!!")
        #print(i_max,v_dc,ld,lq,psi_m,τ_max)
        
        #1257.886559502973 226.2741699796952 1.0598721697550722e-05 2.246350238340217e-05 0.0201600596813075 270
        if(ld>lq):
            t = ld
            ld = lq
            lq = t
        
        
        plot_dct = motor_wiz.motor.run(project_name,i_max,v_dc,ld,lq,psi_m, max_no_load_rpm)
        
        
        response_dict = {}
        
        response_dict["error"] = "none"
        response_dict["status"] = "finished"

        response_dict["geometry"] = motor_wiz.motor.geometry_to_dict()


        response_dict["cost"] = motor_wiz.motor.cost_of_materials_to_dict()


        response_dict["performance"] = motor_wiz.motor.loss_and_efficiency_to_dict()

        response_dict["B_avg"] = float(motor_wiz.B_av)
        
        
        response_dict["ldlq"] = [ld*1e6, lq*1e6, psi_m]

        response_dict["turns"] = round(float(motor_wiz.Turns))

        response_dict["strands"] = round(float(motor_wiz.strands))
        
        response_dict["rph"] = round(float(motor_wiz.motor.R_ph))
        
        #ke = (motor_wiz.motor.no_poles*1000/9.55/2)*psi_m #back emf at 1k RPM
        #kt = (ke/104.72)*1.732
        
        ke = motor_wiz.Ke
        kt = motor_wiz.Kt
        
            
        
        
        response_dict["bemf"] = ke 
        
        response_dict["kt"] = kt
        
        
        response_dict["i_max"] = round(float(i_max))
        
        response_dict["max_no_load_speed"] = round(float(motor_wiz.max_no_load_rpm))
        
        #plots

        response_dict["airgap_plot"] = get_airgap_plot(project_name)

        #base64
        response_dict["density_plot"] = get_density_plot(project_name)

        response_dict["dxf"] = get_dxf(project_name)

        response_dict["winding"] = get_winding(project_name)
        
        response_dict["plots"] = plot_dct

        response_json = json.dumps(response_dict)
        

        #save response
        #json.parse from front end
        msg, status = db.update_motor(user_id, vehicle_name, motor_name, input_data, response_json, "none", "finished")
        
        if(status!=200):
            
            raise ServerException(msg)
            
        else:
            
            #p.join()
            
            shutil.make_archive(dir_name+"/download", 'zip', dir_name+"/download")
            
            file_path = dir_name+"/download.zip"
            msg, status = db.store_zip_file(user_id, vehicle_name, motor_name, file_path)
            
            if(status != 200):
                raise ServerException(msg)

            return response_json

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Server Error: {}".format(str(e))
        response_dict["status"] = "error"

        
        response_json = json.dumps(response_dict)

        msg, status = db.update_motor(user_id, vehicle_name, motor_name, input_data, response_json, "error", "error")
        
        if(status!=200):
            
            raise ValueError(msg)
            
        else:
            
            #p.join()

            return response_json
        
        #p.join()
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = str(e)
    
        response_dict["status"] = "error"

        response_json = json.dumps(response_dict, indent = 4)
        
        #p.join()
        
        return response_json
    


app = Flask(__name__)

cors = CORS(app, resources={r"/create_motor_version": {"origins": "*"},r"/send_email_forgot_password": {"origins": "*"}, r"/verify_otp": {"origins": "*"}, r"/reset_password": {"origins": "*"},r"/application_dynamics": {"origins": "*"}, r"/create_application": {"origins": "*"}, r"/get_all_applications": {"origins": "*"}, r"/get_application_data": {"origins": "*"}, r"/delete_application": {"origins": "*"}, r"/signup": {"origins": "*"},r"/signup_verify_otp": {"origins": "*"}, r"/login": {"origins": "*"}, r"/download": {"origins": "*"}, r"/": {"origins": "*"}, r"/vehicle_dynamics": {"origins": "*"}, r"/motorwiz": {"origins": "*"}, r"/get_all_vehicles": {"origins": "*"}, r"/create_vehicle": {"origins": "*"}, r"/get_vehicle_data": {"origins": "*"},  r"/get_all_motors": {"origins": "*"},  r"/create_motor": {"origins": "*"}, r"/get_motor_data": {"origins": "*"}, r"/create_vehicle_version": {"origins": "*"}, r"/delete_vehicle": {"origins": "*"}, r"/delete_motor": {"origins": "*"}, r"/delete_all_motors": {"origins": "*"},r"/update_profile": {"origins": "*"},r"/get_profile": {"origins": "*"},r"/change_password":{"origin":"*"},r"/get_userdata":{"origin":"*"},r"/update_feedback":{"origin":"*"},r"/update_all_feedback":{"origin":"*"}},support_credentials=True)
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route('/login', methods=['POST'])

def login():
    
    try:
        
        data = request.json
        
        email = data.get("email")
        password = data.get("password")
        
        msg, status = db.login(email, password)
        
        if(status != 200):
            
            raise ServerException(msg)
            
        
        response_dict = {}
        response_dict["error"] = "none"
        response_dict["user_id"] = msg[0]
        response_dict["email"] = msg[1]
        response_dict["user_name"] = msg[2]
        response_dict["valid"] = msg[3]
        

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
        
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Login Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Login Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json


@app.route('/signup', methods=['POST'])
def signup():
    try:
        data = request.json
        email = data.get("email")
        password = data.get("password")
        username = data.get("user_name")
        about = data.get("about")
        location = data.get("location")

        # Check if email already exists in the main database
        msg, status = db.get_email(email)
        if status == 200:
            raise ServerException("Email is already registered")

        otp = otp_db.generate_otp()
        msg, status = otp_db.store_otp1(email, otp)
        
        if status != 200:
            raise ServerException(msg)
        
        msg, status = otp_email.signup_email_otp(email, otp)
        
        if status != 200:
            raise ServerException(msg)

        response_dict = {"message": msg, "error": "none"}
        return json.dumps(response_dict, indent=4)
        
    except ServerException as e:
        response_dict = {"error": f"Sign Up Error: {str(e)}"}
        return json.dumps(response_dict, indent=4)
    
    except ValueError as e:
        response_dict = {"error": f"Sign Up Error: {str(e)}"}
        return json.dumps(response_dict, indent=4)


@app.route('/signup_verify_otp', methods=['POST'])
def signup_verify_otp():
    try:
        data = request.json
        email = data.get("email")
        otp = int(data.get("otp"))
        username = data.get("user_name")
        password = data.get("password")
        about = data.get("about_me")
        location = data.get("location")

        if not otp_db.email_exists_in_otp(email):
            raise ServerException("Email not found in OTP database")
        
        msg, status = otp_db.verify_otp1(email, otp)
        
        if status != 200:
            raise ServerException(msg)

        # Check if email is in the whitelist table
        whitelist_msg, whitelist_status  = db.check_whitelist(email)

        valid = whitelist_status == 200

        # Create the user in the main database
        msg, status = db.create_user(username, email, password, about=about, location=location, valid=valid)
        
        if status != 200:
            raise ServerException(msg)
        
        response_dict = {
            "error": "none",
            "user_id": msg[0],
            "email": msg[1],
            "user_name": msg[2],
            "valid": valid
        }
        return json.dumps(response_dict, indent=4)

    except ServerException as e:
        response_dict = {"error": f"Verify OTP Server Error: {str(e)}"}
        return json.dumps(response_dict, indent=4)
    
    except ValueError as e:
        response_dict = {"error": f"Verify OTP Value Error: {str(e)}"}
        return json.dumps(response_dict, indent=4)

@app.route('/vehicle_dynamics', methods=['POST'])

def vehicle_dynamics_api():

    try:
        #POST all the variables for vehicle dynamics

        data = request.json

        #email = data.get("email")
        user_id = data.get("user_id").strip()
        vehicle_name = data.get("vehicle_name").strip()
        gvw = data.get("gvw")
        gear_ratio = data.get("gear_ratio")
        wheel_radius = data.get("wheel_radius")
        frontal_area = data.get("frontal_area")
        cd = data.get("cd")
        rolling_resistance_coeff = data.get("rolling_resistance_coeff")
        gear_efficiency = data.get("gear_efficiency")
        rated_speed_kmph = data.get("rated_speed_kmph")
        continuous_gradient = data.get("continuous_gradient")
        max_speed = data.get("max_speed")
        gradeability = data.get("gradeability")
        time_to_cross_grade_from_rest = data.get("time_to_cross_grade_from_rest")
        length_of_grade = data.get("length_of_grade")
        acceleration_from_rest_to_speed = data.get("acceleration_from_rest_to_speed")
        time_to_accelerate = data.get("time_to_accelerate")
        v_dc = data.get("v_dc")
        
        factor = float(data.get("factor"))

        #load data
        

        drive_cycle = data['drive_cycle']

        if(type(drive_cycle)==str):
            drive_cycle = json.loads(drive_cycle)


        time_arr = []
        speed_arr = []

        for dct in drive_cycle:

            time_arr.append(dct["x"])
            speed_arr.append(dct["y"])

        

        vehicle = VehicleDynamics(gvw,gear_ratio,wheel_radius,frontal_area,cd,rolling_resistance_coeff,gear_efficiency,rated_speed_kmph,continuous_gradient, max_speed, gradeability, time_to_cross_grade_from_rest, length_of_grade, acceleration_from_rest_to_speed, time_to_accelerate, v_dc)


        #run -> remember to turn off verbose
        vehicle.run(verbose=False, factor=factor, time_arr=time_arr, speed_arr=speed_arr)
        
        
            
        input_data = data
        output_data = vehicle.toDict_modified()
        error = "none"
        desc = vehicle_name
        
        if(not db.is_existing_vehicle(user_id, vehicle_name)):
            
            db.create_vehicle(user_id, vehicle_name, input_data, output_data, error, desc)
            
        else:
            
            db.update_vehicle(user_id, vehicle_name, input_data, output_data, error)

        return vehicle.toJson_modified()

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Server Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Value Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json


@app.route('/motorwiz', methods=['POST'])
def motorwiz_api():
    
    global p

    try:

        #POST all the variables for motorwiz

        data = request.json

        #email = data.get("email")
        user_id = data.get("user_id").strip()
        vehicle_name = data.get("vehicle_name").strip()
        motor_name = data.get("motor_name").strip()
        

        peak_torque = data.get("peak_torque")
        
        rated_speed = data.get("rated_speed")
        rated_torque = data.get("rated_torque")
        
        max_no_load_rpm = data.get("max_no_load_rpm")
        voltage_limit = data.get("voltage_limit")
        trv = data.get("trv")
        aspect_ratio = data.get("aspect_ratio")
        no_of_poles = data.get("no_of_poles")
        no_of_slots = data.get("no_of_slots")
        coil_span = float(data.get("coil_span"))
        g_airgap = float(data.get("g_airgap"))
        L_pm = float(data.get("L_pm"))
        tooth_height = float(data.get("tooth_height"))
        fill_factor = float(data.get("fill_factor"))/100
        ambient_temperature = float(data.get("ambient_temperature"))
        htc_cond_yoke_housing = float(data.get("htc_cond_yoke_housing"))
        htc_cond_wire_yoke = 200.0#data.get("htc_cond_wire_yoke")
        htc_conv = float(data.get("htc_conv"))
        fin_area_factor = float(data.get("fin_area_factor"))
        housing_stack_length_diff= float(data.get("housing_stack_length_diff"))
        magnet = data.get("magnet")
        steel = data.get("steel")
        wire = data.get("wire")
        topology = data.get("topology")
        
        application_type = data.get("application_type")
        
        
        
        htc_cond_slotpaper = htc_cond_yoke_housing
        htc_cond_housing_stator = htc_cond_wire_yoke
        
        #conv_area_factor = 10.5
        #htc_cond_slotpaper = 200
        
        if(application_type=="vehicle"):
            
                    
            if(not db.is_existing_vehicle(user_id, vehicle_name)):
                
                raise ServerException("Vehicle does not exist")
            
        else:

                        
            if(not db.is_existing_application(user_id, vehicle_name)):
                
                raise ServerException("Application does not exist")
            
        
        if(not db.is_existing_motor(user_id, vehicle_name, motor_name)):
            
            raise ServerException("Motor does not exist")
            

    
        input_data = json.dumps(data)
        output_data = json.dumps([])
            
            
        msg, status =  db.update_motor(user_id, vehicle_name, motor_name, input_data, output_data, "none", "started")
        
        if(status!=200):
            
            raise ServerException(msg)
        
        motorwizObj = MotorWiz(peak_torque,max_no_load_rpm, voltage_limit, trv, aspect_ratio, no_of_poles, no_of_slots, coil_span, g_airgap, L_pm, tooth_height, fill_factor,  ambient_temperature, htc_cond_housing_stator, htc_cond_slotpaper, htc_conv, housing_stack_length_diff, fin_area_factor, magnet, steel, wire, topology)

        motor_wiz = motorwizObj.motorwiz_object


        p = Process(target=background_task, args=(user_id, vehicle_name, motor_name, motor_wiz, topology, input_data))
        p.start()
        
        response_dict = {}
        response_dict["status"] = "started"
        response_dict["error"] = "none"
        

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
        
       

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Motorwiz Error: {}".format(str(e))
        response_dict["status"] = "error"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["status"] = "error"
        response_dict["error"] = str(e)

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
        

@app.route('/get_motor_data', methods=['POST'])

def simulation_results():
    

    try:
        data = request.json
    
        user_id = data.get("user_id").strip()
        #email = data.get("email")
        
        vehicle_name = data.get("vehicle_name").strip()
        motor_name = data.get("motor_name").strip()
        
   
        motor_data, status = db.get_motor_data(user_id, vehicle_name, motor_name)
        
        if(status != 200):
            
            raise ServerException(motor_data)
            
            
        
        input_data = motor_data["input"]
        output_data = motor_data["output"]
        
        if(type(output_data)==str):
            output_data = json.loads(output_data)
        
        if(output_data!=[]):
            
            
                        
            response_dict = {}
            response_dict["input"] = input_data
            response_dict["output"] = output_data 
            response_dict["error"] = "none"
            response_dict["status"] = "finished"

            

            response_json = json.dumps(response_dict, indent = 4)
            return response_json
         
    
        
        else:
            
            response_dict = {}
            response_dict["status"] = "running"
            response_dict["error"] = "none"
    
            response_json = json.dumps(response_dict, indent = 4)
            return response_json
        
        
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "simulation results: {}".format(str(e))
        response_dict["status"] = "error"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["status"] = "error"
        response_dict["error"] = "simulation results: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
    

@app.route('/create_vehicle', methods=['POST'])
def create_vehicle():
    
    #os version no db
    try:
        
                
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        vehicle_name = data.get("vehicle_name").strip()
        desc = data.get("desc")
        
        if(desc==None or desc==""):
            desc = vehicle_name
            
        if(db.is_existing_vehicle(user_id, vehicle_name)):
            
            raise ServerException("Vehicle Name Already Exists")
        
        error = "none"

        created_at = datetime.now().strftime("%d/%m/%Y %I:%M %p")
        
        msg, status = db.create_vehicle(user_id, vehicle_name, [], [], error, desc,created_at)
        
        if(status != 200):
            
            raise ServerException(msg)
        
        response_dict = {}
        response_dict["error"] = "none"
        
        response_json = json.dumps(response_dict, indent = 4)
        return response_json
            
        
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "create Vehicles: {}".format(str(e))
        response_dict["status"] = "error"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["status"] = "error"
        response_dict["error"] = "create Vehicles: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
    


@app.route('/get_all_vehicles', methods=['POST'])
def get_vehicles():
    
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        
        msg, status = db.get_all_vehicles(user_id)
        
        if(status!=200):
            
            raise ServerException(msg)
            
        else:
            
            if(type(msg)==str):
                vehicles = json.loads(msg)
            else:
                vehicles=msg
            
            if(vehicles==[]):
                
                response_dict = {}
                response_dict["vehicle_names"] = []
                response_dict["vehicle_desc"] = []
                response_dict["created_at"] = [] 
                response_dict["error"] = "none"
                

                response_json = json.dumps(response_dict, indent = 4)
                return response_json
            
            else:
                
                vehicle_names = []
                desc_lst = []
                create_lst =[]
                
                for v in vehicles:
                    
                    vehicle_names.append(v["vehicle_name"])
                    desc_lst.append(v["desc"])
                    create_lst.append(v["created_at"])
                    
                response_dict = {}
                response_dict["vehicle_names"] = vehicle_names
                response_dict["vehicle_desc"] = desc_lst
                response_dict["created_at"] = create_lst
                response_dict["error"] = "none"
                
    
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
        
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Get All Vehicles: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Get All Vehicles: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        
@app.route('/get_vehicle_data', methods=['POST'])
def get_vehicle_data():
    
    #os version no db
    
    try:
        
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        
        
        vehicle_data, status = db.get_vehicle_data(user_id, vehicle_name)
        
        if(status != 200):
            
            raise ServerException(vehicle_data)
            
        else:
            
            input_data = vehicle_data["input"]
            output_data = vehicle_data["output"]
            
            
            if(input_data==[] or input_data==None or output_data==[] or output_data==None):
                
                #for new vehicle
                
                response_dict = {}
                response_dict["input"] = []
                response_dict["output"] = []
                response_dict["error"] = "none"
                
    
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
            
            else:
                
                response_dict = {}
                response_dict["input"] = input_data
                response_dict["output"] = output_data 
                response_dict["error"] = "none"
                
    
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
                
    
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Get Vehicles data: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Get Vehicles data: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)



@app.route('/create_motor', methods=['POST'])
def create_motor():

    #os version no db
    try:
        
                
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        motor_name = data.get("motor_name").strip()
        desc = data.get("desc")
        
        if(desc==None or desc==""):
            desc = motor_name

        
            
        if(db.is_existing_motor(user_id, vehicle_name, motor_name)):
            
            raise ServerException("Motor Name Already Exists")
        created_at = datetime.now().strftime("%d/%m/%Y %I:%M %p")
        
        motor_data, status = db.create_motor(user_id, vehicle_name, motor_name, [], [], "none", desc, "",created_at)
            
        if(status != 200):
            
            raise ServerException(motor_data)
            
        response_dict = {}
        response_dict["error"] = "none"
        
        response_json = json.dumps(response_dict, indent = 4)
        return response_json
            
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "create motor: {}".format(str(e))
        response_dict["status"] = "error"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["status"] = "error"
        response_dict["error"] = "create motor: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
    


@app.route('/get_all_motors', methods=['POST'])
def get_motors():
    
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        vehicle_name = data.get("vehicle_name").strip()
        
        
        motor_lst, status = db.get_all_motors(user_id, vehicle_name)
        
        if(status!=200):
            
            raise ServerException(motor_lst)
            
        else:
            
            
            if(motor_lst==[]):
                
                response_dict = {}
                response_dict["motor_names"] = []
                response_dict["motor_desc"] = []
                response_dict["created_at"] = []
                response_dict["error"] = "none"
                

                response_json = json.dumps(response_dict, indent = 4)
                return response_json
            
            else:
            
                
                lst = []
                desc = []
                creat = []
                
                for motor_names in motor_lst:
                    
                    lst.append(motor_names["motor_name"])
                    desc.append(motor_names["desc"])
                    creat.append(motor_names["created_at"])
                
                            
                response_dict = {}
                response_dict["motor_names"] = lst
                response_dict["motor_desc"] = desc 
                response_dict["created_at"] = creat
                response_dict["error"] = "none"
                
    
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Get All Motors: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Get All Motors: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)


@app.route('/create_vehicle_version', methods=['POST'])
def create_vehicle_version():
    
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        
        created_at = datetime.now().strftime("%d/%m/%Y %I:%M %p")
        
        msg, status = db.create_vehicle_version(user_id, vehicle_name,created_at)
            
        if(status!=200):
            raise ServerException(msg)
        
        response_dict = {}
        response_dict["error"] = "none"
        
        return response_dict
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "vehicle version: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "vehicle version: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        
        
@app.route('/create_application_version', methods=['POST'])
def create_application_version():
    
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        
        user_id = data.get("user_id").strip()
        
        application_name = data.get("application_name").strip()
        
        created_at = datetime.now().strftime("%d/%m/%Y %I:%M %p")
        msg, status = db.create_application_version(user_id, application_name,created_at)
        
        msg, status = db.create_application_version(user_id, application_name)
            
        if(status!=200):
            raise ServerException(msg)
        
        response_dict = {}
        response_dict["error"] = "none"
        
        return response_dict
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "application version: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "application version: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)


@app.route('/create_motor_version', methods=['POST'])
def create_motor_version():
    
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        
        motor_name = data.get("motor_name").strip()
        
        created_at = datetime.now().strftime("%d/%m/%Y %I:%M %p")
        
        msg, status = db.create_motor_version(user_id, vehicle_name, motor_name,created_at)
            
        if(status!=200):
            raise ServerException(msg)
        
        response_dict = {}
        response_dict["error"] = "none"
        
        return response_dict
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "motor version: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "motor version: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        

@app.route('/delete_vehicle', methods=['POST'])
def delete_vehicle():
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        
        
        if(db.is_existing_vehicle(user_id, vehicle_name)):
            
            msg, status = db.delete_vehicle(user_id, vehicle_name)
            
            if(status!=200):
                raise ServerException(msg)
            else:
                
                response_dict = {}
                response_dict["error"] = "none"
            
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
                
       
        else:
            
            raise ServerException("Vehicle Doesn't exist")
        
            

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "delete vehicle: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "delete vehicle: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
    

@app.route('/delete_motor', methods=['POST'])
def delete_motor():
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        
        motor_name = data.get("motor_name").strip()
        
        if(db.is_existing_motor(user_id, vehicle_name, motor_name)):
            
            msg, status = db.delete_motor(user_id, vehicle_name, motor_name)
            
            if(status!=200):
                raise ServerException(msg)
            else:
                
                response_dict = {}
                response_dict["error"] = "none"
            
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
                
       
        else:
            
            raise ServerException("Motor Doesn't exist")
        
        
    
        
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "delete motor: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "delete motor: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)


@app.route('/delete_all_motors', methods=['POST'])
def delete_all_motor():

    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        
        
        if(db.is_existing_vehicle(user_id, vehicle_name) or db.is_existing_application(user_id, vehicle_name)):
            
            msg, status = db.delete_all_motors(user_id, vehicle_name)
            
            if(status != 200):
                
                raise ServerException(msg)
                
            else:
                
                response_dict = {}
                response_dict["error"] = "none"
            
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
            
        else:
            
            raise ServerException("Vehicle Does Not Exist")
        
    
        
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "delete all motors: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "delete all motors: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        


@app.route('/download', methods=['POST'])
def download():
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        vehicle_name = data.get("vehicle_name").strip()
        
        motor_name = data.get("motor_name").strip()
        
        
        dir_name = "data/{}/{}/{}".format(str(user_id), str(vehicle_name), str(motor_name))
        #shutil.make_archive(dir_name+"/download", 'zip', dir_name+"/download")

        #return send_from_directory(dir_name,"download.zip", as_attachment=True)
        
        file_content, status = db.retrieve_zip_file_from_downloads(user_id, vehicle_name, motor_name)
        
        if(status != 200):
            raise ServerException(file_content)
        
        file_stream = io.BytesIO(file_content)
        file_stream.seek(0)

        return send_from_directory(dir_name,"download.zip", as_attachment=True)

        
        
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Download : {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Download : {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)


@app.route('/', methods=['GET'])
def default_api():
    return "Welcome Motor Mojo!"
 

@app.route('/create_application', methods=['POST'])
def create_application():
    
    #os version no db
    try:
        
                
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        application_name = data.get("application_name").strip()
        desc = data.get("desc")
        
        if(desc==None or desc==""):
            desc = application_name
            
        if(db.is_existing_application(user_id, application_name)):
            
            raise ServerException("Application Name Already Exists")
        
        error = "none"

        created_at = datetime.now().strftime("%d/%m/%Y %I:%M %p")
        
        msg, status = db.create_application(user_id, application_name, [], [], error, desc,created_at)
        
        if(status != 200):
            
            raise ServerException(msg)
        
        response_dict = {}
        response_dict["error"] = "none"
        
        response_json = json.dumps(response_dict, indent = 4)
        return response_json
            
        
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "create Application: {}".format(str(e))
        response_dict["status"] = "error"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["status"] = "error"
        response_dict["error"] = "create Application: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
       
        
       
@app.route('/get_all_applications', methods=['POST'])
def get_applications():
    
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        
        msg, status = db.get_all_applications(user_id)
        
        if(status!=200):
            
            raise ServerException(msg)
            
        else:
            
            if(type(msg)==str):
                applications = json.loads(msg)
            else:
                applications=msg
            
            if(applications==[]):
                
                response_dict = {}
                response_dict["application_names"] = []
                response_dict["application_desc"] = [] 
                response_dict["created_at"] = []
                response_dict["error"] = "none"
                

                response_json = json.dumps(response_dict, indent = 4)
                return response_json
            
            else:
                
                application_names = []
                desc_lst = []
                create_lst = []
                
                for v in applications:
                    
                    application_names.append(v["application_name"])
                    desc_lst.append(v["desc"])
                    create_lst.append(v["created_at"])
                    
                response_dict = {}
                response_dict["application_names"] = application_names
                response_dict["application_desc"] = desc_lst
                response_dict["created_at"] = create_lst
                response_dict["error"] = "none"
                
    
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
        
        
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Get All applications: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Get All applications: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        

        
@app.route('/get_application_data', methods=['POST'])
def get_application_data():
    
    #os version no db
    
    try:
        
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        application_name = data.get("application_name").strip()
        
        
        application_data, status = db.get_application_data(user_id, application_name)
        
        if(status != 200):
            
            raise ServerException(application_data)
            
        else:
            
            input_data = application_data["input"]
            output_data = application_data["output"]
            
            
            if(input_data==[] or input_data==None or output_data==[] or output_data==None):
                
                #for new vehicle
                
                response_dict = {}
                response_dict["input"] = []
                response_dict["output"] = []
                response_dict["error"] = "none"
                
    
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
            
            else:
                
                response_dict = {}
                response_dict["input"] = input_data
                response_dict["output"] = output_data 
                response_dict["error"] = "none"
                
    
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
                
    
    
    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Get applications data: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Get applications data: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        

@app.route('/delete_application', methods=['POST'])
def delete_application():
    #os version no db
    try:
        
        data = request.json
    
        #email = data.get("email")
        user_id = data.get("user_id").strip()
        
        application_name = data.get("application_name").strip()
        
        
        if(db.is_existing_application(user_id, application_name)):
            
            msg, status = db.delete_application(user_id, application_name)
            
            if(status!=200):
                raise ServerException(msg)
            else:
                
                response_dict = {}
                response_dict["error"] = "none"
            
                response_json = json.dumps(response_dict, indent = 4)
                return response_json
                
       
        else:
            
            raise ServerException("application Doesn't exist")
        
            

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "delete application: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "delete application: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)

@app.route('/application_dynamics', methods=['POST'])
def application_dynamics():
    


    try:
        #POST all the variables for vehicle dynamics

        data = request.json

        #email = data.get("email")
        user_id = data.get("user_id").strip()
        application_name = data.get("application_name").strip()
        
        coords = data.get("points")
        
        duty = data.get("duty")
        
        v_dc = float(data.get("input_power"))
        v_type = str(data.get("power_type"))
        
        
        
        if(type(coords) == str):
            
            coords = json.loads(coords)
            
        
        if(type(duty) == str):
            
            coords = json.loads(duty)
            
        print(coords)
        print(duty)

        
        if(v_type=="2"):
            
            v_dc = np.sqrt(2)*v_dc
        elif(v_type == "3"):
            
            v_dc = np.sqrt(6)*v_dc
            
           
        tau_pk = 0
        max_no_load_rpm = 0
            
        for i in range(len(duty)):
            
            coord = coords[i]
            
            if(len(coord) == 2):
                
                tau = coord[1]
                max_no_load = coord[0]
                
            else:
                
                #size == 4
                
                x1,y1,x2,y2 = coord
                
                m = (y1-y2)/(x2-x1)
                
                delta = y2/m
                
                max_no_load = x2+delta
                
                tau = y1
                
            
            if(tau > tau_pk):
                tau_pk = tau
                
            if(max_no_load > max_no_load_rpm):
                max_no_load_rpm = max_no_load
                
        
            
        max_a_rms = (np.sqrt(2)*tau_pk*max_no_load_rpm)/(9.55*np.sqrt(3)*v_dc)
            
        
        

            
        input_data = data
        output_data = {
            "peak_torque": tau_pk,
            "max_no_load_rpm": max_no_load_rpm,
            "i_max": max_a_rms,
            "v_dc": v_dc
        }
        error = "none"
        desc = application_name
        
        if(not db.is_existing_application(user_id, application_name)):
            
            db.create_application(user_id, application_name, input_data, output_data, error, desc)
            
        else:
            
            db.update_application(user_id, application_name, input_data, output_data, error)

            response_json = json.dumps(output_data, indent = 4)
            return response_json

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Server Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Value Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
    
@app.route('/send_email_forgot_password', methods=['POST'])
def send_email_forgot_password():
    
    try:
        #POST all the variables for vehicle dynamics

        data = request.json

        #email = data.get("email")
        email = data.get("email").strip()
        
        if(not otp_db.email_exists(email)):
            raise ServerException("Email not found in user database")
        
        #generate otp
        otp = otp_db.generate_otp()
        
        #store otp in db
        msg, status = otp_db.store_otp(email, otp)
        
        if(status!=200):
            raise ServerException(msg)
            
            
        #send email
        
        msg, status = otp_email.send_otp_email(email, otp)
        
        if(status!=200):
            raise ServerException(msg)
        
        
        response_dict = {}
        response_dict["message"] = msg
        response_dict["error"] = "none"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
        
        
        

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Send Email Server Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Send Email Value Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json

@app.route('/verify_otp', methods=['POST'])
def verify_otp():
    
    try:
        #POST all the variables for vehicle dynamics

        data = request.json

        #email = data.get("email")
        email = data.get("email").strip()
        otp = int(data.get("otp").strip())
        
        if(not otp_db.email_exists(email)):
            raise ServerException("Email not found in user database")
        
       
        #verify otp
        msg, status = otp_db.verify_otp(email, otp)
        
        if(status!=200):
            raise ServerException(msg)
            
        
        response_dict = {}
        response_dict["message"] = msg
        response_dict["error"] = "none"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
        
        
        

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Verify OTP Server Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Verify OTP Value Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json


@app.route('/reset_password', methods=['POST'])
def reset_password():
    
    try:
        #POST all the variables for vehicle dynamics

        data = request.json

        #email = data.get("email")
        email = data.get("email").strip()
        password = data.get("password").strip()
        otp = int(data.get("otp").strip())
        
        if(not otp_db.email_exists(email)):
            raise ServerException("Email not found in user database")
        
       
        #re-verify otp to prevent unauthorized use of api
        msg, status = otp_db.verify_otp(email, otp)
        
        if(status!=200):
            raise ServerException(msg)
            
            
        #update password
        
        msg, status = otp_db.update_password(email, password)
            
        
        response_dict = {}
        response_dict["message"] = msg
        response_dict["error"] = "none"

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
        
        
        
        

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Password Reset Server Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json
    
    except ValueError as e:
        
        response_dict = {}
        response_dict["error"] = "Password Reset Value Error: {}".format(str(e))

        response_json = json.dumps(response_dict, indent = 4)
        return response_json


@app.route('/update_profile', methods=['POST'])
def update_profile_endpoint():
    try:
        data = request.json
        user_id = data.get("user_id").strip()
        new_profile_data = {
            "user_name": data.get("user_name"),
            "email": data.get("email"),
            "about": data.get("about"),
            "location": data.get("location")
        }

        msg, status = db.update_profile(user_id, new_profile_data)
        if status != 200:
            raise ServerException(msg)
        else:
            response_dict = {"message": "Profile updated successfully", "error": "none"}
            response_json = json.dumps(response_dict, indent=4)
            return response_json
    except ServerException as e:
        response_dict = {"error": "Update Profile: {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 500
    except ValueError as e:
        response_dict = {"error": "Update Profile: {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 500

@app.route('/get_profile', methods=['POST'])
def get_profile():
    try:
        data = request.json
        user_id = data.get("user_id").strip()

        profile, status = db.get_user_profile(user_id)

        if status != 200:
            raise ServerException(profile)
    
        response_dict = {}
        response_dict["profile"] = profile
        response_dict["error"] = "none"
        
        response_json = json.dumps(response_dict, indent=4)
        return response_json

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Get Profile: {}".format(str(e))

        response_json = json.dumps(response_dict, indent=4)
        return response_json

    except ValueError as e:
        response_dict = {}
        response_dict["error"] = "Get Profile: {}".format(str(e))

        response_json = json.dumps(response_dict, indent=4)
        return response_json


@app.route('/change_password', methods=['POST'])
def change_password():
    try:
        data = request.json
        user_id = data.get("user_id").strip()
        current_password_hash = data.get("current_password").strip()
        new_password_hash = data.get("new_password").strip()

        result, status = db.change_user_password(user_id, current_password_hash, new_password_hash)
        if status != 200:
            raise ServerException(result)

        response_dict = {"error": "none"}
        response_json = json.dumps(response_dict, indent=4)
        return response_json
    
    except ServerException as e:
        response_dict = {"error": " {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 400

    except ValueError as e:
        response_dict = {"error": " {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 400

    except Exception as e:
        response_dict = {"error": " {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 500

@app.route('/get_userdata', methods=['POST'])
def get_userdata():
    try:
        data = request.json
        user_id = data.get("user_id").strip()

        profile, status = db.get_user_data(user_id)

        if status != 200:
            raise ServerException(profile)
    
        response_dict = {}
        response_dict["user_data"] = profile
        response_dict["error"] = "none"
        
        response_json = json.dumps(response_dict, indent=4)
        return response_json

    except ServerException as e:
        response_dict = {}
        response_dict["error"] = "Get Data: {}".format(str(e))

        response_json = json.dumps(response_dict, indent=4)
        return response_json

    except ValueError as e:
        response_dict = {}
        response_dict["error"] = "Get Data: {}".format(str(e))

        response_json = json.dumps(response_dict, indent=4)
        return response_json


@app.route('/update_feedback', methods=['POST'])
def update_feedback():
    try:
        data = request.json
        user_id = data.get("user_id").strip()
        feedback_data = {
            "feedback": data.get("feedback"),
            "rating": data.get("rating"),
        }

        msg, status = db.update_profile(user_id, feedback_data)
        if status != 200:
            raise ServerException(msg)
        else:
            response_dict = {"message": "Thanks for Rating", "error": "none"}
            response_json = json.dumps(response_dict, indent=4)
            return response_json
    except ServerException as e:
        response_dict = {"error": "Update Rating: {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 500
    except ValueError as e:
        response_dict = {"error": "Update Rating: {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 500


@app.route('/update_all_feedback', methods=['POST'])
def update_all_feedback():
    try:
        data = request.json
        user_id = data.get("user_id").strip()
        feedback_data = {
            "feedback": data.get("feedback"),
            "timestamp": data.get("timestamp"),
        }

        msg, status = db.store_feedback(user_id, feedback_data)
        if status != 200:
            raise ServerException(msg)
        else:
            response_dict = {"message": "Thanks for Your Feedback", "error": "none"}
            response_json = json.dumps(response_dict, indent=4)
            return response_json
    except ValueError as e:
        response_dict = {"error": "Update feedback: {}".format(str(e))}
        response_json = json.dumps(response_dict, indent=4)
        return response_json, 500



if __name__ == '__main__':
    try:
        #app.run()
        #app.run(host='0.0.0.0', port=80, ssl_context=("cert.pem","key.pem"))
        #default port = 5000
        serve(app, host='0.0.0.0', port = 80, threads=2)
        # serve(app, listen="127.0.0.1:5000", threads=2)
    except Exception as e:

        print("Server:",str(e))
