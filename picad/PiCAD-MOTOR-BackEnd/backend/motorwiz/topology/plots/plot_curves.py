from ...trig import *

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.patches import Ellipse

import pandas as pd
import os
import shutil

import warnings
warnings.filterwarnings("ignore")

#important
# all inputs are rms,max voltage is vdc/sqrt(2), 

def curves_to_dict(currentPk,max_voltage,ld,lq,wd0,pole_pairs, project_name, magnet_position, max_no_load_rpm):
    
    wd_rms = wd0/np.sqrt(2)
    
    #change
    wd0 = wd_rms
    
    
    MTPA_dct = getMTPA(currentPk,ld,lq,wd0,pole_pairs,project_name, magnet_position)
    
    id_max = MTPA_dct["id"][-1]
    iq_max = MTPA_dct["iq"][-1]
    torquePk = MTPA_dct["torque"][-1]
    
    
    #print("--------------------------------------")
    #print("\nplot curves: id: {}, iq: {}, i {}, currentPK: {}, tau: {}".format(id_max, iq_max, (id_max**2+iq_max**2)**0.5, currentPk, torquePk))
    
    
    torque_speed_dct = getTorqueSpeedCurve(id_max,iq_max,currentPk,max_voltage,ld,lq,wd0,pole_pairs,max_no_load_rpm,project_name=project_name)
    
    #max_no_load_rpm = torque_speed_dct["operating torque speed"][0][-1]
    
    currentPk = (id_max**2+iq_max**2)**0.5
    current_circle_dct = getCurrentCircle(currentPk)
    
    ellipse_dct = getVoltageEllipse(max_voltage,ld,lq,wd0,currentPk,max_no_load_rpm,pole_pairs)
    
    hyp_dct = getConstantTorqueHyperbola(torquePk,currentPk,ld,lq,wd0,pole_pairs)
    
    curves_dct = {}
    
    curves_dct["PEAK TORQUE"] = torquePk
    curves_dct["MAX NL RPM"] = max_no_load_rpm
    curves_dct["MTPA"] = MTPA_dct
    curves_dct["TORQUE SPEED"] = torque_speed_dct
    curves_dct["CURRENT CIRCLE"] = current_circle_dct
    curves_dct["VOLTAGE ELLIPSE"] = ellipse_dct
    curves_dct["CONSTANT TORQUE"] = hyp_dct
    
    return curves_dct
     


def getTorqueSpeedCurve(id_max,iq_max,currentPk,voltage,ld,lq,wd0,pole_pairs,max_no_load_rpm,project_name=None):
    
    #get id_max and iq_max from MTPA
    torque_speed_dct = {}
    
    id = id_max
    iq = iq_max
    vs = voltage/np.sqrt(3)#phase RMS
    
    current = currentPk
    
    
    t = get_torque_wRMS(id,iq,ld,lq,wd0,pole_pairs)
    
    omega = get_omega_from_voltage(id,iq,ld,lq,wd0,pole_pairs,vs)
    
    
    t_arr = [t]
    operating_t = [t]
    operating_o = [0]
    o_arr = [0]
    
    t_arr.append(t)
    o_arr.append(omega)
    operating_t.append(t)
    operating_o.append(omega)
    
    #for testing purpose
    
    test_dct = {
        "t": [],
        "s": [],
        "id": [],
        "iq": [],
    }

    
        
        
    for i_d in range(int(abs(id))+1,int(current)+1,1):
    
        i_d = -i_d
        
        iq = np.sqrt(current**2-i_d**2)
    
        if(iq!=0):
            theta = taninv(i_d/iq)
        else:
            theta = 90
    
    
        
        t = get_torque_wRMS(i_d,iq,ld,lq,wd0,pole_pairs)
        omega = get_omega_from_voltage(i_d,iq,ld,lq,wd0,pole_pairs,vs)
        
    
    
        if(abs(theta) <= 80):
            #to close curve as magnet gets demagnetized
            operating_t.append(abs(t))
            operating_o.append(omega)
                        
            test_dct["t"].append(abs(t))
            test_dct["s"].append(omega)
            test_dct["id"].append(i_d)
            test_dct["iq"].append(iq)
            
        if(omega > 2*max_no_load_rpm):
            break
    
    
        t_arr.append(abs(t))
        o_arr.append(omega)
        
    
    
    #close the curve for operating points
    operating_t.append(0)
    operating_o.append(operating_o[-1])
    
    torque_speed_dct["torque speed"] = [o_arr,t_arr]
    torque_speed_dct["operating torque speed"] = [operating_o, operating_t]
    
    '''
    if(project_name!=None):
        pd.DataFrame(test_dct).to_csv("{}/download/test_ts.csv".format(project_name))
    '''
    
    return torque_speed_dct
    
    

def getMTPA(currentPk,ld,lq,wd0,pole_pairs,project_name, magnet_position):

    #MTPA

    current = currentPk#/np.sqrt(2)
    
    currents = np.arange(1,int(current)+10,10)
    
    id_arr = [] #x
    iq_arr = [] #y
    tau_arr = []
    current_angles = []
    
    MTPA_dct = {}

 
    
    for i in currents:
    
        tau_prev = None
        id_max = 0
        iq_max = 0
        current_angle = 0
        
        for cangle in range(0,181):
            
            ca = cangle*0.5
    
            id = -i*sin(ca)
            iq = i*cos(ca)
    
            t = get_torque_wRMS(id,iq,ld,lq,wd0,pole_pairs)
            
    
            if(tau_prev==None):
                tau_prev = t
                id_max = id
                iq_max = iq
                current_angle = -ca
    
            else:
    
                if(abs(tau_prev) < abs(t)):
                    tau_prev = t
                    id_max = id
                    iq_max = iq
                    current_angle = -ca
    
    
        id_arr.append(id_max)
        iq_arr.append(iq_max)
        tau_arr.append(tau_prev)
        current_angles.append(current_angle)
    
    
    MTPA_dct["id"] = id_arr
    MTPA_dct["iq"] = iq_arr
    MTPA_dct["torque"] = tau_arr
    MTPA_dct["current angle"] = current_angles
    
    dr = project_name+"/download"
    
    if(not os.path.isdir(dr)):
        os.mkdir(dr)
        
        
    if(magnet_position=="SURFACE"):
        
        if(os.path.isfile(dr+"/MTPA.csv")):
            os.remove(dr+"/MTPA.csv")
            
    else:
        
        pd.DataFrame(MTPA_dct).to_csv("{}/MTPA.csv".format(dr), index=False)
    
    return MTPA_dct


def getCurrentCircle(currentPk):

    #for II Quadrant 
    
    current = currentPk#/np.sqrt(2)
    
    x = -np.arange(0,current,1)
    y = np.sqrt(current**2-x**2)
    
    x = x.tolist()
    y = y.tolist()
    
    x.append(-current)
    y.append(0)
    
    circle_dict = {str(round(current)): [x,y]}
    
    return circle_dict


'''
def getVoltageEllipse(max_voltage,ld,lq,psi_m,max_no_load_rpm,pole_pairs,no_ellipse=5):
    
    #get max no load rpm from MTPA
    wd0 = abs(psi_m)
    
    vs = max_voltage/np.sqrt(2)/np.sqrt(3)

    #ellipse
    center = (-wd0/ld,0)
    
    steps = int(max_no_load_rpm/no_ellipse)
    
    #format -> ellipse_dct = {"omega": [x_arr,y_arr]}
    ellipse_dict = {}
    
    for omega in range(steps,int(max_no_load_rpm),steps):
    
        we = 2*np.pi*pole_pairs*omega/60
        #we = omega
        a = vs/(lq*we)
        b = vs/(ld*we)
    
        voltage_ellipse = Ellipse(xy=center, width=b, height=a, angle=0)
        
        
        path = voltage_ellipse.get_path()
        # Get the list of path vertices
        vertices = path.vertices.copy()
        # Transform the vertices so that they have the correct coordinates
        vertices = voltage_ellipse.get_patch_transform().transform(vertices)
    
        x,y = vertices.T
        
        ellipse_dict[str(omega)] = [x.tolist(),y.tolist()]
        
    
    
    return ellipse_dict
'''    

'''
def getVoltageEllipse(max_voltage,ld,lq,psi_m,max_no_load_rpm,pole_pairs,no_ellipse=6):
    
    #get max no load rpm from MTPA
    wd0 = abs(psi_m) #rms
    
    vs = max_voltage/np.sqrt(3)#vrms phase

    #ellipse
    center = (-wd0/ld,0)
    
    max_no_load_rpm = 5000
    
    steps = int(max_no_load_rpm/no_ellipse)
    
    #format -> ellipse_dct = {"omega": [x_arr,y_arr]}
    ellipse_dict = {}
    
    for omega in range(steps,int(max_no_load_rpm),steps):
    
        we = 2*np.pi*pole_pairs*omega/60
        #we = omega
        a = vs/(lq*we)
        b = vs/(ld*we)
        h = center[0]
        
        step_a = a/100
        step_b = (b+h)/100
        
        x_arr = []
        y_arr = []
        
        #x=0, y = a*(1-h^2/b^2)^0.5
        
        x_arr.append(0)
        
        y = a*np.sqrt(abs(1-(h/b)**2))
        y_arr.append(y)
        
        x=step_b
    
        
        #major axis || x axis
        while x <= b+h:
            
            #equation of ellipse x = -x for II Quad
            y = 1-((-x-h)**2)/(b**2)
            
            y = a*np.sqrt(abs(y))
            
            if(x>b or y>a):
                continue
            
            x_arr.append(-x)
            y_arr.append(y)
                
            x+=step_b
            
        
        #x=-b-h, y=0
        
        x_arr.append(-(abs(b)+abs(h)))
        y_arr.append(0)
        
        ellipse_dict[str(omega)] = [x_arr,y_arr]
        
    return ellipse_dict
        
'''
def getVoltageEllipse(max_voltage,ld,lq,psi_m,i_max,max_no_load_rpm,pole_pairs,no_ellipse=6):
    
    
    steps = int(max_no_load_rpm/no_ellipse)
    
    #format -> ellipse_dct = {"omega": [x_arr,y_arr]}
    ellipse_dict = {}
    
    for omega in range(steps,int(max_no_load_rpm),steps):
        
        x_arr,y_arr = getVoltageEllipsePoint(max_voltage,ld,lq,psi_m,i_max,omega,pole_pairs)
        ellipse_dict[str(omega)] = [x_arr,y_arr]
        
        
    return ellipse_dict 

def getVoltageEllipsePoint(max_voltage,ld,lq,psi_m,i_max,omega,pole_pairs):   
    
    wd0 = abs(psi_m)
    e_centrePt_x,e_centrePt_y = (wd0/ld,0)
    op= {}
    vs_ph = max_voltage/np.sqrt(3)
    
    we = 2*np.pi*pole_pairs*omega/60
    xd = ld*we
    xq = lq*we 
    
    k2 = (vs_ph/xd)**2
    k3 = (vs_ph/xq)**2
    ellipse_pts = []
    pt1,pt2 = (0,-int(i_max))
    
    b = vs_ph/(ld*we)
    
    x,y = [],[]
    
    for Id in range(pt1,pt2,-5):
        k1 = (Id + e_centrePt_x)**2
        ellipsePt_x = Id
        ellipsePt_y1 = (1-(k1/k2)) * k3
        
        if(ellipsePt_y1<0):
            continue
        
        ellipsePt_y1 = np.sqrt(ellipsePt_y1)
        
        if(np.isnan(ellipsePt_x) or np.isnan(ellipsePt_y1)):
            continue
        
        x.append(ellipsePt_x)
        y.append(ellipsePt_y1)
        
    x.append(-(b+e_centrePt_x))
    y.append(0)
        
    return x,y 

def getConstantTorqueHyperbola(torquePk,currentPk,ld,lq,wd0,pole_pairs,no_hyp=3):
    
    #hyperbola
    #get torque peak from MTPA

    i_max = currentPk
    torque_steps =  int(torquePk/no_hyp)
    current_steps = int(i_max/no_hyp)
    
    constant_torques_dct = {}
    

    for torque in range(torque_steps, int(torquePk)+1,torque_steps):
    
        id_arr = -np.array(list(range(current_steps,int(i_max)*2,current_steps)))
    
        iq_arr = get_iq_from_id_rms(id_arr,torque,pole_pairs,wd0,ld,lq)
    
        #for y intercept
        iq = get_iq_from_id_rms(0,torque,pole_pairs,wd0,ld,lq)
        
        id_arr = np.append([0],id_arr)
        iq_arr = np.append([iq],iq_arr)
        
        constant_torques_dct[str(torque)] = [id_arr.tolist(),iq_arr.tolist()]
        

    return constant_torques_dct
    
        
    

def get_iq_from_id_rms(id,torque,pole_pairs,wd0,ld,lq):

    #get iq given id and torque
    A = (1/3)*(torque/pole_pairs)
    
    return A*(1/(abs(wd0)+(ld-lq)*id))     



def get_iq_from_id(id,torque,pole_pairs,wd0,ld,lq):

    #get iq given id and torque
    A = (2/3)*(torque/pole_pairs)
    
    return A*(1/(abs(wd0)+(ld-lq)*id))


def get_torque(id,iq,ld,lq,wd0,pole_pairs):

    tau = (3/2)*pole_pairs*(iq*abs(wd0)+(ld-lq)*id*iq)
    return tau

def get_torque_wRMS(id,iq,ld,lq,wd0,pole_pairs):

    tau = 3*pole_pairs*(iq*abs(wd0)+(ld-lq)*id*iq)
    return tau

def get_omega_from_voltage(id,iq,ld,lq,wd0,pole_pairs,voltage):
    #returns in rpm
    
    vs = voltage#rms phase voltage
    
    A = 1/((ld*id+abs(wd0))**2+(lq*iq)**2)

    omega = (vs/pole_pairs)*np.sqrt(A)

    return omega*9.55


