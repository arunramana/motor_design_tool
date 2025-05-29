import numpy as np
from numpy import pi
from .trig import *


def convertMechAngleToElecAngle(mechAngle,polePairs):
    
    mechAnglesToOneElectricalRotation = 360.0/polePairs
    mechAngleMod = mechAngle%mechAnglesToOneElectricalRotation
    electricalAngle = mechAngleMod*360.0/mechAnglesToOneElectricalRotation
    return electricalAngle



def inv_park(i_d,i_q,rot_angle,polePairs):
    
    #returns i alpha and i beta
    #rot_angle is mechanical angle 
    #convert it to electrical

    rot_angle = convertMechAngleToElecAngle(rot_angle,polePairs)
    
    ialpha = i_d*cos(rot_angle) - i_q*sin(rot_angle)
    
    ibeta = i_q*cos(rot_angle) + i_d*sin(rot_angle)
    
    return ialpha,ibeta


    
def inv_clarke(ialpha,ibeta):
    
    ia = ialpha
    
    ib = (-ialpha + ((3**0.5)*ibeta))/2
    
    ic = (-ialpha - ((3**0.5)*ibeta))/2
    
    return ia,ib,ic

def clarke(wa,wb,wc):
    
    #returns walpha,wbeta
    
    walpha = wa
    wbeta = (wa+2*wb)/(3**0.5)
    
    
    return walpha,wbeta

def park(walpha,wbeta,rot_angle,polePairs):
    
    #return wd and wq
    #rot_angle is mechanical angle 
    #convert it to electrical

    rot_angle = convertMechAngleToElecAngle(rot_angle,polePairs)
    
    wd = walpha*cos(rot_angle) + wbeta*sin(rot_angle)
    
    wq = wbeta*cos(rot_angle) - walpha*sin(rot_angle)
    
    return wd,wq