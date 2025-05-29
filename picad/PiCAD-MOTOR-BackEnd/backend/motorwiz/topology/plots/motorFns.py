import numpy as np
import scipy.optimize as optimize

def getElecFreq(RPM,motorParams):
    return RPM/60.0 * 2*np.pi*motorParams["polePairs"]


def getShortCircuitCurrentPt(motorParams):
    centrePoint_x = motorParams["psiM_Vs_RMS"]/motorParams["Ld"]
    centrePoint_y = 0
    return (centrePoint_x,centrePoint_y)

def getVoltageEllipsePts(RPM,motorParams,controlParams,IdRange):   
    e_centrePt_x,e_centrePt_y = getShortCircuitCurrentPt(motorParams)
    op= {}
    op["Vs_Ph_max"] = controlParams["DCVoltage"]/np.sqrt(2)/np.sqrt(3)
    op["we"] = RPM/60*2*np.pi *motorParams["polePairs"]
    op["Xd"] = motorParams["Ld"] * op["we"]
    op["Xq"] = motorParams["Lq"] * op["we"]
    k2 = (op["Vs_Ph_max"]/op["Xd"])**2
    k3 = (op["Vs_Ph_max"]/op["Xq"])**2
    ellipse_pts = []
    pt1,pt2 = IdRange
    for Id in range(pt1,pt2,-5):
        k1 = (Id + e_centrePt_x)**2
        ellipsePt_x = Id
        ellipsePt_y1 = np.sqrt((1-(k1/k2)) * k3)
        ellipse_pts.append((ellipsePt_x,ellipsePt_y1))
    return ellipse_pts

def getConstantTorqueCurves(torque,motorParams,calcParams,IdRange):
    kt = calcParams["Kt"]
    Iq = torque/kt
    
    hyperBolaPts = []
    pt1,pt2 = IdRange
    for Id in range(pt1,pt2,-10):
        b = optimize.Bounds(0,np.inf)
        cp_OUT = optimize.minimize(solveForIq_fromTorque,Iq,(torque,Id,motorParams),bounds = b)
        if cp_OUT.success:
            Iq = round(cp_OUT.x[0],2)
            hyperBolaPts.append((Id,Iq))
        
    return hyperBolaPts

def getGammaTMax(currentRMS,motorParams):
    #forumula comes from the formula for torque, which we differentiate wrt to gamma and find the max
    #its in the 
    Ld = motorParams["Ld"]
    Lq = motorParams["Lq"]
    psiM = motorParams["psiM_Vs_RMS"]
    Xd = Ld * getElecFreq(1000,motorParams)
    Xq = Lq * getElecFreq(1000,motorParams)
    Eq = psiM * getElecFreq(1000,motorParams)
    #get gammaT max
    deltaV = (Xd - Xq) * currentRMS
    gammaT_max =  np.arcsin(0.25*((Eq/deltaV) + np.sqrt((Eq/deltaV)**2 + 8)))
    return gammaT_max

def getSecondQuadrantCircle(radius):
    circle = np.zeros((90,3))
    circle[:,0] = np.arange(0,90,1)
    circle[:,1] = radius * np.sin(np.deg2rad(circle[:,0]))
    circle[:,2] = radius * np.cos(np.deg2rad(circle[:,0]))   
    return circle[:,1:3]

def solveForCornerPoint(RPM,Id,Iq,motorParams,maxDC):
    w = getElecFreq(RPM,motorParams)
    Xd = motorParams["Ld"]*w
    Xq = motorParams["Lq"]*w
    E = motorParams["psiM_Vs_RMS"] * w
    maxV_LN = maxDC/np.sqrt(2)/np.sqrt(3)
    #without resistance
    Vd = -Xq*Iq
    Vq = Xd*Id + E
    V = np.sqrt(Vd**2 + Vq**2)
    return abs(V - maxV_LN)
  
def solveForId(Id,RPM,motorParams,maxDC,maxRMS):
    w = getElecFreq(RPM,motorParams)
    Xd = motorParams["Ld"]*w
    Xq = motorParams["Lq"]*w
    E = motorParams["psiM_Vs_RMS"] * w
    maxV_LN = maxDC/np.sqrt(2)/np.sqrt(3)
    Iq = np.sqrt(maxRMS**2 - Id**2)
    #without resistance
    Vd = -Xq*Iq
    Vq = Xd*Id + E
    V = np.sqrt(Vd**2 + Vq**2)
    return abs(V - maxV_LN)

def solveForIq_fromTorque(Iq,torque_in,Id,motorParams):
    w = getElecFreq(1000,motorParams) #rpm doesnt matter
    Xd = motorParams["Ld"]*w
    Xq = motorParams["Lq"]*w
    E = motorParams["psiM_Vs_RMS"] * w
    torque = 3*motorParams["polePairs"]*(Iq*E + Id*Iq*(Xd-Xq))/w
    return abs(torque_in - torque)


def solveForTorque(Id,Iq,RPM,motorParams):
    w = getElecFreq(RPM,motorParams)
    Xd = motorParams["Ld"]*w
    Xq = motorParams["Lq"]*w
    E = motorParams["psiM_Vs_RMS"] * w
    torque = 3*motorParams["polePairs"]*(Iq*E + Id*Iq*(Xd-Xq))/w
    return torque

def getTorque(Id,Iq,motorParams):
    psiM = motorParams["psiM_Vs_RMS"]
    Ld = motorParams["Ld"]
    Lq = motorParams["Lq"]
    p = motorParams["polePairs"]
    torque = 3*p*(Iq*psiM + Id*Iq*(Ld-Lq))
    return torque


#--------for later----------------------#
# def getLdLq(current,currentAngle):
#     return 22*1e-6,36*1e-6

# def getResistance(temperature):
#     return 4.2*1e-3

# def getPsiM(temperature):
    
