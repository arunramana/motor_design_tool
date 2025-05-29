import numpy as np
import matplotlib.pyplot as plt
import backend.motorwiz.topology.plots.motorFns as mtr
from pprint import pprint as pprint
import scipy.optimize as optimize
import warnings


def getMotorOperatingRegions(motorParams,controlParams,calcParams,plot_results = 0):
    #return the curves sepreately so that Arun can draw it
    #first get the current Max circle
    maxCurrent = controlParams["RMS_max"]  
    cMax_circle = mtr.getSecondQuadrantCircle(maxCurrent)
    cMax_circle[:,0] = -cMax_circle[:,0]
    #then get the voltage speed ellipses
    #calculate the RPM range by one
    maxNL_rpm = calcParams["NL_RPM_K"] *1000
    RPMS = np.arange(1000,int(maxNL_rpm*1.5),500)
    eArr_out = []
    pts_Id_range = [0,-maxCurrent - 100]
    pts_Iq_max = [maxCurrent + 100]
    for idx in range(RPMS.shape[0]):
        RPM = RPMS[idx]
        ellipsePts = mtr.getVoltageEllipsePts(RPM,motorParams,controlParams,pts_Id_range)
        eArr =np.array(ellipsePts)
        eArr_removedNans = eArr[np.argwhere(~np.isnan(eArr[:,1])).ravel()]
        if eArr_removedNans.shape[0] != 0:
            #again filter out values where the Iq is much greater than the max Current
            mask = eArr_removedNans[:,1] <= pts_Iq_max
            eArr_rangeFiltered = eArr_removedNans[mask]
            if eArr_rangeFiltered.shape[0] != 0:
                eArr_out.append((RPM,eArr_rangeFiltered))
        
    #then get constant torque curves
    dT = calcParams["maxTorque"]*0.2
    torques = np.arange(0.3*calcParams["maxTorque"],calcParams["maxTorque"]+dT,dT)
    tArr_out = []
    pts_Id_range = [0,-maxCurrent - 100]
    pts_Iq_max = [maxCurrent + 100]
    for idx in range(torques.shape[0]):
        torque_in = torques[idx]
        print ("--getting constant torque curve for torque = {}".format(torque_in))
        pts = mtr.getConstantTorqueCurves(torque_in,motorParams,calcParams,pts_Id_range)
        t_Arr =np.array(pts)
        t_Arr_removedNans = t_Arr[np.argwhere(~np.isnan(t_Arr[:,1])).ravel()]
        if t_Arr_removedNans.shape[0] != 0:
            #again filter out values where the Iq is much greater than the max Current
            mask = t_Arr_removedNans[:,1] <= pts_Iq_max
            t_Aarr_rangeFiltered = t_Arr_removedNans[mask]
            if t_Aarr_rangeFiltered.shape[0] != 0:
                tArr_out.append((torque_in,t_Aarr_rangeFiltered))
    
    if plot_results:
        f = plt.figure()
        plt.plot(cMax_circle[:,0],cMax_circle[:,1],'r',label='maxCurrent Circle') #current max circle
        for idx in range(len(eArr_out)):
            RPM,eArr_o = eArr_out[idx]
            plt.plot(eArr_o[:,0],eArr_o[:,1],'-')
            plt.text(eArr_o[-1][0],eArr_o[-1][1],str(int(RPM)),rotation=-30,rotation_mode ='anchor')
        for idx in range(len(tArr_out)):
            torque,tArr_o = tArr_out[idx]
            plt.plot(tArr_o[:,0],tArr_o[:,1],'--x')
            text_idx = 2 if tArr_o.shape[0] >= 2 else 0
            plt.text(tArr_o[text_idx][0],tArr_o[text_idx][1],str(round(torque,1))+"Nm",rotation=10,rotation_mode ='anchor')

        plt.plot([0,-150],[0,0],'k')
        plt.plot([0,0],[0,150],'k')   
        
        ax = plt.gca()
        ax.set_aspect('equal', 'box')
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel("Iq Arms ")
        plt.xlabel("Id Arms ")

        plt.legend()
        plt.title("Motor Operating Regions")
        
        plt.grid()
        plt.show()
        
    #return out data the way Arun wants it
    cMax = {}
    cMax["x"]=cMax_circle[:,0]
    cMax["y"]=cMax_circle[:,1]
    
    ellipses = {}
    for e in eArr_out:
        RPM,eA = e
        ellipses[str(RPM)]= {}
        ellipses[str(RPM)]["x"] = eA[:,0]
        ellipses[str(RPM)]["y"] = eA[:,1]
            
    constantTorques = {}
    for e in tArr_out:
        torque,tA = e
        torque = round(torque,2)
        constantTorques[str(torque)]= {}
        constantTorques[str(torque)]["x"] = tA[:,0]
        constantTorques[str(torque)]["y"] = tA[:,1]
                 
            
    return cMax,ellipses,constantTorques


def getMTPA(motorParams,controlParams,calcParams,plot_results = 0):
    maxCurrent = controlParams["RMS_max"]  
    currents = np.arange(0,maxCurrent,10)
    currents[-1]  = maxCurrent
    out = np.zeros((currents.shape[0],5))
    for idx in range(currents.shape[0]):
        current = currents[idx]
        gammaMax = mtr.getGammaTMax(current,motorParams)
        Id = -current*np.sin(gammaMax)
        Iq = current*np.cos(gammaMax)
        torque  = round(mtr.getTorque(Id,Iq,motorParams),2)
        out[idx] = (current,gammaMax,Id,Iq,torque)
    
    cMax_circle = mtr.getSecondQuadrantCircle(maxCurrent)
    cMax_circle[:,0] = -cMax_circle[:,0]
    if plot_results:
        f = plt.figure()
        plt.plot(cMax_circle[:,0],cMax_circle[:,1],'r',label='maxCurrent Circle') #current max circle
        plt.plot(out[:,2],out[:,3],'b--',label='MTPA curve')
        
        ax = plt.gca()
        ax.set_aspect('equal', 'box')
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()

        plt.legend()
        plt.grid()
        plt.show()
    
    #return out data the way Arun wants it
    cMax = {}
    cMax["x"]=cMax_circle[:,0]
    cMax["y"]=cMax_circle[:,1]
    
    mtpa = {}
    for idx in range(out.shape[0]):
        current,gammaMax,Id,Iq,torque = out[idx]
        mtpa[str(current)] = {"gammaMax" : gammaMax, "Id" : Id , "Iq": Iq , "torque" : torque}
    
    return mtpa,[cMax_circle,out]


def getTorqueSpeedCurve(motorParams,controlParams,calcParams,mtpa_arr,plot_results = 0):
    maxCurrent =  controlParams["RMS_max"]

    cMax_circle = mtr.getSecondQuadrantCircle(maxCurrent)
    cMax_circle[:,0] = -cMax_circle[:,0]
    
    maxCurrentInfo = mtpa_arr[mtpa_arr[:,0] == maxCurrent]
   
    gammaDeg = np.rad2deg(maxCurrentInfo[0][1])
    Id = maxCurrentInfo[0][2]
    Iq = maxCurrentInfo[0][3]
    b = optimize.Bounds(0,np.inf)
    initRPMguess = int(calcParams["NL_RPM_K"]/2*1000)
    
    cp_OUT = optimize.minimize(mtr.solveForCornerPoint,initRPMguess,(Id,Iq,motorParams,controlParams["DCVoltage"]),method='Nelder-Mead',bounds = b)
    
    if cp_OUT.success:
        cornerPointRPM = cp_OUT.x[0]
        print ("cornerPointRPM = {}".format(cornerPointRPM))
        t0 = mtr.solveForTorque(Id,Iq,cornerPointRPM,motorParams)
        #now find the rest of the points after the corner point
        RPMs = np.arange(cornerPointRPM+50,calcParams["NL_RPM_K"]*1000*2,100)
        out = np.zeros((RPMs.shape[0]+2,5))
        out[0] = (0,Id,Iq,gammaDeg,t0)
        out[1] = (cornerPointRPM,Id,Iq,gammaDeg,t0)
        fail = 0
        for idx in range(len(RPMs)):
            RPM = round(RPMs[idx],2)
            b = optimize.Bounds(-np.inf,0)
            # tol1 = 0.01
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cp_OUT = optimize.minimize(mtr.solveForId,Id,(RPM,motorParams,controlParams["DCVoltage"],maxCurrent),method='Nelder-Mead',bounds = b)
                if cp_OUT.success:
                    # print (cp_OUT)
                    Id = cp_OUT.x[0]
                    Iq = np.sqrt(maxCurrent**2 - Id**2)
                    theta = np.rad2deg(np.arctan(abs(Id)/Iq))
                    torque = mtr.solveForTorque(Id,Iq,RPM,motorParams)
                    out[idx+2] = (RPM,Id,Iq,theta,torque)
                    if theta > 89.9:
                        break
                else:
                    #print (cp_OUT)
                    fail = fail + 1
        print ("fail = " + str(fail))
        #remove all zeros rows in out 
        mask = np.any(out,axis = 1)
        out1 = out[mask]
        if plot_results:
            plt.figure()
            plt.plot(cMax_circle[:,0],cMax_circle[:,1],'r',label='maxCurrent Circle') #current max circle
            plt.plot(mtpa_arr[:,2],mtpa_arr[:,3],'b--',label='MTPA curve')
            plt.plot(out1[:,1],out1[:,2],'gx',label='constant voltage mode points')
            
            ax = plt.gca()
            ax.set_aspect('equal', 'box')
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()
            plt.legend()
            plt.grid()
            
            plt.figure()
            plt.plot(out1[:,0],out1[:,4],'b')
            plt.plot(out1[1,0],out1[1,4],'gx') # corner point is second point in array
            plt.plot()
            plt.title("Torque Speed Curve")
            plt.grid()
            
            plt.show()
            
        #return out data the way arun wants it
        #data for torqueSpeed curve
        torqueSpeedCurve = {}
        for idx in range(len(out1)):
            d = out1[idx]
            (RPM,Id,Iq,gammaDeg,t0) = d
            torqueSpeedCurve[str(RPM)] = {}
            torqueSpeedCurve[str(RPM)]["Id"] = Id
            torqueSpeedCurve[str(RPM)]["Iq"] = Iq
            torqueSpeedCurve[str(RPM)]["gammaDeg"] = gammaDeg
            torqueSpeedCurve[str(RPM)]["torque"] = t0
        return torqueSpeedCurve,out1
    else:
        print ("failed finding Corner Point!")
        return -1


def checkWithinValidRegion(x,y,tsArr):
    RPM = tsArr[:,0]
    torques = tsArr[:,4]

    #first see if torque is below max Torque
    maxTorque = np.max(torques)
    if y <= maxTorque:
        #find nearest torque idxs
        d1 = torques - y
        idx1 = np.argmin(abs(d1))
        if idx1 < len(torques)-1:
            idx2 = idx1 + 1
        else:
            idx2= idx1
        #get neighbouring points in the torque speed curve
        t1 = torques[idx1]
        t2 = torques[idx2]
        r1 = RPM[idx1]
        r2 = RPM[idx2]
        #interpolate
        if idx1 == idx2:
            ratio = 1
            rMax = r1
        else:
            ratio = (y - t1)/(t2-t1)
            rMax = r1 + ratio*(r2-r1)
        if x <= rMax:
            return 1
        else:
            return 0
    else:
        return 0

    
    
    
    
    
    
    