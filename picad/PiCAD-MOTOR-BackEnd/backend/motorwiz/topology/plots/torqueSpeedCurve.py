import numpy as np
import motorOperatingParamsController as mo
import motorFns as mtr
from pprint import pprint as pprint

'''
motorParams={"psiM_Vs_RMS":0.0116/np.sqrt(2),
              "polePairs":5,
              "Rs":4.2*1e-3,
              "Ld":14*1e-6,
              "Lq":27*1e-6,
             }
controlParams={"DCVoltage":48,
               "RMS_max":150,
 
                }
'''


def generate_plot_data(motorParams, controlParams):
    
    calcParams={ 
        "BackEMF_VLL_RMS" : motorParams["psiM_Vs_RMS"] * mtr.getElecFreq(1000,motorParams) * np.sqrt(3),
        "BackEMF_VLN_RMS" : motorParams["psiM_Vs_RMS"] * mtr.getElecFreq(1000,motorParams),
        "VoltageMax_LL_RMS" : controlParams["DCVoltage"]/np.sqrt(2),
        "VoltageMax_LN_RMS" : controlParams["DCVoltage"]/np.sqrt(2)/np.sqrt(3)
    }
    
    #update calc Params with some more values
    calcParams["NL_RPM_K"] =(calcParams["VoltageMax_LL_RMS"]/calcParams["BackEMF_VLL_RMS"])
    calcParams["Ke"] = calcParams["BackEMF_VLL_RMS"]/104.72
    calcParams["Kt"] = calcParams["Ke"] * np.sqrt(3)
    calcParams["maxTorque"] = controlParams["RMS_max"]*calcParams["Kt"]
    
    circleDict,ellipseDict,tArrDict = mo.getMotorOperatingRegions(motorParams,controlParams,calcParams,plot_results=0)
    mtpaDict,[c,mtpa_arr] = mo.getMTPA(motorParams,controlParams,calcParams,plot_results=0)
    torqueSpdDict,tsArr = mo.getTorqueSpeedCurve(motorParams,controlParams,calcParams,mtpa_arr,plot_results=0)


