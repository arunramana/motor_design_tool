import numpy as np
from numpy import pi

#some important trignometric functions

def sin(theta):

    return round(np.sin(np.deg2rad(theta)),10)

def cos(theta):

    return round(np.cos(np.deg2rad(theta)),10)

def tan(theta):

    return round(np.tan(np.deg2rad(theta)),10)

def taninv(theta):

    return round(np.arctan(theta)*180/pi,10)
