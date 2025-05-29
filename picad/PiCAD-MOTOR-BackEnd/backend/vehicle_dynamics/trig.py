import numpy as np
from numpy import pi

#some important trignometric functions

def sin(theta):

    return np.round(np.sin(np.deg2rad(theta)),10)

def cos(theta):

    return np.round(np.cos(np.deg2rad(theta)),10)

def tan(theta):

    return np.round(np.tan(np.deg2rad(theta)),10)

def taninv(theta):

    return np.round(np.arctan(theta)*180/pi,10)
