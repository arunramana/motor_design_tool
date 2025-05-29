from dolomites import fnc
from PySide6.QtWidgets import QApplication
import sys, numpy as np

app = QApplication([])
win = fnc.fnc_widget()

# in this script we consider some example of space vectors
# using a current system as study case.

f    = 50  # the frequncy of the current
Im   = 100 # the peak value
t    = np.linspace(0,2/f,201) # we take 2 periods, time discretization

# Find the three-phase quantities in several cases
case = 10 # case choiche

if case == 1: # balsnced symmetrical system (positive sequence)
  a = Im*np.cos(2*np.pi*f*t)
  b = Im*np.cos(2*np.pi*f*t - 2/3*np.pi)
  c = Im*np.cos(2*np.pi*f*t - 4/3*np.pi)

if case == 2: # balanced symmetrical system (negative sequence)
  a = Im*np.cos(2*np.pi*f*t)
  b = Im*np.cos(2*np.pi*f*t + 2/3*np.pi)
  c = Im*np.cos(2*np.pi*f*t + 4/3*np.pi)

if case == 3: # fundamental and 5th order harmonic
  a = Im*np.cos(2*np.pi*f*t)             + 0.2*Im*np.sin(5*2*np.pi*f*t)
  b = Im*np.cos(2*np.pi*f*t - 2/3*np.pi) + 0.2*Im*np.sin(5*2*np.pi*f*t- 5*2/3*np.pi)
  c = Im*np.cos(2*np.pi*f*t - 4/3*np.pi) + 0.2*Im*np.sin(5*2*np.pi*f*t- 5*4/3*np.pi)

if case == 4: # fundamental and 7th order harmonic
  a = Im*np.cos(2*np.pi*f*t)             + 0.2*Im*np.sin(7*2*np.pi*f*t)
  b = Im*np.cos(2*np.pi*f*t - 2/3*np.pi) + 0.2*Im*np.sin(7*2*np.pi*f*t- 7*2/3*np.pi)
  c = Im*np.cos(2*np.pi*f*t - 4/3*np.pi) + 0.2*Im*np.sin(7*2*np.pi*f*t- 7*4/3*np.pi)

if case == 5: # balanced not-symmetrical system
              # (sum of positive and negative sequence)
  a =     Im*np.cos(2*np.pi*f*t)
  b = 0.2*Im*np.cos(2*np.pi*f*t - 2/3*np.pi)
  c = -a-b #just to force ia+ib+ic=0

if case == 6: # fundamental and 3rd order harmonic
  a = Im*(np.cos(2*np.pi*f*t))             + 0.2*Im*np.sin(3*2*np.pi*f*t)
  b = Im*(np.cos(2*np.pi*f*t - 2/3*np.pi)) + 0.2*Im*np.sin(3*2*np.pi*f*t- 3*2/3*np.pi)
  c = Im*(np.cos(2*np.pi*f*t - 4/3*np.pi)) + 0.2*Im*np.sin(3*2*np.pi*f*t- 3*4/3*np.pi)

if case == 7: # square waves
  a = Im*np.sign(np.cos(2*np.pi*f*t))
  b = Im*np.sign(np.cos(2*np.pi*f*t - 2/3*np.pi))
  c = Im*np.sign(np.cos(2*np.pi*f*t - 4/3*np.pi))

if case == 8: # six-step waveform
  t  = np.linspace(0,1/f,121)
  a = np.zeros(121)
  b = np.zeros(121)
  c = np.zeros(121)

  a[0:20]    =    Im/3
  a[20:40]   =  2*Im/3
  a[40:60]   =    Im/3
  a[60:80]   =   -Im/3
  a[80:100]  = -2*Im/3
  a[100:121] =   -Im/3

  b[40:60]   =    Im/3
  b[60:80]   =  2*Im/3
  b[80:100]  =    Im/3
  b[100:121] =   -Im/3
  b[0:20]    = -2*Im/3
  b[20:40]   =   -Im/3

  c[80:100]  =    Im/3
  c[100:121] =  2*Im/3
  c[0:20]    =    Im/3
  c[20:40]   =   -Im/3
  c[40:60]   = -2*Im/3
  c[60:80]   =   -Im/3

if case == 9: # changing amplitude 3-phase system
  a = Im*(1+t*f)*np.cos(2*np.pi*f*t)
  b = Im*(1+t*f)*np.cos(2*np.pi*f*t - 2/3*np.pi)
  c = Im*(1+t*f)*np.cos(2*np.pi*f*t - 4/3*np.pi)

if case == 10:# changing amplitude and frequency 3-phase system
  t    = np.linspace(0,10/f,2001) # we take 10 periods
  a = Im*(1+t*f)*np.cos(2*np.pi*t*(f*t*10))
  b = Im*(1+t*f)*np.cos(2*np.pi*t*(f*t*10) - 2/3*np.pi)
  c = Im*(1+t*f)*np.cos(2*np.pi*t*(f*t*10) - 4/3*np.pi)

# let's execute the widget with the selected data
win.add_series(t,a,b,c)
win.show()
sys.exit(app.exec())
