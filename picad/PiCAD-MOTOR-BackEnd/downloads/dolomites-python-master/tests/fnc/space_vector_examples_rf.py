from dolomites import fnc
from PySide6.QtWidgets import QApplication
import sys, numpy as np
from math import e,pi


app = QApplication([])

win_r = fnc.fnc_widget() # for the rotating reference frame
win_s = fnc.fnc_widget() # for the stationary reference frame

# in this script we consider some example of space vectors
# using a current system as study case.

f    = 50  # the frequncy of the current
wdq  = 2*np.pi*f # speed of the rotating reference frame
Im   = 100 # the peak value
t    = np.linspace(0,2/f,201) # we take 2 periods, time discretization

# Find the three-phase quantities in several cases
case = 7 # case choiche

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

if case == 7: # changing amplitude 3-phase system
  a = Im*(1+t*f)*np.cos(2*np.pi*f*t)
  b = Im*(1+t*f)*np.cos(2*np.pi*f*t - 2/3*np.pi)
  c = Im*(1+t*f)*np.cos(2*np.pi*f*t - 4/3*np.pi)

if case == 8:# changing amplitude and frequency 3-phase system
  t    = np.linspace(0,10/f,2001) # we take 10 periods
  wdq  = 2*np.pi*(f*t*10)
  a = Im*(1+t*f)*np.cos(2*np.pi*t*(f*t*10))
  b = Im*(1+t*f)*np.cos(2*np.pi*t*(f*t*10) - 2/3*np.pi)
  c = Im*(1+t*f)*np.cos(2*np.pi*t*(f*t*10) - 4/3*np.pi)



# generate the sv trajectory in the alpha-beta reference frame
i = 2/3*(a + b*e**(2j*pi/3) + c*e**(4j*pi/3))

# lets convert the solution to the d-q reference frame
i_dq = i*e**(-1j*wdq*t)

win_s.add_sv(t,i.real,i.imag)
win_r.add_sv_dq(t,i_dq.real,i_dq.imag)

# sincronize the two widgets (stationary reference frame is the master)
win_s.chart_view.time_changed.connect(win_r.chart_view.update_time_line)
# win_r.chart_view.time_changed.connect(win_s.chart_view.update_time_line)

# let's execute the widget with the selected data
win_r.show()
win_s.show()
sys.exit(app.exec())
