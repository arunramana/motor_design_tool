from dolomites import fnc
from PySide6.QtWidgets import QApplication
import sys, numpy as np
from math import e, pi
from cmath import *

app = QApplication([])
win_r = fnc.fnc_widget() # for the rotating reference frame
win_s = fnc.fnc_widget() # for the stationary reference frame


# in this script we consider some example of space vectors
# trajectories worth to be considered

f    = 50  # the frequncy of the current
t    = np.linspace(0,5/f,501) # we take 5 periods, time discretization
Um = 100
R  = 0.1
w  = 2*np.pi*f
L  =  0.005

# solution in the alpha-beta reference frame
i =-Um/complex(R,w*L)*np.exp(-t/L*R) + Um/complex(R,w*L)*e**(1j*w*t)

# lets convert the solution to the d-q reference frame
i_dq = i*e**(-1j*w*t)

win_s.add_sv(t,i.real,i.imag)
win_r.add_sv_dq(t,i_dq.real,i_dq.imag)

# sincronize the two widgets (stationary reference frame is the master)
win_s.chart_view.time_changed.connect(win_r.chart_view.update_time_line)
# win_r.chart_view.time_changed.connect(win_s.chart_view.update_time_line)

# let's execute the widget with the selected data
win_r.show()
win_s.show()
sys.exit(app.exec())
