'''
apollo-test-alice.py
(https://gitlab.com/LuigiAlberti/dolomites-python)
'''

# import some libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from pathlib import Path

# import apollo
from dolomites import apollo


# load the data from a file
file_name = 'SynRM_data.txt'
# file_name = 'PMA-SynRM_data.txt'

data = pd.read_csv(file_name, sep='\s+', comment='#', names=["theta_m", "i_d", "i_q", "lambda_d", "lambda_q", "torque"])
mot = apollo.fm(data)
mot.create_maps()

# compute the apparent inductances mot.Ld and mot.Lq
mot.calc_apparent_inductances()

# compute the incremental inductances mot.ldd, mot.ldq, mot.lqd, mot.lqd, mot.lqq, mot.lsigma, mot.ldelta
mot.calc_incremental_inductances(method = 'gradient')


# Compute the mtpa trajectory
mot.calc_MTPA(method = "analytical")

# compute saliency mot.xi and estimation error mot.epsilon
mot.calc_saliency()
mot.calc_sensored_error()

# set the MTPA as reference trajectory
mot.i_d_REF = mot.i_d_MTPA
mot.i_q_REF = mot.i_q_MTPA


# compute sensored trajectory t1
mot.calc_sensored_trajectory()

# compute convergence region
Uh=40
fh=1000
mot.calc_convergence_region(Uh, fh)

# compute the inverse inductances
mot.calc_inverse_incremental_inductances() # mot.gamma_dq, mot.gamma_delta

# Fourier

mot.calc_fourier_inverse_inductances(I=3)


mot.calc_fourier_inverse_inductances(I=4,method='fitting')
mot.calc_fourier_inverse_inductances(I=5,method='fft', k=10,plot=True)
mot.calc_Ihq_few_points(I=5,plot=True)

risultati_5_delta = mot.gamma_delta_fourier_coef
risultati_5_dq    = mot.gamma_dq_fourier_coef

print('')
print(risultati_5_delta)
print(risultati_5_dq)





