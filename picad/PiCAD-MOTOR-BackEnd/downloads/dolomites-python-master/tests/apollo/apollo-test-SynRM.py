'''
SynRM example
(https://gitlab.com/LuigiAlberti/dolomites-python)

In this script, simulated SynRM flux linkage maps are post-processed using apollo.
Inductances, MTPA and self-sensing capabilities are computed according with [1].
SynRM flux linkage maps have been obtained with FEMM simulations.
Maxwell torque is available.

[1] Berto, M.; Alberti, L.; Manzolini, V. & Bolognani, S. Computation of Self-Sensing Capabilities of Synchronous Machines for Rotating High Frequency Voltage Injection Sensorless Control IEEE Transactions on Industrial Electronics, 2021.
https://ieeexplore.ieee.org/document/9403927
'''

# import some libraries
import pandas as pd
import matplotlib.pyplot as plt

# import apollo
from dolomites import apollo

# load the data from a file
file_name = 'SynRM_data.txt'

data = pd.read_csv(file_name, sep='\s+', comment='#', names=["theta_m", "i_d", "i_q", "lambda_d", "lambda_q", "torque"])

mot = apollo.fm(data)
mot.create_maps()




# 3d-plot the flux-linkage and torque characteristics (loaded data)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(mot.i_d, mot.i_q, mot.lambda_d)
ax.view_init(elev=45, azim=-135)
ax.set_xlabel('$i_d$ (A)')
ax.set_ylabel('$i_q$ (A)')
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$\lambda_d (Vs)$', rotation=90)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(mot.i_d, mot.i_q, mot.lambda_q)
ax.view_init(elev=45, azim=-135)
ax.set_xlabel('$i_d$ (A)')
ax.set_ylabel('$i_q$ (A)')
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$\lambda_q (Vs)$', rotation=90)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(mot.i_d, mot.i_q, mot.torque)
ax.view_init(elev=45, azim=-135)
ax.set_xlabel('$i_d$ (A)')
ax.set_ylabel('$i_q$ (A)')
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('torque (Nm)', rotation=90)
plt.show()




# OPTIONAL: fit SynRM maps
# mot.fit_flux_SynRM()
# mot.fit_flux_SynRM(Id_min=0, Id_max=6, n_Id=61, Iq_min=0, Iq_max=6, n_Iq=61)
# # OPTIONAL: save the fitted maps
# mot.save_maps(filename = "fitted_maps.txt")



# compute the apparent inductances mot.Ld and mot.Lq
mot.calc_apparent_inductances()

# compute the incremental inductances mot.ldd, mot.ldq, mot.lqd, mot.lqd, mot.lqq, mot.lsigma, mot.ldelta
mot.calc_incremental_inductances(method = 'gradient')


# Compute the mtpa trajectory
# method: choose between "gradient" and "analytical"
mot.calc_MTPA(method = "analytical")

# # Plot the MTPA trajectory
# plt.plot(mot.i_d_MTPA, mot.i_q_MTPA, 'k')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.xlabel('$i_d$ (A)'); plt.ylabel('$i_q$ (A)')
# plt.show()

# import numpy as np
# plt.plot(mot.i_MTPA, mot.theta_MTPA*180/np.pi, 'k')
# plt.xlabel('Stator current amplitude (A)'); plt.ylabel('Stator current angle (deg)')
# plt.show()


# compute saliency mot.xi and estimation error mot.epsilon
mot.calc_saliency()
mot.calc_sensored_error()


# set the MTPA as reference trajectory
mot.i_d_REF = mot.i_d_MTPA
mot.i_q_REF = mot.i_q_MTPA


# compute sensored trajectory t1
mot.calc_sensored_trajectory()

# compute convergence region
Uh = 40
fh = 1000
mot.calc_convergence_region(Uh, fh)
#mot.calc_convergence_region(Uh, fh, Ihq_star=25e-3)


# legend
# black: reference trajectory REF
# yellow: sensored trajectory t1
# blue: Ihq=0 with negative slope (convergence region)
# red: Ihq=0 with poisitive slope (no convergence)
fig, ax = plt.subplots()
circle2 = plt.Circle((0, 0), 2, fill=False, linewidth=0.1)
circle4 = plt.Circle((0, 0), 4, fill=False, linewidth=0.1)
circle6 = plt.Circle((0, 0), 6, fill=False, linewidth=0.1)
circle8 = plt.Circle((0, 0), 8, fill=False, linewidth=0.1)
circle10 = plt.Circle((0, 0), 10, fill=False, linewidth=0.1)
circle12 = plt.Circle((0, 0), 12, fill=False, linewidth=0.1)
ax.add_artist(circle2)
ax.add_artist(circle4)
ax.add_artist(circle6)
ax.add_artist(circle8)
ax.add_artist(circle10)
ax.add_artist(circle12)
ax.plot(mot.i_d_REF, mot.i_q_REF, 'k', label='REF')
ax.plot(mot.i_d_sensored, mot.i_q_sensored, 'y', label='t1')
ax.contour(mot.i_d, mot.i_q, mot.Ihq, levels=[0], colors='r', linestyles='dashed')
ax.contour(mot.i_d, mot.i_q, mot.Ihq_neg, levels=[0], colors='b')
#ax.contour(mot.i_d, mot.i_q, mot.Ihq_in, levels=[0], colors='g', linestyles='dashed')
#ax.contour(mot.i_d, mot.i_q, mot.Ihq_in_neg, levels=[0], colors='m')
plt.xlabel('$i_d$ (A)')
plt.ylabel('$i_q$ (A)')
ax.set_aspect('equal')
plt.title('Convergence region plot')
plt.legend(loc="lower right")
plt.grid(linewidth=0.2)
plt.show()
# LEGEND:
# black: reference trajectory REF
# yellow: sensored trajectory t1
# blue: Ihq=0 with negative slope (convergence region)
# red: Ihq=0 with positive slope (no convergence)


# # Ihq 3D surface plot
#from matplotlib import cm
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(mot.i_d, mot.i_q, mot.Ihq*1000, cmap=cm.coolwarm, lw=0, rstride=1, cstride=1)
#ax.view_init(elev=75, azim=-135)
#ax.contour(mot.i_d, mot.i_q, mot.Ihq, levels=[0], colors="k")
#plt.xlabel('$i_d$ (A)')
#plt.ylabel('$i_q$ (A)')
#ax.zaxis.set_rotate_label(False)
#ax.set_zlabel('$I_{hq}$ (mA)', rotation=90)
#plt.show()


# # OPTIONAL: compute the inverse incremental inductances gamma
# mot.calc_inverse_incremental_inductances()



# # OPTIONAL: export maps for gnuplot
# from pathlib import Path
# Path("gnuplot").mkdir(parents=True, exist_ok=True)
# mot.export_gnuplot_2d_data(mot.lambda_d, 'gnuplot/FluxD.txt')
# mot.export_gnuplot_2d_data(mot.lambda_q, 'gnuplot/FluxQ.txt')
# mot.export_gnuplot_2d_data(mot.Ld, 'gnuplot/Ld.txt')
# mot.export_gnuplot_2d_data(mot.Lq, 'gnuplot/Lq.txt')
# mot.export_gnuplot_2d_data(mot.ldd, 'gnuplot/ldd.txt')
# mot.export_gnuplot_2d_data(mot.ldq, 'gnuplot/ldq.txt')
# mot.export_gnuplot_2d_data(mot.lqd, 'gnuplot/lqd.txt')
# mot.export_gnuplot_2d_data(mot.lqq, 'gnuplot/lqq.txt')
# mot.export_gnuplot_2d_data(mot.lsigma, 'gnuplot/lsigma.txt')
# mot.export_gnuplot_2d_data(mot.ldelta, 'gnuplot/ldelta.txt')
# mot.export_gnuplot_2d_data(mot.torque, 'gnuplot/torque.txt')
# mot.export_gnuplot_1d_data(mot.i_d_MTPA, mot.i_q_MTPA, 'gnuplot/MTPA_cartesian.txt')
# mot.export_gnuplot_1d_data(mot.i_MTPA, mot.theta_MTPA, 'gnuplot/MTPA_polar.txt')
# mot.export_gnuplot_2d_data(mot.xi, 'gnuplot/xi.txt')
# mot.export_gnuplot_2d_data(mot.epsilon_deg, 'gnuplot/epsilon_deg.txt')
# mot.export_gnuplot_1d_data(mot.i_d_REF, mot.i_q_REF, 'gnuplot/REF_cartesian.txt')
# mot.export_gnuplot_1d_data(mot.i_d_sensored, mot.i_q_sensored, 'gnuplot/sensored.txt')
# mot.export_gnuplot_contour_0_data(mot.Ihq, 'gnuplot/Ihq.txt')
# mot.export_gnuplot_contour_0_data(mot.Ihq_neg, 'gnuplot/Ihq_neg.txt')
# mot.export_gnuplot_2d_data(mot.gamma_dd, 'gnuplot/gamma_dd.txt')
# mot.export_gnuplot_2d_data(mot.gamma_dq, 'gnuplot/gamma_dq.txt')
# mot.export_gnuplot_2d_data(mot.gamma_qd, 'gnuplot/gamma_qd.txt')
# mot.export_gnuplot_2d_data(mot.gamma_qq, 'gnuplot/gamma_qq.txt')



# # OPTIONAL: plot HF ellipses
# # (current response to rotating voltage injection)
# import numpy as np
# x = np.linspace(0, 6, 13) # set i_d query points
# y = np.linspace(0, 6, 17) # set i_q query points
# Uh = 40   # set the amplitude of the HF voltage injection
# fh = 1000 # set the frequency of the HF voltage injection
# k = 1     # set a scaling factor (optional)
# # call plot_ellipses function
# ellipse_plot = mot.plot_ellipses(x, y, Uh, fh, k)
# plt.xlabel('$i_d$ (A)')
# plt.ylabel('$i_q$ (A)')
# plt.show()
