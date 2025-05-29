# [1] Berto, M.; Alberti, L.; Manzolini, V. & Bolognani, S.
# Computation of Self-Sensing Capabilities of Synchronous Machines
# for Rotating High Frequency Voltage Injection Sensorless Control
# IEEE Transactions on Industrial Electronics, 2021.
# https://ieeexplore.ieee.org/document/9403927

# [2] Tinazzi, F.; Bolognani, S.; Calligaro, S.; Kumar, P.; Petrella, R. & Zigliotto, M.
# Classification and review of MTPA algorithms for synchronous reluctance and interior permanent magnet motor drives 
# 2019 21st European Conference on Power Electronics and Applications (EPE '19 ECCE Europe), 2019, P.1-P.10 

# [3] Varatharajan, A.; Cruz, S.; Hadla, H. & Briz, F.
# Predictive torque control of SynRM drives with online MTPA trajectory tracking and inductances estimation 
# 2017 IEEE International Electric Machines and Drives Conference (IEMDC), 2017, 1-7 

# [4] M. Hinkkanen, P. Pescetto, E. Mölsä, S. E. Saarakkala, G. Pellegrino and R. Bojoi,
# "Sensorless Self-Commissioning of Synchronous Reluctance Motors at Standstill Without
# Rotor Locking," in IEEE Transactions on Industry Applications, vol. 53, no. 3, pp.
# 2120-2129, May-June 2017, doi: 10.1109/TIA.2016.2644624.

# [5] M. Berto, L. Alberti, F. Martin and M. Hinkkanen,
# "Online Incremental Inductance Identification for Reluctance Synchronous Motors"
# IECON 2021 - 47th Annual Conference of the IEEE Industrial Electronics Society, 2021, p. 1-6
# https://ieeexplore.ieee.org/document/9589537



'''
LIST OF FUNCTIONS:
self.create_maps()
self.save_maps()
self.export_gnuplot_1d_data()
self.export_gnuplot_2d_data()
self.export_gnuplot_contour_0_data()
self.calc_apparent_inductances()
self.calc_incremental_inductances()
self.calc_inverse_incremental_inductances()
self.calc_MTPA()
self.calc_saliency()
self.calc_sensored_error()
self.calc_sensored_trajectory()
self.calc_convergence_region()
self.plot_ellipses()
self.fit_flux_SynRM()
self.calc_fourier_inverse_inductances()
self.calc_Ihq_few_points()
'''


'''
LIST OF VARIABLES:
self.data
self.i_d
self.i_q
self.lambda_d
self.lambda_q
self.torque
self.lambda_d0
self.lambda_q0
self.Ld
self.Lq
self.ldd
self.ldq
self.lqd
self.lqq
self.lsigma
self.ldelta
self.lneg
self.gamma_dd
self.gamma_dq
self.gamma_qd
self.gamma_qq
self.gamma_sigma
self.gamma_delta
self.gamma_neg
self.torque_d_alpha_ie
self.i_d_MTPA
self.i_q_MTPA
self.i_MTPA
self.theta_MTPA
self.xi
self.epsilon
self.epsilon_deg
self.i_REF
self.theta_REF
self.i_d_REF
self.i_q_REF
self.epsilon_REF
self.i_d_sensored
self.i_q_sensored
self.delta_theta
self.Ihq
self.Ihq_d_alpha_ie
self.Ihq_neg
self.Ihq_in
self.Ihq_in_d_alpha_ie
self.Ihq_in_neg
self.Uh
self.fh
self.alpha_breakpoints
self.gamma_delta_fourier_coef
self.gamma_delta_fourier_rec
self.gamma_dq_fourier_coef
self.gamma_dq_fourier_rec
self.Ihq_polar_I
self.Ihq_fourier
self.Ihq_fourier_rec
self.alpha_breakpoints_NaN
self.gamma_delta_a0
self.gamma_delta_a2
self.gamma_dq_b2
self.Ihq_fourier_2meas
self.gamma_dq_b1
'''



import copy
import numpy as np
import numpy.matlib
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.interpolate import griddata
from scipy.interpolate import interpn
from scipy.fftpack import fft
from scipy.optimize import curve_fit
from scipy.spatial import Delaunay
import warnings
warnings.filterwarnings("ignore")


print('imported apollo (rev 14 January 2024)')


class fm():
    """
    A class to represent a flux map
    fm is a pandas dataset which represents a set of simulations/measures on an electric machine
    The dataset is expected to have the following fields required for the subsequent computations:
    ["theta_m", "i_d", "i_q", "lambda_d", "lambda_q", "torque"]
    Torque is expected to be the "Maxwell" torque.
    """
        
    def __init__(self, data):
        self.data = data  # a pandas dataframe used to store all the available data
      
        # some pivot tables to re-organize the data
        self._lambda_d = pd.pivot_table(self.data, index=['theta_m', 'i_d'], columns=['i_q'], values=['lambda_d'])
        self._lambda_q = pd.pivot_table(self.data, index=['theta_m', 'i_d'], columns=['i_q'], values=['lambda_q'])
        self._torque   = pd.pivot_table(self.data, index=['theta_m', 'i_d'], columns=['i_q'], values=['torque'])
       
    def __repr__(self):

        return self.data.to_string()
        
    def create_maps(self, _theta_m=0):
        """
        Compute the maps for the given rotor position
        The following maps are created as numpy arrays:
        lambda_d(i_d,i_q)
        lambda_q(i_d,i_q)
        torque(i_d,i_q)
        
        After the call of this function all the quantities are matrices
        (lambda_d, lambda_q, torque, i_d, i_q)
        """

        try:
            lambda_d = self._lambda_d.xs(_theta_m)
        except LookupError:
            print('  APOLLO ERROR in create_maps:')
            print('  the given rotor position (%s) is not in the dataset!' %(_theta_m))
            exit()
            return

        self.i_d = lambda_d.index.get_level_values('i_d').to_numpy()
        self.i_q = lambda_d.columns.get_level_values('i_q').to_numpy()
        
        ## hereafter i_d and i_q will represent matrices so that we can manage scattered data
        self.i_d, self.i_q = np.meshgrid(self.i_d, self.i_q)
        self.i_d = self.i_d.transpose()
        self.i_q = self.i_q.transpose()
        
        self.lambda_d = lambda_d.to_numpy()
        
        self.lambda_q = self._lambda_q.xs(_theta_m).to_numpy()
        self.torque = self._torque.xs(_theta_m).to_numpy()
        # Fix for self.torque:
        # count the number of zero- and nan-elements
        counted_zeros = self.torque[np.where(self.torque == 0)].size
        counted_nans = self.torque[np.isnan(self.torque)].size
        # if the sum of zero- and nan-elements is equal to the whole size,
        # generate self.torque from scretch
        if counted_zeros + counted_nans == self.torque.size:
            self.torque = np.zeros(self.i_d.shape)
        
        print('created maps for theta_m = %s' % (_theta_m))
        
        
        ## When a map is scattered it may contain NaNs which have to be interpolated
        ## If any NaN is in lambda_d (and probably so also in the other data) fill the missing values
        if (np.any(np.isnan(self.lambda_d))):

            ## Prepare arrays for griddata: i_d_vec, i_q_vec, lambda_d_vec, lambda_q_vec and torque_vec
            # Convert to one dimensional arrays for griddata
            i_d_vec = self.i_d.flatten()
            i_q_vec = self.i_q.flatten()
            lambda_d_vec = self.lambda_d.flatten()
            lambda_q_vec = self.lambda_q.flatten()
            torque_vec = self.torque.flatten()
            
            # Remove NaNs for griddata
            NaN_mask = ~np.isnan(lambda_d_vec) # indices of the non-NaNs
            i_d_vec = i_d_vec[NaN_mask] # apply mask to keep the non-NaNs, i.e. remove NaNs
            i_q_vec = i_q_vec[NaN_mask]
            lambda_d_vec = lambda_d_vec[NaN_mask]
            lambda_q_vec = lambda_q_vec[NaN_mask]
            torque_vec = torque_vec[NaN_mask]
            
            ## Fill missing values INSIDE the convex hull with griddata
            ## Interpolate with griddata on the original (self.i_d, self.i_q) grid
            points_I = np.vstack((i_d_vec, i_q_vec)).T
            self.lambda_d = griddata(points_I, lambda_d_vec, (self.i_d, self.i_q), method='cubic')
            self.lambda_q = griddata(points_I, lambda_q_vec, (self.i_d, self.i_q), method='cubic')
            self.torque   = griddata(points_I, torque_vec,   (self.i_d, self.i_q), method='cubic')
            
            print('filled missing values in the maps')
     
      
     
     
    def save_maps(self, theta_m=0, filename="apollo_maps.txt"):
        """
        Export maps (for the chosen rotor position) in a text format compatible with apollo.
        If not specified, the output rotor position is 0.
        This is useful to save maps after elaborating them (for example with
        interpolation) for further use
        """
        # flatten output arrays
        theta_m = np.ones(self.i_d.shape) * theta_m
        theta_m = theta_m.flatten()
        i_d = self.i_d.flatten()
        i_q = self.i_q.flatten()
        lambda_d = self.lambda_d.flatten()
        lambda_q = self.lambda_q.flatten()
        torque = self.torque.flatten()
        # define output matrix
        mat_out = np.transpose([theta_m, i_d, i_q, lambda_d, lambda_q, torque])
        # save text
        np.savetxt(filename, mat_out, fmt='%1.10e')
        print('saved maps as "', filename, '"')
     
     
    def export_gnuplot_1d_data(self, vector_1, vector_2, filename):
        """
        This function exports 1d maps in *.dat files compatible with gnuplot
        1d maps are for example the trajectories (REF, MTPA and t1)
        """
        mat = np.vstack((vector_1, vector_2)).T
        np.savetxt(filename, mat, fmt='%1.10e')
        print('created "', filename, '"')

    def export_gnuplot_2d_data(self, map_2d, filename):
        """
        This function exports 2d maps in *.dat files compatible with gnuplot
        (add a blank line at the end of every row)
        2d maps are the fluxes, the torque, the inductances, the inverse incremental inductances...
        """
        X = self.i_d.flatten()
        Y = self.i_q.flatten()
        Z = map_2d.flatten()
        mat = np.transpose([X, Y, Z])

        # stackoverflow 752308
        ## Split data appending infs
        ## The infs will be removed later in order to obtain a blank line
        def split(arr, size):
            arrs = []
            while len(arr) > size:
                pice = arr[:size]
                arrs.append(pice)
                arrs.append(np.array([np.inf, np.inf, np.inf]))
                arr = arr[size:]
            arrs.append(arr)
            return arrs

        np.savetxt(filename, np.vstack(split(mat, self.i_q.shape[1])), fmt='%1.10e')

        # stackoverflow 17140886
        ## Open the saved text file. Replace infs with blanks
        ## The modified file is compatible with gnuplot
        # Read in the file
        with open(filename, 'r') as file:
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('inf', '')
        # Write the file out again
        with open(filename, 'w') as file:
            file.write(filedata)
            
        print('created "', filename, '"')

    def export_gnuplot_contour_0_data(self, map_2d, filename):
        """
        This function exports the isolevel 0 of 2d maps in *.dat files compatible with gnuplot
        It can be used on every 2d map. The level 0 is extracted from matplotlib contour plot.
        """
        cn = plt.contour(self.i_d, self.i_q, map_2d, levels=[0])
        plt.close()
        # ectract the coordinates of all the lines (closed lines and open lines)
        for i in range(len(cn.allsegs[0])):
            modified_filename = filename.replace(".dat", "")
            modified_filename = modified_filename + "_" + str(i) + ".dat"
            # print(cn.allsegs[1][i])
            np.savetxt(modified_filename, cn.allsegs[0][i], fmt='%1.10e')

    def calc_apparent_inductances(self):
        """
        Starting from the flux linkages, compute the apparent inductances
        Refer to the equations in the example
        """
        
        # evaluate lambda_d for id=0
        if 0 in self.i_d: # use the data in Lambda_d if available
            index_id0 = np.where(self.i_d[:,0] == 0)
            lambda_d0 = self.lambda_d[index_id0[0], :]
            lambda_d0 = lambda_d0.flatten()
        else: # otherwise estimate the values. In this case we need negative d-axis currents
            i_d_vec = self.i_d.flatten()
            i_q_vec = self.i_q.flatten()
            lambda_d_vec = self.lambda_d.flatten()
            i_d_vec = i_d_vec[~np.isnan(lambda_d_vec)]
            i_q_vec = i_q_vec[~np.isnan(lambda_d_vec)]
            lambda_d_vec = lambda_d_vec[~np.isnan(lambda_d_vec)]
            points_I = np.vstack((i_d_vec, i_q_vec)).T
            lambda_d0 = griddata(points_I, lambda_d_vec, (0, self.i_q[0,:]), method='cubic')

        # Fix for lambda_d0:
        if np.isnan(lambda_d0).all():
            lambda_d0 = np.zeros(lambda_d0.shape)

        # make lambda_d0 a 2D array
        self.lambda_d0 = np.matlib.repmat(lambda_d0, self.i_d.shape[0], 1)
        
        
        # evaluate lambda_q for iq=0
        if 0 in self.i_q: # use the data in Lambda_q if available
            index_iq0 = np.where(self.i_q[0,:] == 0)
            lambda_q0 = self.lambda_q[:, index_iq0[0]]
            lambda_q0 = lambda_q0.flatten()
        else: # otherwise estimate the values. In this case we need negative q-axis currents
            i_d_vec = self.i_d.flatten()
            i_q_vec = self.i_q.flatten()
            lambda_q_vec = self.lambda_q.flatten()
            i_d_vec = i_d_vec[~np.isnan(lambda_q_vec)]
            i_q_vec = i_q_vec[~np.isnan(lambda_q_vec)]
            lambda_q_vec = lambda_q_vec[~np.isnan(lambda_q_vec)]
            points_I = np.vstack((i_d_vec, i_q_vec)).T
            lambda_q0 = griddata(points_I, lambda_q_vec, (self.i_d[:,0], 0), method='cubic')
        
        # Fix for lambda_q0:
        if np.isnan(lambda_q0).all():
            lambda_q0 = np.zeros(lambda_q0.shape)

        # make lambda_q0 a 2D array
        self.lambda_q0 = np.matlib.repmat(np.array([lambda_q0]).T, 1, self.i_q.shape[1])


        # compute Ld and Lq
        self.Ld = (self.lambda_d - self.lambda_d0) / self.i_d
        self.Lq = (self.lambda_q - self.lambda_q0) / self.i_q

        
        # If any NaN is in Ld fill the missing values (fill the discontinuities)
        if (np.any(np.isnan(self.Ld))):
            # Prepare arrays for griddata: i_d_vec, i_q_vec, Ld_vec
            # Convert to one dimensional arrays for griddata
            i_d_vec = self.i_d.flatten()
            i_q_vec = self.i_q.flatten()
            Ld_vec = self.Ld.flatten()
            # Remove NaNs for griddata
            NaN_mask = ~np.isnan(Ld_vec) # indices of the non-NaNs
            i_d_vec = i_d_vec[NaN_mask] # apply mask to keep the non-NaNs, i.e. remove NaNs
            i_q_vec = i_q_vec[NaN_mask]
            Ld_vec = Ld_vec[NaN_mask]
            # Interpolate with griddata on the original (self.i_d, self.i_q) grid
            points_I = np.vstack((i_d_vec, i_q_vec)).T
            self.Ld = griddata(points_I, Ld_vec, (self.i_d, self.i_q), method='cubic')
        
        # If any NaN is in the first row of Ld (Id=0 axis) fill the missing values
        if (np.all(np.isnan(self.Ld[0,:]))):
            # copy values from second row
            self.Ld[0,:] = self.Ld[1,:]
            
        # If any NaN is in the last row of Ld (Id=0 axis) fill the missing values
        if (np.all(np.isnan(self.Ld[-1,:]))):
            # copy values from second to last row
            self.Ld[-1,:] = self.Ld[-2,:]
        
        # If any NaN is in Lq fill the missing values (fill the discontinuities)        
        if (np.any(np.isnan(self.Lq))):
            # Prepare arrays for griddata: i_d_vec, i_q_vec, Ld_vec
            # Convert to one dimensional arrays for griddata
            i_d_vec = self.i_d.flatten()
            i_q_vec = self.i_q.flatten()
            Lq_vec = self.Lq.flatten()
            # Remove NaNs for griddata
            NaN_mask = ~np.isnan(Lq_vec) # indices of the non-NaNs
            i_d_vec = i_d_vec[NaN_mask] # apply mask to keep the non-NaNs, i.e. remove NaNs
            i_q_vec = i_q_vec[NaN_mask]
            Lq_vec = Lq_vec[NaN_mask]
            # Interpolate with griddata on the original (self.i_d, self.i_q) grid
            points_I = np.vstack((i_d_vec, i_q_vec)).T
            self.Lq = griddata(points_I, Lq_vec, (self.i_d, self.i_q), method='cubic')
        
        # If any NaN is in the first column of Lq (Iq=0 axis) fill the missing values
        if (np.all(np.isnan(self.Lq[:,0]))):
            # copy values from second column
            self.Lq[:,0] = self.Lq[:,1]
        
        # If any NaN is in the last column of Lq (Iq=0 axis) fill the missing values
        if (np.all(np.isnan(self.Lq[:,-1]))):
            # copy values from second to last column
            self.Lq[:,-1] = self.Lq[:,-2]


        print('computed apparent inductances')

    def calc_incremental_inductances(self, method='gradient'):
        """
        Starting from the flux linkages, compute the incremental inductances
        with the given method:
            -) gradient: use the numpy gradient algorithm
            -) diff: use the numpy diff algorithm. Last row and last column of
               the inductances matrices are duplicated in order the mantain the
               same dimension of the flux maps.
        NOTE: 'diff' calculates the differences between adjacent elements, while 'gradient'
        computes the 'central differences'.
        -> 'gradient' should be preferred to make convergence region symmetrical.
        """

        if method == 'gradient':
            (self.ldd, self.ldq) = np.gradient(self.lambda_d, self.i_d[:,0], self.i_q[0,:])
            (self.lqd, self.lqq) = np.gradient(self.lambda_q, self.i_d[:,0], self.i_q[0,:])
            print('computed incremental inductances using gradient')

        elif method == 'diff':
            # compute incremental inductances using diff
            
            self.ldd = ((np.diff(self.lambda_d, axis=0)).T / np.diff(self.i_d[:,0])).T
            # duplicate last row
            self.ldd = np.concatenate((self.ldd, np.matrix(self.ldd[-1, :])), axis=0)
            # numpy matrix to array
            self.ldd = np.array(self.ldd)

            self.ldq = ((np.diff(self.lambda_d, axis=1)) / np.diff(self.i_q[0,:]))
            # duplicate last column
            self.ldq = np.concatenate((self.ldq, np.matrix(self.ldq[:, -1]).T), axis=1)
            # numpy matrix to array
            self.ldq = np.array(self.ldq)

            self.lqd = ((np.diff(self.lambda_q, axis=0)).T / np.diff(self.i_d[:,0])).T
            # duplicate last row
            self.lqd = np.concatenate((self.lqd, np.matrix(self.lqd[-1, :])), axis=0)
            # numpy matrix to array
            self.lqd = np.array(self.lqd)

            self.lqq = ((np.diff(self.lambda_q, axis=1)) / np.diff(self.i_q[0,:]))
            # duplicate last column
            self.lqq = np.concatenate((self.lqq, np.matrix(self.lqq[:, -1]).T), axis=1)
            # numpy matrix to array
            self.lqq = np.array(self.lqq)
            print('computed incremental inductances using diff')

        else:
            print('  APOLLO ERROR in calc_incremental_inductances:')
            print('  impossible to compute the incremental inductances')
            print('  please use "diff" or "gradient" method')
            exit()

        self.lsigma = (self.lqq + self.ldd) / 2
        self.ldelta = (self.lqq - self.ldd) / 2

        # the incremental inductance lsigma is known to be related to the positive-sequence current response.
        # on the other hand, the negative-sequence incremental inductance is [5]:
        self.lneg = np.sqrt(self.ldelta**2 + self.ldq**2)

    def calc_inverse_incremental_inductances(self):
        """
        Starting from the incremental inductances, compute the inverse incremental inductances (gamma)
        """

        # check if incremental inductances have been computed
        try:
            d = np.multiply(self.ldd, self.lqq) - np.multiply(self.ldq, self.lqd)
        except AttributeError:
            print('  APOLLO ERROR in calc_inverse_incremental_inductances:')
            print('  incremental inductances not yet computed!')
            print('  consider to call calc_incremental_inductances() first')
            exit()
            return
        else:
            self.gamma_dd = np.divide(self.lqq, d)
            self.gamma_dq = np.divide(-self.ldq, d)
            self.gamma_qd = np.divide(-self.lqd, d)
            self.gamma_qq = np.divide(self.ldd, d)
            self.gamma_sigma = (self.gamma_dd + self.gamma_qq) / 2
            self.gamma_delta = (self.gamma_dd - self.gamma_qq) / 2
            self.gamma_neg = np.sqrt(self.gamma_delta**2 + self.gamma_dq**2)
            
        print('computed inverse incremental inductances')

    def calc_MTPA(self, method="gradient", quadrant=None):
        """
        Compute the MTPA locus, both in cartesian and polar coordinates, in the
        selected quadrant.
        Two different methods to compute MTPA can be choosen:
        "gradient": tangential gradient of the FEA (Maxwell) torque, if available
        "analytical": incremental inductances equation, based on dq torque
        """
        
        ## if quadrant is not set, auto-detect MTPA quadrant
        if quadrant is None:
            
            # check if incremental inductances have been computed
            try:
                d = np.multiply(self.ldd, self.lqq) - np.multiply(self.ldq, self.lqd)
            except AttributeError:
                print('  APOLLO ERROR in calc_MTPA:')
                print('  impossible to auto-detect the MTPA quadrant')
                print('  incremental inductances not yet computed!')
                print('  consider to call calc_incremental_inductances() first')
                print('  otherwise, set the MTPA quadrant manually')
                exit()
                return
            else:
            
                # auto-detect convention, and then set the MTPA quadrant
                if np.mean(self.ldelta[~np.isnan(self.ldelta)]) > 0:
                    # detected IPM/PMA-SynRM convention. MTPA will be computed on quadrant 2
                    quadrant = 2
                    print('auto-detected MTPA quadrant using IPM convention')
                else:
                    # detected IPM/PMA-SynRM convention. MTPA will be computed on quadrant 1
                    quadrant = 1
                    print('auto-detected MTPA quadrant using SynRM convention')
            
                   
            
        ## compute torque derivative self.torque_d_alpha_ie
        if method == "gradient":
        
            ## test if self.torque contains only zeros or NaNs
            if np.any(np.nan_to_num(self.torque, nan=0)):
                torque_dx = np.gradient(self.torque, self.i_d[:,0], axis=0)
                torque_dy = np.gradient(self.torque, self.i_q[0,:], axis=1)
                self.torque_d_alpha_ie = - torque_dx * self.i_q + torque_dy * self.i_d
            
            else:
                print('  APOLLO ERROR in calc_MTPA:')
                print('  impossible to compute MTPA on quadrant', quadrant, 'using gradient')
                print('  the FEA (Maxwell) torque is not available!')
                print('  consider to use analytical method instead')
                exit()

            
        elif method == "analytical":
            # check if incremental inductances have been computed
            try:
                self.ldd
            except AttributeError:
                print('  APOLLO ERROR in calc_MTPA:')
                print('  impossible to compute MTPA on quadrant', quadrant, 'using analythical')
                print('  incremental inductances not yet computed!')
                print('  consider to call calc_incremental_inductances() first')
                exit()
                return
            else:
                # compute torque derivative using equation (11) by prof Bolognani [2],
                # based on equation (10) by prof Briz [3]
                self.torque_d_alpha_ie = ((self.ldq + self.lqd) * self.i_d) * self.i_q \
                    - (self.ldd * self.i_q**2 + (self.lqq * self.i_d**2)) \
                    + (self.lambda_d * self.i_d) + self.lambda_q * self.i_q

        else:
            print('  APOLLO ERROR in calc_MTPA:')
            print('  impossible to compute the MTPA trajectory')
            print('  please use method "gradient" or "analytical"')
            exit()

        # contour of self.torque_d_alpha_ie=0 (MTPA region)
        cn = plt.contour(self.i_d, self.i_q, self.torque_d_alpha_ie, levels=[0])
        plt.close() # don't show contour plot

        # ectract coordinates of the MTPA region from contour data
        C = np.vstack(cn.allsegs[0])
        Id_MTPA = C[:, 0]
        Iq_MTPA = C[:, 1]

        # occasional out of bounds fix
        Id_MTPA = np.clip(Id_MTPA, np.min(self.i_d[:,0]), np.max(self.i_d[:,0]))
        Iq_MTPA = np.clip(Iq_MTPA, np.min(self.i_q[0,:]), np.max(self.i_q[0,:]))

        if quadrant == 1:
            # select quadrants 1 and 2
            Id_MTPA = Id_MTPA[Iq_MTPA > 0]
            Iq_MTPA = Iq_MTPA[Iq_MTPA > 0]
            # select quadrant 1
            Iq_MTPA = Iq_MTPA[Id_MTPA > 0]
            Id_MTPA = Id_MTPA[Id_MTPA > 0]
            print('computed MTPA on quadrant 1 using', method)

        elif quadrant == 2:
            # select quadrants 1 and 2
            Id_MTPA = Id_MTPA[Iq_MTPA > 0]
            Iq_MTPA = Iq_MTPA[Iq_MTPA > 0]
            # select quadrant 2
            Iq_MTPA = Iq_MTPA[Id_MTPA < 0]
            Id_MTPA = Id_MTPA[Id_MTPA < 0]
            print('computed MTPA on quadrant 2 using', method)

        elif quadrant == 3:
            # select quadrants 3 and 4
            Id_MTPA = Id_MTPA[Iq_MTPA < 0]
            Iq_MTPA = Iq_MTPA[Iq_MTPA < 0]
            # select quadrant 3
            Iq_MTPA = Iq_MTPA[Id_MTPA < 0]
            Id_MTPA = Id_MTPA[Id_MTPA < 0]
            print('computed MTPA on quadrant 3 using', method)

        elif quadrant == 4:
            # select quadrants 3 and 4
            Id_MTPA = Id_MTPA[Iq_MTPA < 0]
            Iq_MTPA = Iq_MTPA[Iq_MTPA < 0]
            # select quadrant 4
            Iq_MTPA = Iq_MTPA[Id_MTPA > 0]
            Id_MTPA = Id_MTPA[Id_MTPA > 0]
            print('computed MTPA on quadrant 4 using', method)

        else:
            print('  APOLLO ERROR in calc_MTPA:')
            print('  wrong quadrant for MTPA computation!')
            print('  set a quadrant between 1 and 4')
            exit()


        # if the coordinates of the extracted MTPA trajectory are empty,
        # select another quadrant
        if (Id_MTPA.size == 0) and (Iq_MTPA.size == 0):
            print('  APOLLO ERROR in calc_MTPA:')
            print('  the considered motor has no MTPA on quadrant ', quadrant)
            print('  please select another quadrant')
            exit()


        # transform cartesian coordinates to polar
        I_MTPA = np.sqrt(Id_MTPA**2 + Iq_MTPA**2)
        theta_MTPA = np.arctan2(Iq_MTPA, Id_MTPA)
        theta_MTPA = theta_MTPA % (2 * np.pi)  # modulus at 2 pi

        # sort vectors with increasing I_MTPA
        order = np.argsort(I_MTPA)
        I_MTPA = I_MTPA[order]
        theta_MTPA = theta_MTPA[order]
        Id_MTPA = Id_MTPA[order]
        Iq_MTPA = Iq_MTPA[order]
        
        # add 0,0 point
        Id_MTPA = np.hstack((0, Id_MTPA))
        Iq_MTPA = np.hstack((0, Iq_MTPA))
        # compute MTPA angle in 0,0 with linear extrapolation
        y = theta_MTPA[0] + (0-I_MTPA[0])/(I_MTPA[-1]-I_MTPA[0]) * (theta_MTPA[-1]-theta_MTPA[0])
        theta_MTPA = np.hstack((y, theta_MTPA))
        I_MTPA = np.hstack((0, I_MTPA))

        # save the results
        self.i_d_MTPA = Id_MTPA
        self.i_q_MTPA = Iq_MTPA
        self.i_MTPA = I_MTPA
        self.theta_MTPA = theta_MTPA

    def calc_saliency(self):
        """
        Compute the hf-saliency (xi)
        """

        # calc hf-saliency
        self.xi = np.divide((self.lsigma + np.sqrt(self.ldelta**2 + self.ldq**2)), (self.lsigma - np.sqrt(self.ldelta**2 + self.ldq**2)))
        print('computed saliency')
                            
    def calc_sensored_error(self):
        """
        Compute the estimation error when the observer is in open loop (epsilon)
        The convention (SynRM or IPM/PMA-SynRM) is automatically detected by default
        """
        # calc estimation error in electrical radiants
        if np.mean(self.ldelta[~np.isnan(self.ldelta)]) > 0:
            # use IPM/PMA-SynRM convention
            self.epsilon = 0.5 * np.arctan2(-self.ldq, self.ldelta)
            print('computed sensored error using IPM convention')
        else:
            # use SynRM convention
            self.epsilon = 0.5 * np.arctan2(self.ldq, -self.ldelta)
            print('computed sensored error using SynRM convention')

        # compute estimation error also in electrical degrees
        self.epsilon_deg = self.epsilon * 180 / np.pi
                
        
    def calc_sensored_trajectory(self):
        """
        Compute the trajectory of the estimated position during a sensored test
        (estimator in open loop), given a current reference trajectory.
        The sensored trajectory is labeled as "t1" in [1]
        """
        
        # compute i_d_REF and i_q_REF if they are not defined
        if (not hasattr(self, 'i_d_REF')) or (not hasattr(self, 'i_q_REF')):
            try:
                self.i_d_REF = self.i_REF * np.cos(self.theta_REF)
                self.i_q_REF = self.i_REF * np.sin(self.theta_REF)
            except :
                print('  APOLLO ERROR in calc_sensored_trajectory:')
                print('  reference trajectory seems not assigned')
                print('  consider to assign i_REF and theta_REF values (or i_d_REF and i_q_REF) to fm object first')
                exit()
                return
                    
        # compute i_REF and theta_REF if they are not defined
        if (not hasattr(self, 'i_REF')) or (not hasattr(self, 'theta_REF')):
            try:
                self.i_REF = np.sqrt(self.i_d_REF**2 + self.i_q_REF**2)
                self.theta_REF = np.arctan2(self.i_q_REF, self.i_d_REF)
            except :
                print('  APOLLO ERROR in calc_sensored_trajectory:')
                print('  reference trajectory seems not assigned')
                print('  consider to assign i_REF and theta_REF values (or i_d_REF and i_q_REF) to fm object first')
                exit()
                return
                
        # sort vectors with increasing i_REF
        order = np.argsort(self.i_REF)
        self.i_REF = self.i_REF[order]
        self.theta_REF = self.theta_REF[order]
        self.i_d_REF = self.i_d_REF[order]
        self.i_q_REF = self.i_q_REF[order]

        # FIX: Cut the REF trajectory if it oversteps the convex hull
        # first, prepare data to compute the convex hull
        # flatten
        i_d_vec = self.i_d.flatten()
        i_q_vec = self.i_q.flatten()
        lambda_d_vec = self.lambda_d.flatten()
        # remove nans
        NaN_mask = ~np.isnan(lambda_d_vec) # indices of the non-NaNs
        i_d_vec = i_d_vec[NaN_mask] # apply mask to keep the non-NaNs, i.e. remove NaNs
        i_q_vec = i_q_vec[NaN_mask]
        lambda_d_vec = lambda_d_vec[NaN_mask]
        # these are the data points inside the convex hull
        points = np.vstack((i_d_vec, i_q_vec)).T
        # compute the convex hull
        hull = Delaunay(points)
        # rearrange the REF trajectory points
        REF_points = np.vstack((self.i_d_REF.flatten(), self.i_q_REF.flatten())).T
        # check which points of the REF trajectory are located inside the convex hull
        in_hull = hull.find_simplex(REF_points)>=0
        # keep only the points of the REF trajectory that are located inside the convex hull
        self.i_d_REF = self.i_d_REF[in_hull]
        self.i_q_REF = self.i_q_REF[in_hull]
        self.i_REF = self.i_REF[in_hull]
        self.theta_REF = self.theta_REF[in_hull]

        
        # evaluate epsilon along the reference trajectory
        try:
            epsilon_REF = interpn((self.i_d[:,0], self.i_q[0,:]), self.epsilon, (self.i_d_REF, self.i_q_REF))
        except ValueError:
            # if one of the requested xi is out of bounds in dimension 1:
            # remove the point of the REF trajectory with bigger amplitude
            self.i_REF = self.i_REF[0:-2]
            self.theta_REF = self.theta_REF[0:-2]
            self.i_d_REF = self.i_d_REF[0:-2]
            self.i_q_REF = self.i_q_REF[0:-2]
            # compute again the interpolation
            epsilon_REF = interpn((self.i_d[:,0], self.i_q[0,:]), self.epsilon, (self.i_d_REF, self.i_q_REF))
        
        # store the result
        self.epsilon_REF = epsilon_REF
        
        # plt.plot(self.i_REF, epsilon_REF, label="epsilon_REF")
        # plt.plot(self.i_REF, self.theta_REF, label="theta_REF")
        # plt.plot(self.i_REF, self.theta_REF + epsilon_REF, label="theta_REF + epsilon_REF = t1")
        # plt.legend()
        # plt.show()
        
        # add epsilon to the REF angle -> obtain the angle of sensored trajectory t1
        theta_sensored = self.theta_REF + epsilon_REF
        
        # set output (sensored trajectory) in cartesian coordinates
        self.i_d_sensored = self.i_REF * np.cos(theta_sensored)
        self.i_q_sensored = self.i_REF * np.sin(theta_sensored)
        
        # FIX: Cut the trajectory t1 if it oversteps the convex hull
        # rearrange the REF trajectory points
        t1_points = np.vstack((self.i_d_sensored.flatten(), self.i_q_sensored.flatten())).T
        # check which points of the trajectory t1 are located inside the convex hull
        in_hull = hull.find_simplex(t1_points)>=0
        # keep only the points of the trajectory t1 that are located inside the convex hull
        self.i_d_sensored = self.i_d_sensored[in_hull]
        self.i_q_sensored = self.i_q_sensored[in_hull]

        
        print('computed sensored trajectory')
        

    def calc_convergence_region(self, Uh, fh, Ihq_star=0):
        """
        Compute the signal Ihq, input of the position observer,
        considering a given reference trajectory as current reference.
        The convergence region is the locus Ihq=0 when, for increasing angles,
        the slope is negative (Ihq_neg)
        Input:
        Uh: amplitude of the HF voltage injection (V)
        fh: amplitude of the HF voltage injection (Hz)
        """
        # compute i_d_REF and i_q_REF if they are not defined
        if (not hasattr(self, 'i_d_REF')) or (not hasattr(self, 'i_q_REF')):
            try:
                self.i_d_REF = self.i_REF * np.cos(self.theta_REF)
                self.i_q_REF = self.i_REF * np.sin(self.theta_REF)
            except :
                print('  APOLLO ERROR in calc_convergence_region:')
                print('  reference trajectory seems not assigned')
                print('  consider to assign i_REF and theta_REF values (or i_d_REF and i_q_REF) to fm object first')
                exit()
                return
                    
        # compute i_REF and theta_REF if they are not defined
        if (not hasattr(self, 'i_REF')) or (not hasattr(self, 'theta_REF')):
            try:
                self.i_REF = np.sqrt(self.i_d_REF**2 + self.i_q_REF**2)
                self.theta_REF = np.arctan2(self.i_q_REF, self.i_d_REF)
            except :
                print('  APOLLO ERROR in calc_convergence_region:')
                print('  reference trajectory seems not assigned')
                print('  consider to assign i_REF and theta_REF values (or i_d_REF and i_q_REF) to fm object first')
                exit()
                return
                
                

        # Obtain coordinate matrices in polar form
        I = np.sqrt(self.i_d**2 + self.i_q**2)
        alpha_ie = np.arctan2(self.i_q, self.i_d)
        alpha_ie = alpha_ie % (2 * np.pi)  # modulus at 2 pi

        # compute Ihq amplitude
        amplitude = Uh / (2 * np.pi * fh) * np.divide(np.sqrt(self.ldelta**2 + self.ldq**2), self.ldd * self.lqq - self.ldq * self.lqd)
        
        # compute delta_theta
        self.delta_theta = alpha_ie - np.interp(I, self.i_REF, self.theta_REF)
        
        # compute Ihq
        self.Ihq = - amplitude * np.sin(2 * self.delta_theta - 2 * self.epsilon)

        # calc Ihq_d_alpha_ie, tangential derivative (with respect to angle) of Ihq
        Ihq_dx = np.gradient(self.Ihq, self.i_d[:,0], axis=0)
        Ihq_dy = np.gradient(self.Ihq, self.i_q[0,:], axis=1)
        Ihq_d_alpha_ie = - Ihq_dx * self.i_q + Ihq_dy * self.i_d
        self.Ihq_d_alpha_ie = Ihq_d_alpha_ie

        # extract Ihq locus with negative slope
        self.Ihq_neg = copy.deepcopy(self.Ihq)
        self.Ihq_neg[Ihq_d_alpha_ie > 0] = np.nan


        # add Ihq_star compensation

        # compute Ihq_in
        self.Ihq_in = self.Ihq + Ihq_star

        # calc Ihq_in_d_alpha_ie, tangential derivative (with respect to angle) of Ihq_in
        Ihq_in_dx = np.gradient(self.Ihq_in, self.i_d[:,0], axis=0)
        Ihq_in_dy = np.gradient(self.Ihq_in, self.i_q[0,:], axis=1)
        Ihq_in_d_alpha_ie = - Ihq_in_dx * self.i_q + Ihq_in_dy * self.i_d
        self.Ihq_in_d_alpha_ie = Ihq_in_d_alpha_ie

        # extract Ihq_in locus with negative slope
        self.Ihq_in_neg = copy.deepcopy(self.Ihq_in)
        self.Ihq_in_neg[Ihq_in_d_alpha_ie > 0] = np.nan
        
        print('computed convergence region with respect to trajectory REF')
        
        self.Uh=Uh
        self.fh=fh


    def plot_ellipses(self, x, y, Uh, fh, k=1):
        """
        Superimpose a plot of the hf ellipses for various points in the dq current plane
        The ellipses represent the current response to a rotating voltage injection
        Input:
        x: i_d query points (A)
        y: i_q query points (A)
        Uh: amplitude of the HF voltage injection (V)
        fh: amplitude of the HF voltage injection (Hz)
        k: scaling factor
        """
        # clip query bounds
        x = np.clip(x, np.amin(self.i_d[:,0]), np.amax(self.i_d[:,0]))
        y = np.clip(y, np.amin(self.i_q[0,:]), np.amax(self.i_q[0,:]))
        
        # flatten query breakpoints
        x, y = np.meshgrid(x, y)
        x = x.flatten()
        y = y.flatten()      
            
        # interpolate on the query grid points
        lsigma = interpn((self.i_d[:,0], self.i_q[0,:]), self.lsigma, (x, y), method="nearest")
        ldelta = interpn((self.i_d[:,0], self.i_q[0,:]), self.ldelta, (x, y), method="nearest")
        ldq = interpn((self.i_d[:,0], self.i_q[0,:]), self.ldq, (x, y), method="nearest")
        epsilon_deg = interpn((self.i_d[:,0], self.i_q[0,:]), self.epsilon_deg, (x, y), method="nearest")

        # compute major and minor semi-axes (and scale them with k)
        SM = Uh / (2 * np.pi * fh) * np.divide((lsigma + np.sqrt(ldelta**2 + ldq**2)), (lsigma**2 - ldelta**2 - ldq**2)) * k
        sm = Uh / (2 * np.pi * fh) * np.divide((lsigma - np.sqrt(ldelta**2 + ldq**2)), (lsigma**2 - ldelta**2 - ldq**2)) * k
        
        # define ellipses
        ells = [Ellipse(xy=(x[i], y[i]), width=SM[i], height=sm[i], angle=epsilon_deg[i], edgecolor='k', fc='None', lw=1) for i in range(x.size)]
        
        # plot ellipses
        fig = plt.figure(0)
        ax = fig.add_subplot(111, aspect='equal')
        for e in ells:
            ax.add_artist(e)
        ax.set_xlim(np.amin(x) - 2*np.amax(SM[~np.isnan(SM)]), np.amax(x) + 2*np.amax(SM[~np.isnan(SM)]))
        ax.set_ylim(np.amin(y) - 2*np.amax(SM[~np.isnan(SM)]), np.amax(y) + 2*np.amax(SM[~np.isnan(SM)]))
        
        # return the HF ellipses plot
        plot_out = plt.gcf()
        return plot_out
        
        # close plt
        plt.close()
                    

    def fit_flux_SynRM(self, method="average", S=None, T=None, U=None, V=None,
                        Id_min=None, Id_max=None, n_Id=None, Iq_min=None, Iq_max=None, n_Iq=None):
        """
        Fit the flux linkages using the functions in [4].
        The fitting has a smoothing effect on the maps.
        NOTE: the fitting procedure can be used for RELUCTANCE MOTOR ONLY
        (no permanent magnet).
        
        Two methods are available:
        "average": Two separate fitting are computed, on d and q-axis.
                   The average between the common parameters adq, U, V
                   is computed.
        "following": The fitting on axis d is done at first. The parameters
                   adq, U, V are stored and used as bounds for the following
                   fitting on axis q.
                   
        You can set manually the coefficients S, T, U and V. For example:
        mot.fit_reluctance(method = "average", S=6, T=1, U=2, V=0)
        
        You can also set the output Id-Iq grid. Otherwise, the default output grid
        is the same of the input.
        """
        
        # Define fitting functions.
        # func1 is for axis d; eq (12a) in [4]
        # func2 is for axis q; eq (12b) in [4]
        def func1(x, ad0, add, adq, S, U, V):
            return (ad0 + add * np.power(np.absolute(x[0]), S) + adq / (V + 2) *
                 np.power(np.absolute(x[0]), U) * np.power(np.absolute(x[1]), V + 2)) * x[0]

        def func2(x, aq0, aqq, adq, T, U, V):
            return (aq0 + aqq * np.power(np.absolute(x[1]), T) + adq / (U + 2) *
                 np.power(np.absolute(x[0]), U + 2) * np.power(np.absolute(x[1]), V)) * x[1]
        
                 
        # prepare arrays for curve_fit (flatten and remove NaNs)
        # flatten
        i_d_vec = self.i_d.flatten()
        i_q_vec = self.i_q.flatten()
        lambda_d_vec = self.lambda_d.flatten()
        lambda_q_vec = self.lambda_q.flatten()
        # remove NaNs
        i_d_vec = i_d_vec[~np.isnan(lambda_q_vec)]
        i_q_vec = i_q_vec[~np.isnan(lambda_q_vec)]
        lambda_d_vec = lambda_d_vec[~np.isnan(lambda_q_vec)]
        lambda_q_vec = lambda_q_vec[~np.isnan(lambda_q_vec)]   
        
        # define the xdata
        xdata = np.vstack((lambda_d_vec, lambda_q_vec))
        
        # define the y data
        ydata1 = i_d_vec
        
        # Define the bounds
        if S is None:
            S_min = 0
            S_max = 20
        else:
            S_min = np.clip(S - 1e-6, 0, 20)
            S_max = np.clip(S + 1e-6, 0, 20)
        
        if U is None:
            U_min = 0
            U_max = 20
        else:
            U_min = np.clip(U - 1e-6, 0, 20)
            U_max = np.clip(U + 1e-6, 0, 20)
            
        if V is None:
            V_min = 0
            V_max = 20
        else:
            V_min = np.clip(V - 1e-6, 0, 20)
            V_max = np.clip(V + 1e-6, 0, 20)

        bounds1 = ((0, 0, 0, S_min, U_min, V_min), (np.Inf, np.Inf, np.Inf, S_max, U_max, V_max))
        
        # Compute the first fitting on axis d
        popt1, pcov1 = curve_fit(func1, xdata, ydata1, bounds = bounds1, maxfev=2000)
        #yfit1 = func1(xdata, *popt1)
        
        #print(popt1)
        ad0 = popt1[0]
        add = popt1[1]
        adq1 = popt1[2]
        S = popt1[3]
        U1 = popt1[4]
        V1 = popt1[5]
        
        if method == "average":
            # interpolate separately on both d- and q-axis
            # and average the coefficients
        
            # define the ydata
            ydata2 = i_q_vec
            
            # Define the bounds
            if T is None:
                T_min = 0
                T_max = 20
            else:
                T_min = np.clip(T - 1e-6, 0, 20)
                T_max = np.clip(T + 1e-6, 0, 20)
                
            if U is None:
                U_min = 0
                U_max = 20
            else:
                U_min = np.clip(U - 1e-6, 0, 20)
                U_max = np.clip(U + 1e-6, 0, 20)
                
            if V is None:
                V_min = 0
                V_max = 20
            else:
                V_min = np.clip(V - 1e-6, 0, 20)
                V_max = np.clip(V + 1e-6, 0, 20)
            
            bounds2 = ((0, 0, 0, T_min, U_min, V_min), (np.Inf, np.Inf, np.Inf, T_max, U_max, V_max))
            
            # Compute the second fitting on axis q
            popt2, pcov2 = curve_fit(func2, xdata, ydata2, bounds = bounds2, maxfev=2000)
            #yfit2 = func2(xdata, *popt2)

            #print(popt2)
            aq0 = popt2[0]
            aqq = popt2[1]
            adq2 = popt2[2]
            T = popt2[3]
            U2 = popt2[4]
            V2 = popt2[5]

            # average coefficients
            adq = (adq1+adq2)/2
            U = (U1+U2)/2
            V = (V1+V2)/2

        if method == "following":
            # fit the d-axis flux linkage and keep those values for the q-axis interpolation
            
            # define the ydata
            ydata2 = i_q_vec
            
            # Define the bounds
            adq_min = adq1 - 1e-6
            adq_max = adq1 + 1e-6
            if T is None:
                T_min = 0
                T_max = 20
            else:
                T_min = np.clip(T - 1e-6, 0, 20)
                T_max = np.clip(T + 1e-6, 0, 20)
            U_min = U1 - 1e-6
            U_max = U1 + 1e-6
            V_min = V1 - 1e-6
            V_max = V1 + 1e-6
            bounds2 = ((0, 0, adq_min, T_min, U_min, V_min), (np.Inf, np.Inf, adq_max, T_max, U_max, V_max))
            
            # Compute the second fitting on axis q
            popt2, pcov2 = curve_fit(func2, xdata, ydata2, bounds = bounds2, maxfev=2000)         
            #yfit2 = func2(xdata, *popt2)

            #print(popt2)
            aq0 = popt2[0]
            aqq = popt2[1]
            adq = popt2[2]
            T = popt2[3]
            U = popt2[4]
            V = popt2[5]

        
        # print resulting fitting parameters
        print('fitted maps using',method,':')
        print("ad0=", '%.2f'%ad0, "  add=", '%.2f'%add, "  adq=", '%.2f'%adq,
        "  aq0=", '%.2f'%aq0, "  aqq=", '%.2f'%aqq, "      S=", '%.2f'%S,
        "  T=", '%.2f'%T, "  U=", '%.2f'%U, "  V=", '%.2f'%V)
        
        # We obtained the fitting coefficients ad0,add,adq,aq0,aqq,S,T,U,V
        # With these coefficients we can express the currents
        # i_d,i_q as function of the fluxes lambda_d,lambda_q
        
        # The fitted currents will be computed on a FluxD-FluxQ grid.
        # Automatically set the FluxD-FluxQ grid properties:
        FluxD_min = -3 * np.max(np.abs(lambda_d_vec))
        FluxD_max =  3 * np.max(np.abs(lambda_d_vec))
        n_FluxD = self.i_d[:,0]
        if self.i_d.shape[0] < 251:
            n_FluxD = 251
        FluxQ_min = -3 * np.max(np.abs(lambda_q_vec))
        FluxQ_max =  3 * np.max(np.abs(lambda_q_vec))
        n_FluxQ = self.i_q[0,:]
        if self.i_q.shape[1] < 251:
            n_FluxQ = 251
            
        # Create the provisional FluxD-FluxQ grid
        FluxD_breakpoints_fg = np.linspace(FluxD_min, FluxD_max, n_FluxD)
        FluxQ_breakpoints_fg = np.linspace(FluxQ_min, FluxQ_max, n_FluxQ)
        FluxD_fg, FluxQ_fg = np.meshgrid(FluxD_breakpoints_fg, FluxQ_breakpoints_fg)
        FluxD_fg = FluxD_fg.transpose()
        FluxQ_fg = FluxQ_fg.transpose()

        # Compute the fitted currents on that grid
        Id_fg = (ad0 + add * np.power(np.absolute(FluxD_fg), S) + adq / (V + 2) *
                 np.power(np.absolute(FluxD_fg), U) * np.power(np.absolute(FluxQ_fg), V + 2)) * FluxD_fg
        Iq_fg = (aq0 + aqq * np.power(np.absolute(FluxQ_fg), T) + adq / (U + 2) *
                 np.power(np.absolute(FluxD_fg), U + 2) * np.power(np.absolute(FluxQ_fg), V)) * FluxQ_fg
                 
        if np.any(np.isnan(Id_fg)) or np.any(np.isnan(Iq_fg)):
            print('  APOLLO ERROR in fit_reluctance:')
            print('  impossible to fit the flux maps')
            if method == 'average':
                print('  please use "following" method')
            if method == 'following':
                print('  please use "average" method')
            exit()
        
        # We computed the fitted currents on a flux grid.
        # Our aim is to invert the current maps (flux grid):
        # i_d(lambda_d,lambda_q), i_q(lambda_d,lambda_q)
        # into flux maps (current grid):
        # lambda_d(i_d,i_q), lambda_q(i_d,i_q)
        
        # set the properties of the output Id-Iq grid
        if Id_min is None:
            Id_min = np.min(self.i_d[:,0])
        if Id_max is None:
            Id_max = np.max(self.i_d[:,0])
        if n_Id is None:
            n_Id = self.i_d.shape[0]
        if Iq_min is None:
            Iq_min = np.min(self.i_q[0,:])
        if Iq_max is None:
            Iq_max = np.max(self.i_q[0,:])
        if n_Iq is None:
            n_Iq = self.i_q.shape[1]
        
        # create the output Id-Iq grid
        Id_breakpoints_cg = np.linspace(Id_min, Id_max, n_Id)
        Iq_breakpoints_cg = np.linspace(Iq_min, Iq_max, n_Iq)
        Id_cg, Iq_cg = np.meshgrid(Id_breakpoints_cg, Iq_breakpoints_cg)
        Id_cg = Id_cg.transpose()
        Iq_cg = Iq_cg.transpose()

        # Invert the current maps to flux maps
        points_I = np.vstack((Id_fg.flatten(), Iq_fg.flatten())).T
        values_lambda_d = FluxD_fg.flatten()
        values_lambda_q = FluxQ_fg.flatten()
        lambda_d_fit = griddata(points_I, values_lambda_d, (Id_cg, Iq_cg), method='cubic')
        lambda_q_fit = griddata(points_I, values_lambda_q, (Id_cg, Iq_cg), method='cubic')
        
        # interpolate Maxwell torque on the new grid
        points_I_torque = np.vstack((self.i_d.flatten(), self.i_q.flatten())).T
        values_torque = self.torque.flatten()
        torque_fit = griddata(points_I_torque, values_torque, (Id_cg, Iq_cg), method='cubic')
        
        # save fitting results
        self.i_d = Id_cg
        self.i_q = Iq_cg
        self.lambda_d = lambda_d_fit
        self.lambda_q = lambda_q_fit
        self.torque = torque_fit
        
        
        
    def calc_fourier_inverse_inductances(self, I='', k=2, method='fitting', steps=1000, plot=False):
        """
        This function calculates the Fourier coefficients of the inverse incremental
        inductances gamma_delta and gamma_dq in polar coordinates for a current modulus I.
        It also reconstructed Ihq signal with the Fourier series of the inverse incremental inductances
        
        You set manually the current modulus value I in which this function calculates
        the fourier coefficients and the Ihq, the harmonic order k and the method by which it
        calculates the coefficients.
        
        Two methods are available:
        
        "fitting": fit the inverse incremental inductances using fist terms of
                   the Fourier series as a function, default
        
        "fft": compute fast Fourier trasform on the inverse incremental inductances 
               to obtain the coefficients of the Fourier series
        
        
        steps: the number of breakpoints to create the angle vector, default is 1000
        
        this function stores the results in the following variables:
        
        self.gamma_delta_fourier_coef: the computed Fourier series coefficients (k+1 elements)
        self.gamma_delta_fourier_rec: the reconstructed gamma_delta using the Fourier series
        self.gamma_dq_fourier_coef: the computed Fourier series coefficients (k+1 elements)
        self.gamma_dq_fourier_rec: the reconstructed gamma_dq using the Fourier series
        self.alpha_breakpoints: angles for polar coordinates
        self.Ihq_fourier: the reconstructed Ihq using the Fourier series of inverse incremental inductances

        """
        # trasform into polar coordinates
        alpha_array = np.arctan2(self.i_q, self.i_d)
        alpha_array = np.mod(alpha_array, 2*np.pi) 
        I_array = np.sqrt(self.i_d**2 + self.i_q**2)

        # Define self.alpha_breakpoints
        self.alpha_breakpoints = math.pi/180*np.linspace(0,360,steps)
        
        # convert two-dimensional vector into one-dimensional 
        alpha_array = alpha_array.flatten()
        I_array = I_array.flatten()
        gamma_delta = self.gamma_delta.flatten()
        gamma_dq = self.gamma_dq.flatten()
        
        # remove NaN
        NaN_mask = ~np.isnan(gamma_delta) 
        I_array = I_array[NaN_mask] 
        alpha_array = alpha_array[NaN_mask] 
        gamma_delta = gamma_delta[NaN_mask]
        gamma_dq = gamma_dq[NaN_mask]
        
        
        points_I = np.vstack((alpha_array, I_array)).T
        
        # check the consistency of input data
        if method != 'fft' and method != 'fitting' :
            print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
            print('  impossible to compute the fourier coefficients of the inverse incremental inductance')
            print('  please use "fft" or "fitting" method')
            exit()
            
        if np.min(alpha_array) >= 0.2 or np.max(alpha_array) <= 6.1 :
            print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
            print('  impossible to compute the fourier coefficients of inverse incremental inductance')
            print('  Necessary to have the measurements of the inverse inductances along the entire period from 0 to 2pi')
            exit()
            
    
        if I=='' :
            print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
            print('  impossible to compute the fourier coefficients of inverse incremental inductance')
            print('  please set the current modulus value')
            exit()
            
        elif I < 0 or I > np.min([np.max(np.max(self.i_d)),np.max(np.max(self.i_q))]) :
            print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
            print('  impossible to compute the fourier coefficients of inverse incremental inductance')
            print('  please set I values between 0 and i_d_max or i_q_max')
            exit()
            
        # compute griddata for the current modulus value I
        gamma_delta_polar_I = griddata(points_I, gamma_delta, (self.alpha_breakpoints, [I] * steps), method='cubic')
        gamma_dq_polar_I    = griddata(points_I, gamma_dq, (self.alpha_breakpoints, [I] * steps), method='cubic')
        print('I=',I)
    
        if method == 'fft': #compute fast Fourier transform (fft)
        
            # coefficients of Fourier series for gamma_delta
            N = len(gamma_delta_polar_I)
            # replace Nan with values obtained through piecewise cubic spline interpolation
            gamma_delta_polar_I = np.concatenate((gamma_delta_polar_I[math.ceil(N/2):],gamma_delta_polar_I[:math.ceil(N/2)]),axis=None)
            gamma_delta_polar_I = pd.DataFrame(gamma_delta_polar_I).interpolate(method='cubic').values.ravel()
            gamma_delta_polar_I = np.concatenate((gamma_delta_polar_I[math.ceil(N/2):],gamma_delta_polar_I[:math.ceil(N/2)]),axis=None)
            gamma_delta_polar_I_mask = np.isnan(gamma_delta_polar_I)
            if np.any(gamma_delta_polar_I_mask) == True :
                print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                print('  impossible to compute the fourier coefficients of the delta inverse incremental inductance with fft method')
                print('  NaN found during the fft computation; try to use a higher current value I')
                exit()
            else:
                if k < 2 or k >((N-1)/2):
                    print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                    print('  impossible to compute the fourier coefficients of the delta inverse incremental inductance with fft')
                    print('  please use a different value of k')
                    exit()
                else:
                    Z = fft(gamma_delta_polar_I)
                    Z = Z/N
                    #create vector that will contain the first k coefficients of the Fourier series: 
                    #self.gamma_delta_fourier_coef = [a0 a1 b1 a2 b2 a3 b3 ...]
                    #a_k for cosines and b_k for sines where k is the harmonic order.
                    self.gamma_delta_fourier_coef = np.zeros (2*k+1)
                    self.gamma_delta_fourier_coef[0] = np.real(Z[0])
                    self.gamma_delta_fourier_rec = np.zeros ( ((k+1),(len(self.alpha_breakpoints))) )
                    self.gamma_delta_fourier_rec[0,:]=self.gamma_delta_fourier_coef[0] 
                    print('Gamma_delta')
                    print('a',0,'=','%.4f'%self.gamma_delta_fourier_coef[0])
                    for i in range (1,k+1):
                        self.gamma_delta_fourier_coef[2*i-1] = self.gamma_delta_fourier_coef[2*i-1] + 2 * np.real(Z[i])
                        self.gamma_delta_fourier_coef[2*i] = self.gamma_delta_fourier_coef[2*i] - 2 * np.imag(Z[i])
                        print('a',i,'=','%.4f'%self.gamma_delta_fourier_coef[2*i-1],'b',i,'=','%.4f'%self.gamma_delta_fourier_coef[2*i])
                        self.gamma_delta_fourier_rec[i,:] = self.gamma_delta_fourier_coef[2*i-1]*np.cos(i*self.alpha_breakpoints) + self.gamma_delta_fourier_coef[2*i]*np.sin(i*self.alpha_breakpoints)
                    self.gamma_delta_fourier_rec = np.sum(self.gamma_delta_fourier_rec[:,],axis=0)
            
            # coefficients of Fourier series for gamma_dq
            N = len(gamma_dq_polar_I)
            # replace Nan with values obtained through piecewise cubic spline interpolation
            gamma_dq_polar_I = np.concatenate((gamma_dq_polar_I[math.ceil(N/2):],gamma_dq_polar_I[:math.ceil(N/2)]),axis=None)
            gamma_dq_polar_I = pd.DataFrame(gamma_dq_polar_I).interpolate(method='cubic').values.ravel()
            gamma_dq_polar_I = np.concatenate((gamma_dq_polar_I[math.ceil(N/2):],gamma_dq_polar_I[:math.ceil(N/2)]),axis=None)
            gamma_dq_polar_I_mask = np.isnan(gamma_dq_polar_I)
            if np.any(gamma_dq_polar_I_mask) == True :
                print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                print('  impossible to compute the fourier coefficients of the dq inverse incremental inductance with fft method')
                print('  NaN found during the fft computation; try to use a higher current value I') 
                exit()    
            else:
                if k<2 and k>((N-1)/2):
                    print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                    print('  impossible to compute the fourier coefficients of the dq inverse incremental inductance with fft')
                    print('  please use differet value of k')
                    exit()           
                else:
                    Z = fft(gamma_dq_polar_I)
                    Z = Z/N
                    #create vector that will contain the first k coefficients of the Fourier series: 
                    #self.gamma_dq_fourier_coef = [a0 a1 b1 a2 b2 a3 b3 ...]
                    #a_k for cosines and b_k for sines where k is the harmonic order
                    self.gamma_dq_fourier_coef = np.zeros (2*k+1)
                    self.gamma_dq_fourier_coef[0] = np.real(Z[0])
                    self.gamma_dq_fourier_rec = np.zeros ( ((k+1),(len(self.alpha_breakpoints))) )
                    self.gamma_dq_fourier_rec[0,:]=self.gamma_dq_fourier_coef[0]
                    print('Gamma_dq')
                    print('a',0,'=','%.4f'%self.gamma_dq_fourier_coef[0])
                    for i in range (1,k+1):
                        self.gamma_dq_fourier_coef[2*i-1] = self.gamma_dq_fourier_coef[2*i-1] + 2 * np.real(Z[i])
                        self.gamma_dq_fourier_coef[2*i] = self.gamma_dq_fourier_coef[2*i] - 2 * np.imag(Z[i])
                        print('a',i,'=','%.4f'%self.gamma_dq_fourier_coef[2*i-1],'b',i,'=','%.4f'%self.gamma_dq_fourier_coef[2*i])
                        self.gamma_dq_fourier_rec[i,:] = self.gamma_dq_fourier_coef[2*i-1]*np.cos(i*self.alpha_breakpoints) + self.gamma_dq_fourier_coef[2*i]*np.sin(i*self.alpha_breakpoints)
                    self.gamma_dq_fourier_rec = np.sum(self.gamma_dq_fourier_rec[:,],axis=0)
                    
        #Fit the inverse incremental inductances using first terms of the Fourier series as a function  
        elif method == 'fitting':      
        
            # coefficients of Fourier series for gamma_delta
            #Remove NaN for fitting
            NaN_mask = ~np.isnan(gamma_delta_polar_I) 
            alpha_breakpoints_delta = self.alpha_breakpoints[NaN_mask] 
            gamma_delta_polar_I_nn = gamma_delta_polar_I[NaN_mask]
            if len(alpha_breakpoints_delta) == 0 :
                print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                print('  impossible to compute the fourier coefficients of the delta inverse incremental inductance with fitting method')
                print('  NaN found during the fitting computation; try to use a higher current value I')
                exit()
            else:
                if k == 2 :
                    def fun(x, a0, a1, b1, a2, b2):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x)
                    self.gamma_delta_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_delta, gamma_delta_polar_I_nn,maxfev=2000)
                    print('Gamma_delta',"a0=", '%.4f'%self.gamma_delta_fourier_coef[0], "a1=", '%.4f'%self.gamma_delta_fourier_coef[1],"b1=", '%.4f'%self.gamma_delta_fourier_coef[2], \
                    "a2=", '%.4f'%self.gamma_delta_fourier_coef[3],"b2=", '%.4f'%self.gamma_delta_fourier_coef[4])
                    self.gamma_delta_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_delta_fourier_coef)
                    
                elif k == 3 :
                    def fun(x, a0, a1, b1, a2, b2, a3, b3):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x) + a3*np.cos(3*x) + b3*np.sin(3*x)
                    self.gamma_delta_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_delta, gamma_delta_polar_I_nn,maxfev=2000)
                    print('Gamma_delta',"a0=", '%.4f'%self.gamma_delta_fourier_coef[0], "a1=", '%.4f'%self.gamma_delta_fourier_coef[1],"b1=", '%.4f'%self.gamma_delta_fourier_coef[2], \
                            "a2=", '%.4f'%self.gamma_delta_fourier_coef[3],"b2=", '%.4f'%self.gamma_delta_fourier_coef[4],"a3=", '%.4f'%self.gamma_delta_fourier_coef[5],"b3=", '%.4f'%self.gamma_delta_fourier_coef[6])
                    self.gamma_delta_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_delta_fourier_coef)
                    
                elif k ==4 :
                    def fun(x, a0, a1, b1, a2, b2, a3, b3, a4, b4):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x) + a3*np.cos(3*x) + b3*np.sin(3*x)+ a4*np.cos(4*x) + b4*np.sin(4*x)
                    self.gamma_delta_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_delta, gamma_delta_polar_I_nn,maxfev=2000)
                    print('Gamma_delta',"a0=", '%.4f'%self.gamma_delta_fourier_coef[0], "a1=", '%.4f'%self.gamma_delta_fourier_coef[1],"b1=", '%.4f'%self.gamma_delta_fourier_coef[2], \
                    "a2=", '%.4f'%self.gamma_delta_fourier_coef[3],"b2=", '%.4f'%self.gamma_delta_fourier_coef[4],"a3=", '%.4f'%self.gamma_delta_fourier_coef[5],"b3=", '%.4f'%self.gamma_delta_fourier_coef[6],\
                    "a4=", '%.4f'%self.gamma_delta_fourier_coef[7],"b4=", '%.4f'%self.gamma_delta_fourier_coef[8])
                    self.gamma_delta_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_delta_fourier_coef)
                    
                elif k==5:
                    def fun(x, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x) + a3*np.cos(3*x) + b3*np.sin(3*x)+ a4*np.cos(4*x) + b4*np.sin(4*x) + a5*np.cos(5*x) + b5*np.sin(5*x)
                    self.gamma_delta_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_delta, gamma_delta_polar_I_nn,maxfev=2000)
                    print('Gamma_delta',"a0=", '%.4f'%self.gamma_delta_fourier_coef[0], "a1=", '%.4f'%self.gamma_delta_fourier_coef[1],"b1=", '%.4f'%self.gamma_delta_fourier_coef[2],\
                    "a2=", '%.4f'%self.gamma_delta_fourier_coef[3],"b2=", '%.4f'%self.gamma_delta_fourier_coef[4],"a3=", '%.4f'%self.gamma_delta_fourier_coef[5],"b3=", '%.4f'%self.gamma_delta_fourier_coef[6],\
                    "a4=", '%.4f'%self.gamma_delta_fourier_coef[7],"b4=", '%.4f'%self.gamma_delta_fourier_coef[8],"a5=", '%.4f'%self.gamma_delta_fourier_coef[9],"b5=", '%.4f'%self.gamma_delta_fourier_coef[10])
                    self.gamma_delta_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_delta_fourier_coef)
                    
                else :
                    print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                    print('  impossible to compute the fourier coefficients of the delta inverse incremental inductance with fitting method')
                    print('  please use k between 2 and 5')
                    exit()
                    
            # coefficients of Fourier series for gamma_dq
            #Remove NaN for fitting
            NaN_mask = ~np.isnan(gamma_dq_polar_I) 
            alpha_breakpoints_dq = self.alpha_breakpoints[NaN_mask]
            gamma_dq_polar_I_nn = gamma_dq_polar_I[NaN_mask]
            if len(alpha_breakpoints_dq) == 0 :
                print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                print('  impossible to compute the fourier coefficients of the dq inverse incremental inductance with fitting method')
                print('  NaN found during the fitting computation; try to use a higher current value I')
                exit()
            else:    
                if k == 2 :
                    def fun(x, a0, a1, b1, a2, b2):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x)
                    self.gamma_dq_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_dq, gamma_dq_polar_I_nn,maxfev=2000)
                    print('Gamma_dq',"a0=", '%.4f'%self.gamma_dq_fourier_coef[0], "a1=", '%.4f'%self.gamma_dq_fourier_coef[1],"b1=", '%.4f'%self.gamma_dq_fourier_coef[2],\
                    "a2=", '%.4f'%self.gamma_dq_fourier_coef[3],"b2=", '%.4f'%self.gamma_dq_fourier_coef[4])
                    self.gamma_dq_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_dq_fourier_coef)
                    
                elif k == 3 :
                    def fun(x, a0, a1, b1, a2, b2, a3, b3):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x) + a3*np.cos(3*x) + b3*np.sin(3*x)
                    self.gamma_dq_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_dq, gamma_dq_polar_I_nn,maxfev=2000)
                    print('Gamma_dq',"a0=", '%.4f'%self.gamma_dq_fourier_coef[0], "a1=", '%.4f'%self.gamma_dq_fourier_coef[1],"b1=", '%.4f'%self.gamma_dq_fourier_coef[2], \
                    "a2=", '%.4f'%self.gamma_dq_fourier_coef[3],"b2=", '%.4f'%self.gamma_dq_fourier_coef[4],"a3=", '%.4f'%self.gamma_dq_fourier_coef[5],"b3=", '%.4f'%self.gamma_dq_fourier_coef[6])
                    self.gamma_dq_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_dq_fourier_coef)
                    
                elif k == 4 :
                    def fun(x, a0, a1, b1, a2, b2, a3, b3, a4, b4):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x) + a3*np.cos(3*x) + b3*np.sin(3*x)+ a4*np.cos(4*x) + b4*np.sin(4*x)
                    self.gamma_dq_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_dq, gamma_dq_polar_I_nn,maxfev=2000)
                    print('Gamma_dq',"a0=", '%.4f'%self.gamma_dq_fourier_coef[0], "a1=", '%.4f'%self.gamma_dq_fourier_coef[1],"b1=", '%.4f'%self.gamma_dq_fourier_coef[2], \
                    "a2=", '%.4f'%self.gamma_dq_fourier_coef[3],"b2=", '%.4f'%self.gamma_dq_fourier_coef[4],"a3=", '%.4f'%self.gamma_dq_fourier_coef[5],"b3=", '%.4f'%self.gamma_dq_fourier_coef[6],\
                    "a4=", '%.4f'%self.gamma_dq_fourier_coef[7],"b4=", '%.4f'%self.gamma_dq_fourier_coef[8])
                    self.gamma_dq_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_dq_fourier_coef)
                    
                elif k == 5 :
                    def fun(x, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5):
                        return a0 + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x) + a3*np.cos(3*x) + b3*np.sin(3*x)+ a4*np.cos(4*x) + b4*np.sin(4*x) + a5*np.cos(5*x) + b5*np.sin(5*x)
                    self.gamma_dq_fourier_coef, pcov1 = curve_fit(fun, alpha_breakpoints_dq, gamma_dq_polar_I_nn,maxfev=2000)
                    print('Gamma_dq',"a0=", '%.4f'%self.gamma_dq_fourier_coef[0], "a1=", '%.4f'%self.gamma_dq_fourier_coef[1],"b1=", '%.4f'%self.gamma_dq_fourier_coef[2], \
                    "a2=", '%.4f'%self.gamma_dq_fourier_coef[3],"b2=", '%.4f'%self.gamma_dq_fourier_coef[4],"a3=", '%.4f'%self.gamma_dq_fourier_coef[5],"b3=", '%.4f'%self.gamma_dq_fourier_coef[6],\
                    "a4=", '%.4f'%self.gamma_dq_fourier_coef[7],"b4=", '%.4f'%self.gamma_dq_fourier_coef[8],"a5=", '%.4f'%self.gamma_dq_fourier_coef[9],"b5=", '%.4f'%self.gamma_dq_fourier_coef[10])
                    self.gamma_dq_fourier_rec = fun(self.alpha_breakpoints, *self.gamma_dq_fourier_coef)
                  
                    
                else :
                    print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
                    print('  impossible to compute the fourier coefficients of the dq inverse incremental inductance with fitting method')
                    print('  please use k between 2 and 5')
                    exit()
            
        # plot
        # if you need to plot outside, it is necessary to store alpha_breakpoints_delta
        # and gamma_delta_polar_I (there are some Nan to take into account)
        if plot:
            fig, ax = plt.subplots()
            ax.plot(self.alpha_breakpoints, gamma_delta_polar_I, lw=2, label="$\gamma_\Delta$ input")
            ax.plot(self.alpha_breakpoints, self.gamma_delta_fourier_rec, label="$\gamma_\Delta$ Fourier")
            ax.set(xlabel='alpha (rad)', ylabel='$\gamma_\Delta (H^{-1})$')
            plt.legend(loc="upper right")
            plt.show()
            fig, ax = plt.subplots()
            ax.plot(self.alpha_breakpoints, gamma_dq_polar_I, lw=2, label="$\gamma_{dq}$ input")
            ax.plot(self.alpha_breakpoints, self.gamma_dq_fourier_rec, label="$\gamma_{dq}$ Fourier")
            ax.set(xlabel='alpha (rad)', ylabel='$\gamma_{dq} (H^{-1})$')
            plt.legend(loc="upper right")
            plt.show()
        
           
        # Reconstruction of signal Ihq
        try : 
            self.Ihq
        except: 
            print('  APOLLO ERROR in calc_fourier_inverse_inductances:')
            print('  impossible to compute the Ihq reconstruction ')
            print('  please compute function calc_convergence_region first')
            exit()
            
        # trasform Ihq into polar coordinates
        Ihq = self.Ihq.flatten()
        NaN_mask = ~np.isnan(Ihq) 
        Ihq=Ihq[NaN_mask]
        self.alpha_breakpoints = math.pi/180*np.linspace(0,360,1000) 
        self.Ihq_polar_I = griddata(points_I, Ihq, (self.alpha_breakpoints, [I] * 1000), method='cubic')
        theta_star = np.interp(I, self.i_MTPA, self.theta_MTPA) #MTPA angle for I current

        
        # Reconstruct Ihq with Fourier series of inverse incremental inductances
        if  np.mean(self.ldelta[~np.isnan(self.ldelta)]) < 0: 
            print('computed Ihq_fourier using REL convention')
            a=(self.gamma_dq_fourier_coef[4]+self.gamma_delta_fourier_coef[3])/2 
            c=self.gamma_delta_fourier_coef[0]
            e=(-self.gamma_dq_fourier_coef[4]+self.gamma_delta_fourier_coef[3])/2
            self.Ihq_fourier=self.Uh/(2*math.pi*self.fh)*(a*np.sin(-2*theta_star)+c*np.sin(2*self.alpha_breakpoints-2*theta_star)+e*np.sin(4*self.alpha_breakpoints-2*theta_star))
            
        else :
            print('computed Ihq_fourier using PMAREL convention')
            a=(-self.gamma_dq_fourier_coef[4]-self.gamma_delta_fourier_coef[3])/2 
            b=-self.gamma_dq_fourier_coef[2]/2
            c=-self.gamma_delta_fourier_coef[0]
            d=-b
            e=(self.gamma_dq_fourier_coef[4]-self.gamma_delta_fourier_coef[3])/2
            self.Ihq_fourier=self.Uh/(2*math.pi*self.fh)*(a*np.sin(-2*theta_star)+b*np.sin(self.alpha_breakpoints-2*theta_star)+\
            c*np.sin(2*self.alpha_breakpoints-2*theta_star)+d*np.sin(3*self.alpha_breakpoints-2*theta_star)+e*np.sin(4*self.alpha_breakpoints-2*theta_star))
            
        #plot real and approximated Ihq
        if plot:
            plt.plot(self.alpha_breakpoints, self.Ihq_polar_I, label="$I_{hq}$ input")
            plt.plot(self.alpha_breakpoints, self.Ihq_fourier, label="$I_{hq}$ Fourier")
            plt.plot(self.alpha_breakpoints, np.zeros(self.Ihq_polar_I.size), 'gray') # Ihq=0
            plt.plot(np.interp(I, self.i_MTPA, self.theta_MTPA), 0, 'ko') #MTPA angle for I current
            plt.xlabel('alpha (rad)'); plt.ylabel('$I_{hq}$ (A)')
            plt.legend(loc="upper right")
            plt.show() 
        
        self.I=I

    def calc_Ihq_few_points(self, I='', steps=1000, plot=False):
    
        """
        This function calculates the reconstructed Ihq signal with only two gamma_delta and 
        gamma_dq measurements for a current modulus I.
        It uses the Ihq reconstructed function obtained with gamma_delta and gamma_dq Fourier analysis.
        Gamma_delta and gamma_dq measurements are for alpha = 45° and alpha =90°
        
        You set manually:
        
        I: the current modulus value I in which this function calculates the Ihq
        
        steps: the number of breakpoints to create the angle vector, default is 1000
        
        This function stores the results in the following variables:
        self.alpha_breakpoints: angles for polar coordinates
        self.Ihq_fourier_rec: the reconstructed Ihq using the Fourier series of inverse incremental inductances

        """
    
        # trasform into polar coordinates
        alpha_array = np.arctan2(self.i_q, self.i_d)
        alpha_array = np.mod(alpha_array, 2*np.pi) 
        I_array = np.sqrt(self.i_d**2 + self.i_q**2)

        # Define self.alpha_breakpoints
        self.alpha_breakpoints = math.pi/180*np.linspace(0,360,steps)
        
        # convert two-dimensional vector into one-dimensional 
        alpha_array = alpha_array.flatten()
        I_array = I_array.flatten()
        gamma_delta = self.gamma_delta.flatten()
        gamma_dq = self.gamma_dq.flatten()
        
        # remove NaN
        NaN_mask = ~np.isnan(gamma_delta) 
        I_array = I_array[NaN_mask]
        alpha_array = alpha_array[NaN_mask] 
        gamma_delta = gamma_delta[NaN_mask]
        gamma_dq = gamma_dq[NaN_mask]
        
        points_I = np.vstack((alpha_array, I_array)).T    
    
        if I=='' :
            print('  APOLLO ERROR in calc_Ihq_few_points:')
            print('  impossible to compute the Ihq with few measurements')
            print('  please set the current modulus value')
            exit()
            
        elif I < 0 or I > np.min([np.max(np.max(self.i_d)),np.max(np.max(self.i_q))]) :
            print('  APOLLO ERROR in calc_Ihq_few_points:')
            print('  impossible to compute the Ihq with few measurements')
            print('  please set I values between 0 and i_d_max or i_q_max')
            exit()
            
        # compute griddata for the current modulus value I
        gamma_delta_polar_I = griddata(points_I, gamma_delta, (self.alpha_breakpoints, [I] * steps), method='cubic')
        gamma_dq_polar_I    = griddata(points_I, gamma_dq, (self.alpha_breakpoints, [I] * steps), method='cubic')
        print('I=',I)
        
        # remove NaN
        NaN_mask = ~np.isnan(gamma_delta_polar_I) 
        gamma_delta_polar_I = gamma_delta_polar_I[NaN_mask]
        gamma_dq_polar_I = gamma_dq_polar_I[NaN_mask]
        self.alpha_breakpoints_NaN = self.alpha_breakpoints[NaN_mask] 
        
        try :
            self.Ihq
        except : 
            print('  APOLLO ERROR in calc_Ihq_few_points:')
            print('  impossible to compute the Ihq with few measurements')
            print('  please compute function calc_convergence_region first')
            exit() 
            
        # trasform Ihq into polar coordinates
        Ihq = self.Ihq.flatten()
        NaN_mask = ~np.isnan(Ihq) 
        Ihq=Ihq[NaN_mask]
        Ihq_polar_I = griddata(points_I, Ihq, (self.alpha_breakpoints, [I] * steps), method='cubic')
        # Ihq_polar_I =Ihq_polar_I [NaN_mask]
        theta_star = np.interp(I, self.i_MTPA, self.theta_MTPA) #MTPA angle for I current  
        
        # Reconstruct Ihq with two measurements of the inverse incremental inductances gamma_delta and gamma_dq
        if  np.mean(self.ldelta[~np.isnan(self.ldelta)]) < 0: 
            print('computed Ihq_fourier_2meas using REL convention')
        
            gamma_delta_90=np.interp(math.pi/2,self.alpha_breakpoints_NaN,gamma_delta_polar_I)          
            gamma_delta_45=np.interp(math.pi/4,self.alpha_breakpoints_NaN,gamma_delta_polar_I)          
            gamma_dq_45=np.interp(math.pi/4,self.alpha_breakpoints_NaN,gamma_dq_polar_I)
            
            self.gamma_delta_a0=gamma_delta_45                 
            self.gamma_delta_a2=gamma_delta_45-gamma_delta_90
            self.gamma_dq_b2=gamma_dq_45

            a_rec=(self.gamma_delta_a2+self.gamma_dq_b2)/2 
            c_rec=self.gamma_delta_a0
            e_rec=(-self.gamma_dq_b2+self.gamma_delta_a2)/2

            self.Ihq_fourier_2meas=self.Uh/(2*math.pi*self.fh)*(a_rec*np.sin(-2*theta_star)+c_rec*np.sin(2*self.alpha_breakpoints-2*theta_star)+e_rec*np.sin(4*self.alpha_breakpoints-2*theta_star))
        
        else:
            print('computed Ihq_fourier_2meas using PMAREL convention')
        
            gamma_delta_90 = np.interp(math.pi/2,self.alpha_breakpoints_NaN,gamma_delta_polar_I)
            gamma_delta_135 = np.interp(3*math.pi/4,self.alpha_breakpoints_NaN,gamma_delta_polar_I)
            gamma_dq_90 = np.interp(math.pi/2,self.alpha_breakpoints_NaN,gamma_dq_polar_I)
            gamma_dq_135 = np.interp(3*math.pi/4,self.alpha_breakpoints_NaN,gamma_dq_polar_I)
            
            self.gamma_delta_a0=gamma_delta_135
            self.gamma_delta_a2=gamma_delta_135-gamma_delta_90
            self.gamma_dq_b1=gamma_dq_90
            self.gamma_dq_b2=-(gamma_dq_135-(gamma_dq_90*np.sin(3*math.pi/4)))

            a_rec=(-self.gamma_delta_a2-self.gamma_dq_b2)/2 
            b_rec=-self.gamma_dq_b1/2
            c_rec=-self.gamma_delta_a0
            d_rec=-b_rec
            e_rec=(-self.gamma_delta_a2+self.gamma_dq_b2)/2

            self.Ihq_fourier_2meas=self.Uh/(2*math.pi*self.fh)*(a_rec*np.sin(-2*theta_star)+b_rec*np.sin(self.alpha_breakpoints-2*theta_star)+c_rec*np.sin(2*self.alpha_breakpoints-2*theta_star)
            +d_rec*np.sin(3*self.alpha_breakpoints-2*theta_star)+e_rec*np.sin(4*self.alpha_breakpoints-2*theta_star))
            
        # plot
        if plot:
              
            try :
                self.Ihq_fourier
                if I==self.I :
                    plt.plot(self.alpha_breakpoints, Ihq_polar_I, label="$I_{hq}$ input")
                    plt.plot(self.alpha_breakpoints, self.Ihq_fourier, label="$I_{hq}$ Fourier")
                    plt.plot(self.alpha_breakpoints, self.Ihq_fourier_2meas, label="$I_{hq}$ TwoMeas")
                    plt.plot(self.alpha_breakpoints, np.zeros(Ihq_polar_I.size), 'gray') # Ihq=0
                    plt.plot(np.interp(I, self.i_MTPA, self.theta_MTPA), 0, 'ko') #MTPA angle for I current
                    plt.xlabel('alpha (rad)'); plt.ylabel('$I_{hq}$ (A)')
                    plt.legend(loc="upper right")
                    plt.show()                
                else:
                    print('  APOLLO ERROR in calc_Ihq_few_points:')
                    print('  impossible to show plots with Ihq_fourier')
                    print('  please set the same current amplitude I')
                    exit() 
 
            except :
                plt.plot(self.alpha_breakpoints, Ihq_polar_I, label="$I_{hq}$ input")
                plt.plot(self.alpha_breakpoints, self.Ihq_fourier_2meas, label="$I_{hq}$ TwoMeas")
                plt.plot(self.alpha_breakpoints, np.zeros(Ihq_polar_I.size), 'gray') # Ihq=0
                plt.plot(np.interp(I, self.i_MTPA, self.theta_MTPA), 0, 'ko') #MTPA angle for I current
                plt.xlabel('alpha (rad)'); plt.ylabel('$I_{hq}$ (A)')
                plt.legend(loc="upper right")
                plt.show() 

  


    
    
    

## TODO
## 1. check comments
## 2. esempi con infittimento dei dati e interpolazione su griglia non regolare
## 3. save_maps è utile per poter salvare (e ricaricare) una mappa prima di modificarla
##    aggiungere esempio con mappa di partenza e salvataggio dopo interpolazione
##    aggiungere esempio di salvataggio prima e dopo interpolazione "Marko"
## 4. aggiungere esempi di figure con gnuplot (e funzioni esportazione da apollo)
## 5. aggiungere negli esempi jupyter (e anche in appendice della tesi?)
##    le equazioni implementate nelle varie funzioni
##    calc_apparent_inductances
##    calc_incremental_inductances
##    calc_inverse_incremental_inductances
##    calc_MTPA
##    calc_saliency
##    calc_sensored_error
##    ...più altre mancanti
## 6. Aggiungere snipet di esempio per lambda_d(id) per diverse iq 
##    e lambda_q(iq) per diverse id 
## 7. MTPV ?
## 8. coppia-velocità? potenza-velocità? Mappe efficienza?   
