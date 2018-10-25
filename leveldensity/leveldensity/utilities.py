################################################################################
# Author: Alec Golas                                                           #
# Date: September 24th, 2018                                                   #
# Title: utilities.py                                                          #
# Contains some math utilities                                                 #
################################################################################
from __future__ import print_function
import matplotlib
matplotlib.use('agg')

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
################################################################################
class Math:

    """Numerical operations used in the level density model calculations and
    data analysis"""

    def __init__(self, x_data=None, y_data=None, array=None, value=None):

        self.x_data = x_data
        self.y_data = y_data
        self.array = array
        self.value = value

    def dfdx_1d(self):

        if self.x_data is None:
            raise TypeError('Independent variable data required')
        if self.y_data is None:
            raise TypeError('Dependent variable data required')

        grid = self.x_data
        func = self.y_data

        assert np.shape(grid) == np.shape(func)

        step = np.zeros(np.shape(grid))
        df = np.zeros(np.shape(func))
        dfdx = np.zeros(np.shape(func))

        step[0] = grid[1] - grid[0]
        step[-1] = grid[-1] - grid[-2]
        step[1:-1] = grid[2:] - grid[:-2]

        df[0] = func[1] - func[0]
        df[-1] = func[-1] - func[-2]
        df[1:-1] = func[2:] - func[:-2]

        dfdx = df/step
        dx = step
        dx[1:-1] = step[1:-1]/2
        x = np.cumsum(dx)
        deriv_array = np.asarray([x, dfdx])

        return deriv_array

    def integral(self, grid, func):
        ## If spacings is not linear, then do integral_nonlinear#
        pass

    def integral_linear(self):

        if self.x_data is None:
            raise TypeError('Independent variable data required')
        if self.y_data is None:
            raise TypeError('Dependent variable data required')

        grid = self.x_data
        func = self.y_data

        dx = grid[1] - grid[0]


        trapezoid = func*dx
        trapezoid[0] = trapezoid[0]/2.0
        trapezoid[-1] = trapezoid[-1]/2.0

        integral = np.cumsum(trapezoid)
        return integral

    def nearest_value_index(self, array=None, value=None):

        if array is None:
            array=self.array
        if value is None:
            value=self.value

        if len(array)==1:
            index = 0
        else:
            array = np.asarray(array)
            index = np.argmin(np.abs(array - value))

        return index

    def nearest_value(self, array=None, value=None):

        if array is None:
            array=self.array
        if value is None:
            value = self.value

        array = np.asarray(array)
        index = np.argmin(np.abs(array - self.value))

        return array[index]

    def normalize(self, array=None):

        if array is None:
            array = self.array
        array = np.asarray(array)
        sum_array = np.sum(array)
        norm = array/sum_array
        
        return norm

################################################################################
class CurveFitting(Math):

    def __init__(self, x_data, y_data, global_params):
        Math.__init__(self, x_data, y_data)
        self.global_params = global_params


    def smoothing_2d(self, window, mode='polynomial'):

        if self.x_data is None:
            raise TypeError('Independent variable data required')
        if self.y_data is None:
            raise TypeError('Dependent variable data required')
        popt=None
        if mode == 'polynomial':
            curve = self.polynomial_fit
            params = np.zeros((len(self.y_data), 3))
        elif mode == 'exponential':
            curve = self.exponential_fit
            params = np.zeros((len(self.y_data), 3))
        else:
            raise TypeError('Mode must be "polynomial" or "exponential"')

        half_win = int((window - 1)/2)

        inner_bounds = range(half_win, len(self.x_data) - half_win)
        y_out = np.zeros(np.shape(self.y_data))

        y_out[:half_win] = self.y_data[:half_win]
        y_out[-half_win:] = self.y_data[-half_win:]

        popt = [0.,0.,0.]

        for i in inner_bounds:

            x_win = self.x_data[i - half_win:i + half_win]
            y_win = self.y_data[i - half_win:i + half_win]

            try:
                popt, pcov = curve_fit(curve, x_win, y_win, p0=p0,bounds=(-2,2), maxfev=1e9)

            except:
                y_est = curve(self.x_data[i], *popt)
                y_out[i] = y_est

            p0 = popt
            y_est = curve(self.x_data[i], *popt)
            params[i] = popt
            y_out[i] = y_est
        
        return (np.asarray([self.x_data, y_out]) , params)

    def polynomial_fit(self, x, A, B, C):

        A_global, B_global, C_global = self.global_params
        A = A_global + A_est
        B = B_global # + B_est
        C = C_global # + C_est

        return A*x**2 + B*x + C

    def exponential_fit(self, x, A_est, B_est, C_est):

        A_global, B_global, C_global = self.global_params

        A = A_global + A_est
        B = B_global #+ B_est
        C = C_global #+ C_est

        return A*np.exp(B*(x-C))


################################################################################
class Plotting:

    """
    Utility to plot data from the different modules
    """

    def __init__(self, x, y,):
        pass

####################### End of utilities.py ####################################
