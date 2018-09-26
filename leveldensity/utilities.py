################################################################################
# Author: Alec Golas                                                           #
# Date: September 24th, 2018                                                   #
# Title: utilities.py                                                          #
# Contains some math utilities                                                 #
################################################################################
import numpy as np

################################################################################
class Math:

    def dfdx_1d(self, grid, func):

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

        return dfdx

    def nearest_value_index(self, array, value):
        array = np.asarray(array)
        index = np.argmin(np.abs(array - value))

        return index

    def nearest_value(self, array, value):
        array = np.asarray(array)
        index = np.argmin(np.abs(array - value))

        return array[index]

    def summation(self, func):
        sum = np.zeros(np.shape(func))
        for i in range(len(func)):
            sum[i] = np.sum(func[0:i])
        return sum

################################################################################
class Plotting:
    pass



################################################################################
####################### End of utilities.py ####################################
