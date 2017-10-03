# the purpose of this module is to interpolate stopping power data
# obtained from NIST and provide a function to numerically integrate
# energy loss due to scattering through a known distance using this
# interpolated function. it is a class instead of a module so that you
# can use it on any dataset. to use, call constructor and pass in
# arrays (energy, stopping power).

import numpy as np
import pandas as pd 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import libjacob.jutils as jutils



TEST_STOPPING_POWER_INERPOLATION = 1


class stopping_power_interpolation( object ):


    # constructor takes x (presumably energy) and y (stopping power).
    def __init__( self, x, y, xbounds=None ):

        if len(x) != len(y):
            print 'ERROR: lengths of input arrays must be the same.'
            self = None
            return self


        # this is used to check whether input data is in interpolated region.
        self.xbounds = xbounds

        # perform cut if necessary
        if xbounds is not None:
            y = xcut( x, y, xbounds )
            x = xcut( x, x, xbounds )
            
        # these are used for plotting.
        self.xdata = x
        self.ydata = y
        self.xmin = min(x)
        self.xmax = max(x) 

        # find the interpolation function
        interp = interp1d( x, y, kind='cubic' )

        # remove possibility of unwanted side effects by not allowing
        # more arguments to be passed (which is possible before this).
        self._interp = lambda x: interp( x )


        
    # apply the interpolation. return none if error.
    def interp( self, x ):

        # if x < self.xmin or x > self.xmax:
        #     print 'ERROR: x is out of interpolation range.'
        #     return None
        
        return self._interp(x) 
        

    
    # plot the interpolation on test data. it is built to dominate the plot
    def plot( self, ax=None ):

        ax_not_supplied = ax is None

        # create a new axis if not supplied.
        if ax_not_supplied:
            plt.clf()
            ax = plt.axes()

        # add data to plot and spline with 100 * density of data.
        ax.scatter( self.xdata, self.ydata )
        newx = np.linspace( self.xmin, self.xmax, len(self.xdata) * 100 )
        ax.plot( newx, self.interp(newx), '-r' )

        # presumably if you pass your own axes then you will do this on your own.
        if ax_not_supplied:
            plt.show()



    # input: initial energy E0, distance traveled.  output: final
    # energy Ef, calculated assuming that the distance traveled is
    # small enough that no numerical integration is required.
    def energy_loss_naive( self, E0, dist ):
        return self.interp(E0) * dist
        


    # input: initial energy E0, distance traveled.  output: final
    # energy Ef, calculated using numerical integration.
    def energy_loss( self, E0, dist ):
        pass

    

    @classmethod
    def from_nist_file( cls, filename, density, xbounds=None ):
        data =  pd.read_table( filename, delim_whitespace=1, skiprows=11, usecols=[0,3],
                               names=['energy', 'stopping_power'] )

        
        data['energy'] *= 1000  # change from MeV to keV
        data['stopping_power'] *= density  # nist files are normalized to density.
        
        return cls( data['energy'], data['stopping_power'], xbounds )



# toggle to test the class on some test data when this script is run. mainly intended for
# developmnet purposes. when using this for an interpolation you should plot it to make sure
# the interpolation is good, not using this test data. 
_TEST_STOPPING_POWER_INERPOLATION = 0


def _test_stopping_power_interpolation():
    energies = np.array( [0,2,4,6,8] )
    stopping_powers = np.sin( energies )
    x = stopping_power_interpolation( energies, stopping_powers )
    x.plot()


    
if _TEST_STOPPING_POWER_INERPOLATION:
    _test_stopping_power_interpolation()
