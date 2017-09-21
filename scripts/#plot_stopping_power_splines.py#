# the purpose of this script is to plot cubic splines of the
# stopping power data that are used for estimation of dead layer
# energy loss in order to make sure that they look good.

import libjacob.jacob_pyplot as jplt
import deadlayer_helpers.stopping_power_interpolation as stop

import matplotlib.pyplot as plt
import pandas as pd


# config: stopping power files 

file_si = '../data/stopping_power_data/alpha_in_si_stopping_power.txt'
density_si = 2.328  # g/mL

def _run():

    plt.clf()
    ax = plt.axes()
        
    stop_interp_si = stop.stopping_power_interpolation.from_nist_file(
        file_si, density_si )

    stop_interp_si.plot( ax )

    # will look at interpolations in this region.
    bounds = [ 3750, 6250 ]
    ax.set_xlim( bounds )
    ax.set_ylim( [ 500, 2000 ] )
    
    # set labels
    ax.set_title( 'Stopping Power Interpolation' )
    ax.set_xlabel( 'Energy (keV)' )
    ax.set_ylabel( 'Interpolated Stopping Power' )
    

    plt.show()


_run()
