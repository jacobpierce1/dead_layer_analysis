import numpy as np

# import deadlayer_helpers.geometry as geom
# import deadlayer_helpers.data_handler as data
import bpt

import matplotlib.pyplot as plt

import jspectroscopy as spec
import jutils 
# import libjacob.jmeas as meas 

import scipy.interpolate



plot_fstrip = 0




# name = 'full_bkgd_tot'
name = 'alpharun11-19'
detnum = 1


if plot_fstrip : 
    fstrip = 4
    bstrips = [ 10, 20, 30]
    num_strips = len( bstrips ) 
    channels = [ 1400, 3100 ]
else :
    bstrip = 15
    fstrips = [ 3,29 ]
    num_strips = len( fstrips ) 
    channels = [ 0, 3000 ] 

analysis_mgr = spec.dssd_analysis_manager( name, '../../storage', (32,32), [ 1,1,1] ) 


cols = 'bgrk'

peakzoom = None




plt.figure( figsize=(10,5) )

ax = plt.axes()

for i in range( num_strips ) :

    if plot_fstrip : 
        data = analysis_mgr.get_data( detnum, fstrip, bstrips[i] )
        label = bstrips[i]
    else :
        data = analysis_mgr.get_data( detnum, fstrips[i], bstrip, mode = 'eback' )
        label = fstrips[i]
        
        
    x = np.arange( channels[0], channels[1], 2 ) 
    y, tmp = np.histogram( data, bins = x )
    x = x[:-1] 

    
    ax.semilogy( x, y, c = cols[i], label = label,
                 ls='steps-mid', linewidth = 0.3 )

    # ax.set_title( '%s: det=%s, fstrip=%d, bstrips=%s' % ( name, str( detnum ),
    #                                                       fstrip, str(bstrips) ) ) 
    
    # print( ( f,b ) )
    # print( 'mu values:' )
    # print( [ channels[i][j][ f, b ] for i in range(3) for j in range(2) ] )

           
ax.legend

ax.set_xlim( channels )

ax.legend( loc = 'best' ) 

plt.show() 
