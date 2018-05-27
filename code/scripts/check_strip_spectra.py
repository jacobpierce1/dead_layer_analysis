import numpy as np

import deadlayer_helpers.geometry as geom
import deadlayer_helpers.data_handler as data


import matplotlib.pyplot as plt

import jspectroscopy as spec
import libjacob.jutils as jutils
import libjacob.jmeas as meas 

import scipy.interpolate





# create and populate the histogram array 
def data_fetcher( name, detnum, x, y, mode  ) :

    xaxis = np.arange( 5000 )

    if detnum is not None : 
        infile = ( '../../data/extracted_ttree_data/'
                   + name + '/%s_%d_%d_%d.bin' % ( name, detnum, x, y ) )

    else :
        infile = ( '../../data/extracted_ttree_data/'
                   + name + '/%s_%d_%d.bin' % ( name, x, y ) )
    
    # print( infile ) 
    
    efront_histo = np.zeros( xaxis.size )

    if not data.construct_histo_array( infile, efront_histo, mode ) :
        print( 'error: couldnt open file' ) 

    dy = np.sqrt( efront_histo )
    dy[ dy==0 ] = 1 
        
    return ( xaxis, efront_histo, dy ) 






# db_name = 'alpharun20-30'
db_name = 'alpharun11-19' 
detnum = 3
channels = [ 1400, 3100 ]
f = 15
bstrips = [ 7,8, 9, 10 ]
mode = 0


# db_name = 'det3_cent'
# detnum = None
# channels = [2600,3100]
# f = 15
# bstrips = [ 13, 15, 17 ]
# mode = 0


cols = 'bgrk'

peakzoom = None


# db_path = '../../storage/databases/'
# db = spec.spectrum_db( db_path + db_name )
# channels = db.read_values( '../../storage/mu_values/' + db_name
#                               + '_mu_values.bin' ) 



plt.figure( figsize=(10,5) )

ax = plt.axes()

for i in range( len( bstrips ) ) :

    b = bstrips[i]

    x, y, dy = data_fetcher( db_name, detnum, f, b, mode  )

    ax.semilogy( x, y, c = cols[i], label = str( b ),
                 ls='steps-mid', linewidth = 0.3 )

    ax.set_title( '%s: det=%s, fstrip=%d, bstrips=%s' % ( db_name, str( detnum ),
                                                          f, str(bstrips) ) ) 
    
    # print( ( f,b ) )
    # print( 'mu values:' )
    # print( [ channels[i][j][ f, b ] for i in range(3) for j in range(2) ] )

           
ax.legend

ax.set_xlim( channels )

ax.legend( loc = 'best' ) 

plt.show() 
