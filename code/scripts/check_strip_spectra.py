import numpy as np

import deadlayer_helpers.geometry as geom
import deadlayer_helpers.data_handler as data


import matplotlib.pyplot as plt

import jspectroscopy as spec
import libjacob.jutils as jutils
import libjacob.jmeas as meas 

import scipy.interpolate




# create and populate the histogram array 
def data_fetcher( name, x, y ) :

    xaxis = np.arange( 5000 )
        
    infile = ( '../../data/extracted_ttree_data/'
               + name + '/%s_%d_%d.bin' % ( name, x, y ) )

    # print( infile ) 
    
    efront_histo = np.zeros( xaxis.size )

    if not data.construct_histo_array( infile, efront_histo ) :
        print( 'error: couldnt open file' ) 

    dy = np.sqrt( efront_histo )
    dy[ dy==0 ] = 1 
        
    return ( xaxis, efront_histo, dy ) 






db_name = 'angled'
f = 16
bstrips = [ 2, 16, 31 ]
cols = [ 'r', 'g', 'b' ] 

db_path = '../../storage/databases/'
db = spec.spectrum_db( db_path + db_name )
channels = db.read_mu_values( '../../storage/mu_values/' + db_name
                              + '_mu_values.bin' ) 



plt.figure( figsize=(10,5) )

ax = plt.axes()

for i in range( len( bstrips ) ) :

    b = bstrips[i]

    x, y, dy = data_fetcher( db_name, f, b  )

    ax.semilogy( x, y, c = cols[i], label = str( b ),
                 ls='steps-mid', linewidth = 0.3 )

    # print( ( f,b ) )
    # print( 'mu values:' )
    # print( [ channels[i][j][ f, b ] for i in range(3) for j in range(2) ] )

           
ax.legend

ax.set_xlim( [2550,3100] )

ax.legend( loc = 'best' ) 

plt.show() 
