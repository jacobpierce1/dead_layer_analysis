import numpy as np

import bpt

import matplotlib.pyplot as plt

# import jspectroscopy as spec
import jutils 

import scipy.interpolate





# db_name = 'alpharun20-30'
# db_name = 'full1071-1078' 
db_name = 'full_bkgd_tot'
detnum = 1
#channels = [ 1400, 3100 ]
channels = [ 2900, 6000 ]
f = 28
bstrips = [ 12]



# db_name = 'det3_cent'
# detnum = None
# channels = [2600,3100]
# f = 15
# bstrips = [ 13, 15, 17 ]
# mode = 0


cols = 'bgrk'

peakzoom = None


db = spec.spectrum_db( '../../.../storage/', db_name )


channels = db.load_dill( 'primary_peaks' ) 



plt.figure( figsize=(10,5) )

ax = plt.axes()

for i in range( len( bstrips ) ) :

    b = bstrips[i]

    x, y, dy = bpt.data_fetcher( '../../../bpt-data/extracted_root_tree_data/',
                                 db_name,
                                 detnum, f, b )

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
