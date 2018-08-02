import jspectroscopy as spec
import numpy as np
import matplotlib.pyplot as plt
import colorcet
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import os 
plt.rcParams["font.family"] = "Arial"
plt.rcParams[ 'mathtext.default' ] = 'regular' 


sys.path.append('../../')
import new_deadlayer_estimator as dl_estimator

savepath = '../../../../plots_important/sample_hitmaps.eps'


analysis_mgr = spec.dssd_analysis_manager( 'full_bkgd_tot', '../../../../storage',
                                           (32,32) , [1,1,1] )


det = 0
fstrip = 25

data = analysis_mgr.load_dill( 'means' )
secants = analysis_mgr.load_dill( 'secant_matrices' ) 


ndata = 32
# 24.6542711918525
redchisqr = 0.8805096854233037
params = np.array( [ 2.08919995, -77.49714507,  28.61558523,  26.36277313 ] ) 
a, b = params[:2]
dl = params[2:]



actual_energies = [  [ 3182.690 ] ,
                     [ 0.231 * 5762.64  + 0.769 * 5804.77 ] ]


source_names = [ r'$^{148}Gd$', r'$^{264}Cm$' ] 


f, axarr = plt.subplots( 2, 1, figsize = ( 4, 8 ) ) 

f.suptitle( 'Source Means', fontsize = 20 ) 

f.subplots_adjust( hspace = 0.5 ) 

for i in range( 2 ) :
    idx = i + 1 
    axarr[i].set_title( source_names[i], fontsize = 15 )
    tmp_data = data[idx][0][0][25]

    x = secants[i][0][25]
    y = tmp_data.x
    dy = tmp_data.dx
    
    axarr[i].errorbar( x, y, dy, ls = 'none' ) 

    mask = ~ ( np.isnan( x ) | np.isnan( y ) )
                    
    min_sectheta = min( x[ mask ] ) 
    max_sectheta = max( x[ mask ] ) 
    min_chan = ( ( actual_energies[i][0] - dl[i] * min_sectheta ) - b ) / a
    max_chan = ( ( actual_energies[i][0] - dl[i] * max_sectheta ) - b ) / a
    axarr[i].plot( [ min_sectheta, max_sectheta ], [ min_chan, max_chan ] ) 
                    

axarr[0].text( 0.7, 0.9, r'$\tilde \chi ^2 = %.2f$' % redchisqr,
               transform = axarr[0].transAxes,
               fontsize = 12 )
    

# plt.savefig( savepath, dpi = 100 ) 

plt.show() 
