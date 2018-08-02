import jspectroscopy as spec
import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats as st
from sklearn.neighbors.kde import KernelDensity
import statsmodels.api as sm



def get_peakpos( bins, hist, guess, bounds ) :
    tmp = hist.astype( float ) 
    mask = ( bins < guess + bounds[0] ) | ( bins > guess + bounds[1] )
    tmp[ mask ] = np.nan 
    return np.nanargmax( tmp ) 
    

def place_label( ax, text, guess, bins, hist ) :

    peakpos = get_peakpos( bins, hist, guess, [ -30, 30 ] )
    peak_height = hist[ peakpos ] 
    
    ylims = ax.get_ylim()
    annotation_scale_factor = ( ylims[1] / ylims[0] )**( 0.03 )
    x = bins[ peakpos ] 
    y = annotation_scale_factor * peak_height

    print( x,y) 
    ax.text( x, y, text, fontsize = 16, ha = 'center' ) 



labels = [ r'$^8$Li$\to\alpha\alpha\beta^-$', r'$^{148}$Gd', r'$^{244}$Cm' ]
guesses = [ 763, 1600, 2903 ]     
num_labels = len( labels ) 

y_bottom = 1e-6
hist_bin_size = 10

kde_bandwidth = 1

savepath = '../../../../plots_important/'

db = spec.spectrum_db( 'full_bkgd_tot', '../../../../storage' )

# coords = (0, 27, 16 )
coords = ( 0, 16, 16 ) 

data = db.get_data( *coords ) 
x = np.arange( 0, 3000, hist_bin_size ) 
y, bin_edges = np.histogram( data, x ) 
x = x[0:-1]

# ax = plt.axes() 

f = plt.figure( figsize = ( 6,7  ) )
f.suptitle( 'Uncalibrated $^8$Li Decay Spectrum\nwith Calibration Sources', fontsize = 20  ) 

ax = f.add_subplot( 111 ) 
ax.semilogy( x, y, ls = 'steps-mid' )
ax.set_xlim( right = 3300 )
ax.set_ylim( top = 1e6 )
ax.set_xlabel( 'Channels', fontsize = 16 )
ax.set_ylabel( 'Counts', fontsize = 16  ) 
ax.tick_params(axis='both', which='major', labelsize = 14)

for i in range( num_labels ) :
    place_label( ax, labels[i], guesses[i], x, y ) 

    
# plt.savefig( savepath + 'exp2_example_spectrum.eps', dpi = 2000 )
plt.show() 
