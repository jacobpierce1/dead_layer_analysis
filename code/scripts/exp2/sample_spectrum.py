import jspectroscopy as spec
import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats as st
from sklearn.neighbors.kde import KernelDensity
import statsmodels.api as sm

# data_name = 'full_bkgd_tot_cut'
data_name = 'full_bkgd_tot'
coords = ( 0, 16, 16 ) 


y_bottom = 1e-6
hist_bin_size = 10

kde_bandwidth = 1


mgr = spec.dssd_analysis_manager( data_name, '../../../storage', (32,32), [1,1,1] )

# coords = (0, 27, 16 )


data = mgr.get_data( *coords ) 

print( data ) 

x = np.arange( 0, 3000, hist_bin_size ) 
y, bin_edges = np.histogram( data, x ) 
x = x[0:-1]

# ax = plt.axes() 

f = plt.figure( figsize = ( 6,7  ) )

ax = f.add_subplot( 211 )

ax.semilogy( x, y, ls = 'steps-mid' )
ax.set_xlim( right = 3300 )
# ax.set_ylim( top = 1e6 )


    
plt.show() 
