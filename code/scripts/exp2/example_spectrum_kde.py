import jspectroscopy as spec
import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats as st
from sklearn.neighbors.kde import KernelDensity
import statsmodels.api as sm


use_kde = 0
y_bottom = 1e-6
hist_bin_size = 20

kde_bandwidth = 1

db = spec.spectrum_db( 'full_bkgd_tot', '../../../storage' )

coords = (0 , 16, 16 )

if not use_kde : 
    data = db.get_data( *coords ) 
    x = np.arange( 0, 3000, hist_bin_size ) 
    y, bin_edges = np.histogram( data, x ) 
    x = x[0:-1]
    
    # ax = plt.axes() 
    
    f = plt.figure( figsize = ( 12, 10 ) )

    ax1 = f.add_subplot( 211 )
    ax1.semilogy( x, y, ls = 'steps-mid' )
    ax1.set_xlim( right = 3300 ) 


    ax2 = f.add_subplot( 223 )
    ax2.semilogy( x, y, ls = 'steps-mid' )
    ax2.set_xlim( left = 1325, right = 1630 ) 
    ax2.set_title( 'Gd Peak' ) 

    
    ax3 = f.add_subplot( 224 )
    ax3.semilogy( x, y, ls = 'steps-mid' )
    ax3.set_xlim( left = 2420, right = 2930 ) 
    ax3.set_title( 'Cm Peak' ) 


else :
    x = np.linspace( 0, 3000, 10000 )
    data = db.get_data( *coords )

    # kde = st.gaussian_kde( data, kde_bandwidth )
    # y = kde( x )
    # print( kde.factor )

    x = x.reshape( -1 , 1 ) 
    data = data.reshape( -1, 1  )
    kde = KernelDensity(kernel='gaussian', bandwidth = kde_bandwidth ).fit( data )
    y = np.exp( kde.score_samples( x ) )
    x = x[:,0]

    print( y )
    print( y.shape ) 
    print( x.shape ) 
    print( kde.bandwidth ) 
    # kde = sm.nonparametric.KDEUnivariate( data )
    # kde.fit( kernel = 'uni', fft = 0 )
    # y = np.vectorize( kde.evaluate )( x ) 

    f = plt.figure( figsize = ( 12, 10 ) )

    ax1 = f.add_subplot( 211 )
    ax1.semilogy( x, y, ls = 'steps-mid' )
    ax1.set_xlim( right = 3300 ) 
    ax1.set_ylim( bottom = y_bottom, top = 1 )

    ax2 = f.add_subplot( 223 )
    ax2.semilogy( x, y, ls = 'steps-mid' )
    ax2.set_xlim( left = 1325, right = 1630 ) 
    ax2.set_title( 'Gd Peak' ) 
    ax2.set_ylim( bottom = y_bottom, top = 1 )
    
    ax3 = f.add_subplot( 224 )
    ax3.semilogy( x, y, ls = 'steps-mid' )
    ax3.set_xlim( left = 2420, right = 2930 ) 
    ax3.set_title( 'Cm Peak' )
    ax3.set_ylim( bottom = y_bottom, top = 1 )
    
    # ax.semilogy( x, y )
    
    

plt.show() 
