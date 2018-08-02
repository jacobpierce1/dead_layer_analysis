import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import jspectroscopy as spec
import numpy as np
import colorcet



analysis_mgr = spec.dssd_analysis_manager( 'full_bkgd_tot', '../../../../storage',
                                           (32,32) , [1,1,1] )

redchisqr = [ 0.905, 1.11 ]

hitmaps = analysis_mgr.load_dill( 'hitmaps' )

hitmaps_resid = analysis_mgr.load_dill( 'hitmaps_resid' ) 

f, axarr = plt.subplots( 1, 2, figsize = (10,5), squeeze = 0 ) 

axarr[0,0].tick_params(axis='both', which='major', labelsize = 14)
axarr[0,1].tick_params(axis='both', which='major', labelsize = 14)

f.subplots_adjust( wspace = 0.5 )

source_names = [ '$^{148}$Gd Hitmap', '$^{264}$Cm' ]

det = 3

title = 'Det %d: $^{148}$Gd Hitmap and Point Source Fit' % ( det )

f.suptitle( title, fontsize = 20 )

for i in range(1) :
    idx = i + 1 
    hitmap = hitmaps[idx][0][det].T
    resid = hitmaps_resid[idx][0][det].T 

    cmap = colorcet.m_rainbow
    cmap.set_bad('black',1.)
    im = axarr[ i, 0 ].imshow( hitmap, cmap = cmap, origin = 'lower'  )
    divider = make_axes_locatable( axarr[i,0] )
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im, cax=cax)
    cax.tick_params(axis='both', which='major', labelsize = 14)
    
    axarr[i,0].set_title( source_names[i], fontsize = 16 ) 
            
        
    cmap = colorcet.m_diverging_bkr_55_10_c35
    cmap.set_bad('white',1.)
    im = axarr[ i,1 ].imshow( resid, cmap = cmap, origin = 'lower' )
    divider = make_axes_locatable( axarr[i,1] )
    cax = divider.append_axes("right", size="5%", pad=0.1)
    f.colorbar(im, cax=cax)
    cax.tick_params(axis='both', which='major', labelsize = 14)
        
    axarr[i,1].set_title( r'Fit Residuals: $\tilde \chi^2 = %.2f$'
                          % ( redchisqr[i] ), fontsize = 16 )


savepath = '../../../../plots_important/exp2_hitmaps.eps'
plt.savefig( savepath, format = 'eps', dpi = 1000 )
    
plt.show() 

