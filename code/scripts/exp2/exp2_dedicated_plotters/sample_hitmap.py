import jspectroscopy as spec
import numpy as np
import matplotlib.pyplot as plt
import colorcet
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["font.family"] = "Arial"
plt.rcParams[ 'mathtext.default' ] = 'regular' 


savepath = '../../../../plots_important/sample_hitmaps.eps'


analysis_mgr = spec.dssd_analysis_manager( 'full_bkgd_tot', '../../../../storage',
                                           (32,32) , [1,1,1] )

hitmaps = analysis_mgr.load_dill( 'hitmaps' )

source_names = [ r'$^{148}Gd$', r'$^{264}Cm$' ] 


f, axarr = plt.subplots( 2, 1, figsize = ( 4, 8 ) ) 

f.suptitle( 'Source Hitmaps', fontsize = 20 ) 

f.subplots_adjust( hspace = 0.5 ) 

for i in range( 2 ) :
    idx = i + 1 
    axarr[i].set_title( source_names[i], fontsize = 15 )
    im = axarr[i].imshow( hitmaps[idx][0][3], cmap = colorcet.m_rainbow )
    divider = make_axes_locatable( axarr[i] )
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im, cax=cax)




plt.savefig( savepath, dpi = 100 ) 

plt.show() 
