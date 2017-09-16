## includes 

from jacob_math import xcut

import array
import numpy as np
import matplotlib.pyplot as plt

import sys
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)




# plot a histogram with my usual settings 
def plot_histo( ax, x, histo_array, plot_bounds=None, xlabel="", ylabel="", title="", logscale=0 ):
    
    # clever
    if( plot_bounds==None ):
        plot_bounds = [ f(x) for f in [min,max] ] 
    
    # ensure that we are working with arrays.
    x = np.array(x)
    histo_array = np.array(histo_array)
    
    ax.errorbar ( x, histo_array, yerr=np.sqrt(histo_array), errorevery=1)
#    ax.semilogy( range(len(histo_array)), histo_array, yerr=np.sqrt(histo_array), errorevery=1)
    # ax.errorbar ( range(len(histo_array)), histo_array, yerr=np.sqrt(histo_array), fmt='none', errorevery=10)

    if( logscale ):
        ax.set_yscale('log')

    if(ylabel):
        ax.set_ylabel (ylabel, fontsize=18)
    if(xlabel):
        ax.set_xlabel ("Energy (keV)", fontsize=18)
    if(title):
        ax.set_title (title, fontsize=22)
    
    yexpand = 3 if logscale else 1.1
    # ax.legend(loc='upper right', fancybox=True, shadow=True)
    ax.axis((plot_bounds[0], plot_bounds[1], 0, max(histo_array) * yexpand))
    # plt.show()



def add_legend(ax):
    ax.legend(loc='upper right', fancybox=True, shadow=True)



    
# set the scale so that the fractional space between the edges of the plot
# are away from the farthest data points such that xbuf, ybuf are the fractional
# size of the plot occupied by empty space.
def set_linear_scale_plot_bounds( ax, x, y, xbuf=0.20, ybuf=0.20 ):

    x1 = min(x)
    y1 = min(y)
    x2 = max(x)
    y2 = max(y)

    left = min(x) - xbuf * (x2 - x1)
    right = max(x) + xbuf * (x2 - x1)
    bottom = min(y) - ybuf * (y2 - y1)
    top = max(y) + ybuf * (y2 - y1)

    ax.axis( (left, right, bottom, top) )
  
        
        
# add fit to the plot 
def add_fit_to_plot( ax, x, fit_bounds, p, perr, fitfunc ):
    newx = xcut( x, x, fit_bounds )
    fit = ax.plot( newx, fitfunc(p, newx), '-r', label="Fit")
    plt.setp(fit[0], linewidth=2)



# save in high quality 
def saveplot_high_quality( directory, fname ):
    plt.savefig( directory + '/' + fname + '.eps', format='eps', dpi=2000)

def saveplot_med_quality( directory, fname ):
    plt.savefig( directory + '/' + fname + '.eps', format='eps', dpi=500)


def saveplot_low_quality( directory, fname ):
    plt.savefig( directory + '/' + fname + '.png', format='png') #, dpi=2000)



