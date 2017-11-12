## includes 

from .jmath import xcut

import array
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

import sys
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)




# plot a histogram with my usual settings 
def plot_histo( ax, x, histo_array, plot_bounds=None, xlabel="", ylabel="", title="", logscale=0 ):
    
    # clever
    if plot_bounds is 'minmax':
        plot_bounds = [ f(x) for f in [min,max] ] 

        
    # ensure that we are working with arrays.
    x = np.array(x)
    histo_array = np.array(histo_array)

    ax.errorbar( x, histo_array, yerr=np.sqrt(histo_array), errorevery=1)
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

    if plot_bounds is not None:
        ax.axis((plot_bounds[0], plot_bounds[1], 0, max(histo_array) * yexpand))



        

def add_legend(ax, loc=0):
    if loc==1:
        locstr = 'upper right'
    elif loc==2:
        locstr = 'lower right'
    elif loc==-1:
        locstr = 'upper left'
    elif loc==-2:
        locstr = 'lower left'
    else:
        print( 'USAGE for jacob_pyplot.add_legend : loc=1 -> upper right, loc=2 -> lower right, ' + 
                   'loc=-1 -> upper left, loc=-2 -> lower left' )

        return 0 

    # remove duplicates
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    # plt.legend(
    
    ax.legend(by_label.values(), by_label.keys(), loc=locstr, fancybox=True, shadow=True)
    return 1




    
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
def add_fit_to_plot( ax, x, fit_bounds, fitfunc, linewidth=3, color='-r' ):
    xfit = xcut( x, x, fit_bounds )
    fit = ax.plot( xfit, fitfunc( xfit ), color, label="Fit")
    plt.setp( fit[0], linewidth = linewidth )


    
    

# default plot params 
def plot ( ax, x, y, xerr=None, yerr=None,
           xlabel=None, ylabel=None,
           title=None, color='r', errorevery=1,
           leglabel = None,
           legloc = None, logscale=0 ):

    # default: don't plot errorbars.
    plot_errorbar = 0
    
    if xerr is not None or yerr is not None:

        plot_errorbar = 1

        # omit fmt parameter for python 3, use fmt = 'none' for python 2.
        ax.errorbar ( x, y, xerr = xerr, yerr = yerr,
                      errorevery = errorevery, fmt = 'none',
                      ecolor = color, label = leglabel )

    else:
        ax.plot( x, y, color=color )


    # handle the title / etc
    if( logscale ):
        ax.set_yscale('log')

    if ylabel is not None:
        ax.set_ylabel (ylabel, fontsize=18)

    if xlabel is not None:
        ax.set_xlabel (xlabel, fontsize=18)

    if title is not None:
        ax.set_title( title, fontsize=22 )

    if not plot_errorbar:
        set_linear_scale_plot_bounds( ax, x, y )

    else:
        if xerr is not None:
            xbounds = [ min(x - xerr), max( x + xerr) ]

        else:
            xbounds = [ min(x), max(x) ]
            
        if yerr is not None:
            ybounds = [ min(y - yerr), max( y + yerr) ]

        else:
            ybounds = [ min(y), max(y) ] 

    if legloc != None:
        add_legend( ax, legloc ) 
        
            
#     set_linear_scale_plot_bounds( ax, xbounds, ybounds  )
        
    # add_legend( ax, legloc )    



    
            
# save in high quality 
def saveplot_high_quality( directory, fname ):
    plt.savefig( directory + '/' + fname + '.eps', format='eps', dpi=2000)

def saveplot_med_quality( directory, fname ):
    plt.savefig( directory + '/' + fname + '.eps', format='eps', dpi=500)


def saveplot_low_quality( directory, fname ):
    plt.savefig( directory + '/' + fname + '.png', format='png') #, dpi=2000)



