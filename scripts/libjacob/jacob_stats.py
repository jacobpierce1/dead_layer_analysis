import libjacob.jacob_math as jmath
import libjacob.jacob_pyplot as jplt
import matplotlib.pyplot as plt
import numpy as np




def linear_calibration( x, y, dy, p0, print_fit_data=0, ax=None ):
    
    # moment of truth 
    linear_fit = lambda p, x_: p[0]*x_ + p[1]

    result = jmath.jacob_least_squares( x, y, dy, p0, linear_fit )
    if result is not None:
        reduc_chisq, dof, pf, pferr = result
    else:
        print 'ERROR: fit failed to converge.'
        return 0

    # print the fit data
    if print_fit_data:
        messages = [ 'reduc_chisq', 'dof', 'pf', 'pferr' ]
        for i in range(len(result)):
            print messages[i] + ': ' + str(result[i])
    
    # plot the fit
    fit_bounds = [ min(x), max(x) ]
    x_forfit = np.linspace( fit_bounds[0], fit_bounds[1], 200*y.size )
    jplt.add_fit_to_plot( ax, x_forfit, fit_bounds, pf, pferr, linear_fit )

    # add data
    ax.errorbar( x, y, dy, fmt='None', errorevery=1 )
    
    jplt.set_linear_scale_plot_bounds( ax, x, y )
    jplt.add_legend(ax) 

    return 1



