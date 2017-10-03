import libjacob.jacob_math as jmath
import libjacob.jacob_pyplot as jplt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd




def linear_calibration( x, y, dy, p0, print_fit_data=0, ax=None, invert=0, scatter=0 ):
    
    linear_fit = lambda p, x_: p[0]*x_ + p[1]
    return calibration( linear_fit, x, y, dy, p0, print_fit_data, ax, invert, scatter )



def quadratic_calibration( x, y, dy, p0, print_fit_data=0, ax=None, invert=0, scatter=0 ):
    
    quadratic_fit = lambda p, x_: p[0]*x_**2.0 + p[1]*x_ + p[2]
    return calibration( quadratic_fit, x, y, dy, p0, print_fit_data, ax, invert, scatter )





# general calibration, pass linear / quadraditic etc function as fitfunc 
# often least_squares is called for a calibration, this function does that
# and the usual things i do when performing a calibration.
# no error handling done other than returning None if there is a failed fit.
# if invert is 1 then the linear equation that is returned is inverted.
# use scatter=1 to do a scatter instead of an errorbar
def calibration( fitfunc, x, y, dy, p0, print_fit_data=0, ax=None, invert=0, scatter=0 ):
    result = jmath.jacob_least_squares( x, y, dy, p0, fitfunc )
    if result is not None:
        reduc_chisq, dof, pf, pf_delta = result
    else:
        return None

    # print the fit data
    if print_fit_data:
        messages = [ 'reduc_chisq', 'dof', 'pf', 'pf_delta' ]
        for i in range(len(result)):
            print messages[i] + ': ' + str(result[i])
    
    # plot the fit
    if ax is not None:
        fit_bounds = [ min(x), max(x) ]
        x_forfit = np.linspace( fit_bounds[0], fit_bounds[1], 200*y.size )
        jplt.add_fit_to_plot( ax, x_forfit, fit_bounds, pf, pf_delta, fitfunc )

        # add data
        if not scatter:
            ax.errorbar( x, y, dy, fmt='None', errorevery=1 )
        else:
            ax.scatter( x, y )
            jplt.set_linear_scale_plot_bounds( ax, x, y )
        jplt.add_legend(ax, 2) 

    # invert everything if invert is supplied 
    if invert:
        pf_new = [ 1.0 / pf[0], 0 - pf[1] / pf[0] ]
        pf_delta_new = [ pf_delta[0] / pf[0]**2,
                         ( pf[1] / pf[0] * np.linalg.norm( [ pf_delta[i] / pf[i] for i in range(2) ] ) ) ] 

        pf[:] = pf_new[:]
        pf_delta[:] = pf_delta_new[:]
    

    calibration_function = lambda x_: fitfunc( pf, x_ ) 
        
    # return the calibrated function pointer so that it can be used. 
    return pd.Series( [calibration_function, reduc_chisq, dof, pf, pf_delta, invert ],
                    index=['f', 'reduc_chisq', 'dof', 'pf', 'pf_delta', 'inverted' ] )
