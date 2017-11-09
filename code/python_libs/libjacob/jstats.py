import libjacob.jmeas as meas
import libjacob.jmath as jmath
import libjacob.jpyplot as jplt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit.models import LinearModel



def linear_calibration( x, y, dy = None, dx = None, m_guess = None, b_guess = None,
                        print_results = 0, ax = None, invert = 0,
                        scatter = 0 ):

    model = LinearModel() 
            
    if m_guess is None:

        xmin = min(x)
        xmax = max(x)
        ymin = min(y)
        ymax = max(y)

        m_guess = ( ymax - ymin ) / ( xmax - xmin )

    if b_guess is None:
        b_guess = ( y[0] - m_guess * x[0] ) 
    
    model.make_params( m = m_guess, b = b_guess )


    # construct weights based on the inputs for dx and dy.
    if dy is not None:

        for i in range( len( dy ) ):
            if dy[i] == 0:
                y[i] = np.nan
                x[i] = np.nan
                dy[i] = np.nan

        if dx is None:
            weights = 1.0 / np.asarray( dy )
        else:
            raise NotImplementedError()
            # weights = 1.0 / np.sqrt( 
    else:
        if dx is None:
            weights = None
        else:
            raise NotImplementedError()

    print( x )
    print( y )
    print( weights )

    print( len( x ) )
    print( len( y ) )
    print( len( weights ) ) 
        
    model_result = model.fit( y, x = x, weights = weights, nan_policy = 'omit' )

    if print_results:
        print( model_result.fit_report() )

    # determine if we had a succussful fit. ier
    # is the return code of scipy.optimize.leastsq.
    # if this is 1, then lmfit thinks we succeeded. 
    successful_fit = ( model_result.success and ( model_result.ier < 4 ) )
                       # and ( model_result.redchi < reduc_chisq_max ) )

    if not successful_fit:
        return None

    mparam = model_result.params['slope']
    bparam = model_result.params['intercept']

    m = meas.meas( mparam.value, mparam.stderr )
    b = meas.meas( bparam.value, bparam.stderr ) 

    # error analysis built into these calculations via meas.
    if invert: 
        b = - b / m
        m = 1 / m 


    if ax is not None:
        graph_func = lambda _x : m.x * _x + b.x 
        ax_xaxis = np.array( [ xmin, xmax ] ) 
        
        ax.plot( ax_xaxis, graph_func( ax_xaxis ), label = 'Fit', color = 'r' )
        
    
    return m, b, model_result
    
    # return model_result if successful_fit else None 





    # linear_fit = lambda p, x_: p[0]*x_ + p[1]
    # return calibration( linear_fit, x, y, dy, p0, print_fit_data, ax, invert, scatter )



# def quadratic_calibration( x, y, dy, p0, print_fit_data=0, ax=None, invert=0, scatter=0 ):
    
#     quadratic_fit = lambda p, x_: p[0]*x_**2.0 + p[1]*x_ + p[2]
#     return calibration( quadratic_fit, x, y, dy, p0, print_fit_data, ax, invert, scatter )





# general calibration, pass linear / quadraditic etc function as fitfunc 
# often least_squares is called for a calibration, this function does that
# and the usual things i do when performing a calibration.
# no error handling done other than returning None if there is a failed fit.
# if invert is 1 then the linear equation that is returned is inverted.
# use scatter=1 to do a scatter instead of an errorbar

def calibration( fitfunc, x, y, dy, p0, print_fit_data=0, ax=None, invert=0, scatter=0 ):
    
    result = jmath.jleast_squares( x, y, dy, p0, fitfunc )
    if result is not None:
        reduc_chisq, dof, pf, pf_delta = result
    else:
        return None

    # print the fit data
    if print_fit_data:
        messages = [ 'reduc_chisq', 'dof', 'pf', 'pf_delta' ]
        for i in range(len(result)):
            print( messages[i] + ': ' + str(result[i]) )
    
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
