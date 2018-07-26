import numpy as np
import scipy.optimize 
import scipy.interpolate 
import numdifftools
import jspectroscopy as spec
import sys 
from jutils import meas




gd_energy =  3182.690

cm_energies = [ 5762.64,  5804.77 ] 
cm_probs = [ 0.231,  0.769 ]

analysis_mgr = spec.dssd_analysis_manager( 'full_bkgd_tot', '../../../storage', (32,32),
                                           [1,1,1] )

data = analysis_mgr.load_dill( 'agg_dl_sums' ) 

print( data )


density_si = 2.328 # g / cm^3

def construct_si_stopping_power_interpolation( plot = 0 ) :

    # data = np.loadtxt( '../../data/stopping_power_data/alpha_stopping_power_si.txt',
    #                   skiprows = 10, unpack = 1 )

    energy, stopping_power = np.loadtxt(
        '../../../data/stopping_power_data/alpha_stopping_power_si.txt',
        usecols = [0,3], unpack = 1 )

    energy *= 1000 
    stopping_power *= density_si * 1000 * 100 / 1e9

    print( 'stopping power interp data:' )
    print( 'energy: ', energy )
    print( 'stopping: ', stopping_power )
    
    interp = scipy.interpolate.interp1d( energy, stopping_power, kind = 'cubic' )

    if plot :
        ax = plt.axes()
        interp_axis = np.linspace( min( energy ), max( energy ), 100 )
        ax.scatter( energy, stopping_power, color='r' )
        ax.plot( interp_axis,
                 interp( interp_axis ),
                 color = 'b' )
        plt.show()
        return 1
    
    return interp
    




stop = construct_si_stopping_power_interpolation()

stopping_powers = np.array( [ stop( gd_energy ),
                              np.dot( cm_probs, stop( cm_energies ) ) ] )


source_ids = np.array(
[
    [
        [ 0, 1, 0, 0 ],
        [ 0, 1, 1, 1 ]
    ],
    [
        [ 0, 0, 1, 0 ],
        [ 1, 1, 0, 0 ]
    ]
] )



bad_data_mask = np.array(
[
    [
        [ 0, 0, 0, 0 ],
        [ 0, 0, 0, 0 ]
    ],
    [
        [ 1, 0, 0, 0 ],
        [ 1, 0, 0, 0 ]
    ]
], dtype = bool )




data[ bad_data_mask ] = meas.nan




print( data )


def objective( params ) :

    deadlayers = params[:4]
    sources = params[4:]
    sources = sources.reshape((2,2))
    
    resid = np.zeros((2,2,4))

    for i in range(2) :
        for j in range(2) : 
            for d in range(4) :
                source_id = source_ids[i,j,d]
                resid[i,j,d] = ( deadlayers[d] * stopping_powers[i]
                                 + sources[i, source_id] - data[i,j,d].x  )

    resid /= data.dx

    # print( resid ) 

    resid = resid[ ~np.isnan( resid ) ].flatten()
    # print( resid ) 
    return np.sum(  resid ** 2 )
    # return resid


params_guess = np.array( [ 100.0, 100.0, 100.0, 100.0, 15.0, 15.0, 15.0, 15.0 ] ) 



# ret = scipy.optimize.minimize( objective, params_guess, method = 'nelder-mead',
#                                       options = { 'maxiter' : 1e7, 'xatol' : 1e-12, 'fatol' : 1e-12 } ) 
# print( ret ) 

# hessian_func = numdifftools.Hessian( objective, full_output=True )
# hessian, info = hessian_func( ret.x ) 

# try: 
#     cov = np.linalg.inv( hessian ) 
# except( np.linalg.LinAlgError ) :
#     self.cov = None 
#     self.success = 0
    
        
# params_result = ret.x
# chisqr = ret.fun

# nfree = len( data.flatten() ) - len( params_result ) 
# redchisqr = chisqr / nfree
# print( 'chisqr: ', chisqr ) 
# print( 'redchisqr: ', redchisqr ) 

# params_result_errors = np.sqrt( np.diag( cov ) * redchisqr )
# print( 'result: ', params_result )
# print( 'errors: ', params_result_errors ) 








# ret = scipy.optimize.leastsq( objective, params_guess,
#                               full_output = 1, ftol = 1e-12,
#                               xtol = 1e-12, maxfev = int(1e7) )

# params_result, cov, info, msg, status = ret

# chisqr = np.sum( info['fvec']**2 )
# nfree = len( data ) - len( params_guess ) 
# redchisqr = chisqr / nfree

# params_result_errors = np.sqrt( np.diag( cov ) * redchisqr )


# print( redchisqr ) 
# print( params_result )
# print( params_result_errors ) 




ret = scipy.optimize.basinhopping( objective, params_guess, disp = 0, niter = 100,
                                   minimizer_kwargs = { 'ftol' : 1e-12 } ) 
print( ret ) 

cov = ret.lowest_optimization_result.hess_inv
    
            
params_result = ret.x
chisqr = ret.fun

nfree = len( data.flatten() ) - len( params_result ) 
redchisqr = chisqr / nfree
print( 'chisqr: ', chisqr ) 
print( 'redchisqr: ', redchisqr ) 

params_result_errors = np.sqrt( np.diag( cov ) * redchisqr )
print( 'result: ', params_result )
print( 'errors: ', params_result_errors ) 



print( data ) 




# ret = scipy.optimize.leastsq( objective, params_guess, ftol = 1e-12 )
# print( ret ) 


# try: 
#     cov = np.linalg.inv( hessian ) 
# except( np.linalg.LinAlgError ) :
#     self.cov = None 
#     self.success = 0
    
        
# params_result = ret.x
# chisqr = ret.fun

# nfree = len( data.flatten() ) - len( params_result ) 
# redchisqr = chisqr / nfree
# print( 'chisqr: ', chisqr ) 
# print( 'redchisqr: ', redchisqr ) 

# params_result_errors = np.sqrt( np.diag( cov ) * redchisqr )
# print( 'result: ', params_result )
# print( 'errors: ', params_result_errors ) 
