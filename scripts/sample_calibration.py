## includes 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import json

# my lib 
import jacob_math
import jacob_pyplot
import sql_db_manager
import deadlayer_analysis
import jacob_utils




# configurable test file for example()
test_dir = '../extracted_ttree_data/'
test_file = test_dir + '/deadlayerdet3rt/deadlayerdet3rt_16_16.bin'



# input filename of database, return DataFrame containing the DB  
def read_db_into_df( db_name ):
    with sqlite3.connect( db_name ) as conn:
        return pd.read_sql_query( 'SELECT * from ' + sql_db_manager.tablename, conn )
    return None



# input: dataframe in the format of the sql tables storing data
def get_mu_values( df, coords ):
    return get_values( df, coords, 'pf' )
            


def get_mu_delta_values( df, coords ):
    return get_values( df, coords, 'pferr' )


# function for retrieveing arrays of mu values or mu uncertainty values.
# currently cannot handle anything else because of dependence on get_fitnum_and_index.
# used to construct get_mu_values and get_mu_delta_values 
def get_values( df, coords, col_name ):
    # read in data: the 5 mu values for a particular fit
    values = []

    # for each peak: look up value and appned to mu_values either mu or np.nan 
    for peaknum in np.arange(5):
        fitnum, index_in_pf = deadlayer_analysis.get_fitnum_and_index_in_pf( peaknum )
        successful_fit = df.successful_fit.values[ coords[0]*32 + coords[1] + fitnum ]

        if successful_fit:
            pf = json.loads( df[col_name].loc[ coords[0]*32 + coords[1] + fitnum ] ) 
            values.append( pf[ index_in_pf ] )
        else:
            values.append( np.nan )

    return np.asarray(values)



# do an example on the test data 
def sample_calibration():

    # need db / coords for which every fit converged.
    db_name = sql_db_manager.centered_db
    df = read_db_into_df( db_name )
    coords = (15,15) 
    mu_values = get_mu_values( df, coords )
    mu_delta_values = get_mu_delta_values( df, coords )

    # this is a constant array of the peak energies, we have to flattne though. 
    energies = np.asarray( jacob_utils.flatten_list( deadlayer_analysis.peak_energies ) )

    # check if nan 
    if any( np.isnan( mu_values ) ):
        print 'ERROR: one of the fits did not converge, use different pixel.'
        return 0

    # print energies
    # print mu_values
    # print mu_delta_values

    plt.clf()
    ax = plt.axes()
    ax.errorbar( x=energies, y=mu_values, yerr=mu_delta_values, fmt='None', errorevery=1 )
    plt.show()
    
    p0 = [ 0.5, -118.0 ]

    # moment of truth 
    linear_fit = lambda p, x_: p[0]*x_ + p[1]

    result = jacob_math.jacob_least_squares( energies, mu_values, mu_delta_values, p0, linear_fit )
    if result is not None:
        reduc_chisq, dof, pf, pferr = result
    else:
        print 'ERROR: fit failed to converge.'
        return 0


    # print the data
    messages = [ 'reduc_chisq', 'dof', 'pf', 'pferr' ]
    for i in range(len(result)):
        print messages[i] + ': ' + str(result[i])

  
    # # plot the fit
    # fit_bounds = [ min(mu_values), max(mu_values) ]
    # x_forfit = np.linspace( fit_bounds[0], fit_bounds[1], 200*mu_values.size )
    # jacob_pyplot.add_fit_to_plot( ax, x_forfit, fit_bounds, pf, pferr, linear_fit )

    
    
    # ax.plot( x_forfit, linear_fit( pf, x ), 

    jacob_pyplot.set_linear_scale_plot_bounds( ax, energies, mu_values )
    jacob_pyplot.add_legend(ax) 



sample_calibration()



    
