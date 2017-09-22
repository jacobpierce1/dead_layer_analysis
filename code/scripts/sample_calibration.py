## includes 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# my lib 
import libjacob.jacob_math as jmath
import libjacob.jacob_pyplot as jplt
import libjacob.jacob_utils as jutils 
import libjacob.jacob_stats as jstats
import deadlayer_analysis as dl

import sql_db_manager
import deadlayer_analysis




# configurable test file for example
test_dir = '../data/extracted_ttree_data/'
test_file = test_dir + '/deadlayerdet3rt/deadlayerdet3rt_16_16.bin'





# do an example on the test data 
def sample_calibration():
        
    # need db / coords for which every fit converged.
    db_name = sql_db_manager.centered_db
    df = dl.read_db_into_df( db_name )
    coords = (15,15) 
    mu_values = dl.get_mu_values( df, coords )
    mu_delta_values = dl.get_mu_delta_values( df, coords )

    # check if nan 
    if any( np.isnan( mu_values ) ):
        print 'ERROR: one of the fits did not converge, use different pixel.'
        return 0

    # this is a constant array of the peak energies, we have to flattne though. 
    energies = np.asarray( jutils.flatten_list( deadlayer_analysis.peak_energies ) )

    plt.clf()
    ax = plt.axes()
    
    p0 = [ 0.5, -118.0 ]
    
    jstats.linear_calibration( energies, mu_values, mu_delta_values, p0, print_fit_data=0, ax=ax )

    # customize 
    plt.show()       

    return 1



sample_calibration()



    
