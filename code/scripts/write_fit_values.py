import matplotlib
matplotlib.use( 'agg' ) 


import exp1.exp1_geometry
import exp2.exp2_geometry 

import jspectroscopy as spec
import numpy as np 



# db_names = [ 'moved' ]
#db_names = [ 


# db_names = [ 'det3_cent', 'det3_moved' ] 


# db_names = [ 'alpharun11-19', 'alpharun20-30', 'full_bkgd_tot',
#              'angled', 'moved', 'centered', 'flat', 'det3_cent', 'det3_moved' ] 


# geometry
exp1_secant_matrices = exp1.exp1_geometry.get_secant_matrices()
exp2_secant_matrices = exp2.exp2_geometry.get_secant_matrices()

print( exp2_secant_matrices[0][3][0] )


# actual_energies = [ [ 3182.690 ], [ 5762.64, 5804.77 ] ]
peak_indices = [ [1], [0,1] ]

# db_names = [ 'full_bkgd_tot' ]
# source_names = ['Gd 148', 'Cm 244']




db_names = [ 'centered', 'angled', 'flat', 'moved' ]

peak_indices = [ [1,2], [1,2], [1] ]
source_names = [ 'Pu 238', 'Pu 240', 'Cf 249' ]




for name in db_names :

    if name in [ 'centered', 'angled', 'flat', 'moved' ] :
        secant_matrices = exp1_secant_matrices[ name ] 
        
    else :
        secant_matrices = exp2_secant_matrices[1:]

        
    db = spec.spectrum_db( name, '../../storage/' ) 

    # db.write_peak_values( 1 )
    db.write_all_fit_params()
    db.package_all_fit_params()
    # db.calibrate_pixels( peak_indices, actual_energies )
    # db.write_calibrated_params()
    # db.write_peakdetect()

    db.save_dill( secant_matrices, 'secant_matrices' )
    # db.plot_heatmap( 'secant_matrices', source_names  )
    db.plot_all_params( source_names, secant_matrices  ) 

    db.disconnect() 

    # mu_values = test_db.read_mu_values( path )

    # print( mu_values ) 
