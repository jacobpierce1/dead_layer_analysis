import numpy as np
import deadlayer_delpes.deadlayer_estimator as dl_estimator


    




    

# dbs = [ dbmgr.angled ]
# dbs = [ dbmgr.centered, dbmgr.angled, dbmgr.moved, dbmgr.flat ] 
# dbs = [ dbmgr.flat ]
# dbs = dbmgr.all_dbs 
# dbs = [ dbmgr.angled ]

# dbs = dbmgr.det3_dbs
dbs = [ dbmgr.det3_moved ] 

# source_indices = [ [0,1], [0,1], [1] ]
source_indices = [ [0,1], [0,1], [1] ]



subtitle = ''

fstrips = np.arange( 2, 30 )

# print( fstrips )
# fix_source_deadlayers = { 'source_constant_0_0' :   13.1292106,
#                           'source_constant_0_1' :   12.5279746,
#                           'source_constant_1_0' :   3.90229562,
#                           'source_constant_1_1' :   2.89311618,
#                           'source_constant_2_1' :   5.11531913 }

fix_source_deadlayers = { 'source_constant_0' :   8.52309058,
                          'source_constant_1' :   3.34767770,
                          'source_constant_2' :   8.91838293  }


model_params = dead_layer_model_params( vary_det_deadlayer = 1,
                                        quadratic_source = 0,
                                        quadratic_det = 0,
                                        calibrate_each_pixel = 0,
                                        interp_stopping_power = 1,
                                        mu = 0,
                                        average_over_source = 0,
                                        pulse_height_defect = 0,
                                        fstrips_requested = fstrips,
                                        bstrips = np.arange( 2, 30 ),
                                        different_fstrip_deadlayers = 0,
                                        different_pixel_deadlayers = 0,
                                        fix_source_deadlayers = fix_source_deadlayers,
                                        one_source_constant = 1 )
                                        # ignore_outer_pixels = 0 )



linear_calibration_on_each_x_strip( dbs, source_indices,
                                    model_params,
                                    cut_high_sectheta = 0,
                                    annotate = 0,
                                    # view_pixel = [ dbmgr.moved, 5 ],
                                    subtitle = 'Det 3: Source Params Extrapolated from Det 1\n',
                                    reset_angles = None,
                                    residual_scatter_plot = 0,
                                    plot_3d = 1 ) #  ,
                                    # savefig_dir = '../../../deadlayer_analysis_paper/images/' )



# secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1,
#                                                 average_over_source = 0 )

# print( secant_matrices[ 'cf_249' ][1].x ) 

# secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1,
#                                                 average_over_source = 1 )
# print( secant_matrices[ 'cf_249' ][1].x ) 



