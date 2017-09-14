# my includes 
import sql_db_manager
import jacob_file_parsers

## includes 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import json


# this script is meant to be run after parse_all_data.py. that script writes all the data
# about all the fits (including failed ones) to a json file. that file is then read 
# here and analyzed. there is no data processing performed in the previous script 
# other than fitting the functions. in this script we reconstruct the fit functinos 
# using the fit parameters in the json file and pick up the analysis from there.
# this is because the curve fitting process takes a while but is final, the data
# analysis phase can take much longer. 

# modifiable config
FILE_SILICON_STOPPING_POWER = '../stopping_power_data/alpha_in_si_stopping_power.txt'


# constants, do not touch
NUM_FITS_PER_PIXEL = 3
NUM_PEAKS_PER_FIT = [ 2, 2, 1 ]



# here are the main peaks:
#240 PU:
#   http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=240PU&unc=nds	
#   5021.23 15       0.00422 15 
#   5123.68 23 	    27.10 % 10 	  1.389 5 
#   5168.17 15 	    72.80 % 10 	  3.762 5 

# 238 Pu: 
#   http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=238PU&unc=nds
#   5456.3 3 	    28.98 % 10 	  1.581 5 
#   5499.03 20 	    70.91 % 10 	  3.899 6 

# 249 Cf: 
#     5813.3 10 	    82.2 % 5 	  4.78 3 


# energies in keV of main alpha peaks of 240 Pu, 238 Pu, 249 Cf 
peak_energies = [ [ 5123.68, 5168.17 ], [ 5456.3, 5499.03 ], [5813.10] ]
peak_energies_delta = [ [0.23, 0.15], [0.3, 0.20], [0.10] ]

# https://physics.nist.gov/cgi-bin/Star/compos.pl?mode=text&refer=ap&matno=014 
density_silicon = 2.33000E+00  

# construct a pd.Series from the stopping power file.
# si_stopping_power = pd.Series( np.loadtxt( FILE_SILICON_STOPPING_POWER, skiprows=11, usecols=(0,3), unpack=1 ) )


def energy_lost_in_si_dead_layer( initial_energy, distance_traveled ):
    pass 

    
        
def analysis():
    
    # attempt to load our sql databases into a DataFrame    
    databases = [ sql_db_manager.rotated_db, sql_db_manager.centered_db ]
    dataframes = [0] * len(databases) 
    
    for i in range(len(databases)):
        with sqlite3.connect( databases[i] ) as sql_conn:
            # df = pd.read_sql_table( sql_db_manager.tablename, sql_conn )
            dataframes[i] = pd.read_sql_query( 'SELECT * from ' + sql_db_manager.tablename, sql_conn )
    # plot_mu_variations( dataframes[0] )
    
    # todo: regression on  each plot 
    
    plt.clf()
    
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    # X, Y, Z = axes3d.get_test_data(0.05)
    # plot_mu_differences_grid_3d( ax, 4, dataframes[0], dataframes[1] )
    # plot_mu_differences_grid_3d_with_contours( ax, 4, dataframes[0], dataframes[1] )
    ax = plt.axes()
    # plot_mu_variations( ax, 1, dataframes [0] )
    plot_mu_differences_grid( ax, 1, dataframes[1], dataframes[0] )
    
    plt.show()


# fitnum is the fit that this peak belongs to and index_in_pf is 
# the index of the peak within pf. a dedicated function is necessary because of 
# the convoluted way in which the data is stored. 
def get_fitnum_and_index_in_pf( peaknum ):
    return ( peaknum // 2, 5 + 2* peaknum % 2 )



def plot_mu_differences_grid_3d_with_contours( ax, peaknum, df1, df2 ):
    X, Y = np.meshgrid( range(32), range(32) )
    mu_differences_grid = get_mu_differences_grid( peaknum, df1, df2 )
    imshow_3d_with_projections( ax, X, Y, mu_differences_grid ) 
    
    
def plot_mu_differences_grid_3d( ax, peaknum, df1, df2 ):
    X, Y = np.meshgrid( range(32), range(32) )
    mu_differences_grid = get_mu_differences_grid( peaknum, df1, df2 )
    imshow_3d( ax, X, Y, mu_differences_grid )
    
    # im = ax.imshow( mu_differences_grid, cmap='plasma', interpolation='none' )
    # plt.colorbar(im)


## in this case, avearge over the 
#def plot_average_mu_differences( ax, peaknum, df1, df2 ):
#    pass


def plot_mu_differences_grid( ax, peaknum, df1, df2 ):
    mu_differences_grid = get_mu_differences_grid( peaknum, df1, df2 )
    im = ax.imshow( mu_differences_grid, cmap='plasma', interpolation='none' )
    plt.colorbar(im)
    


# get (32,32) array with values of mu(df1) - mu(df2) 
def get_mu_differences_grid( peaknum, df1, df2 ):
    df1_mu_grid = get_mu_grid_where_valid( peaknum, df1 )
    df2_mu_grid = get_mu_grid_where_valid( peaknum, df2 )
    return df1_mu_grid - df2_mu_grid 
    

    
    
def get_mu_grid_where_valid( peaknum, df ):
    fitnum, peaknum_in_pf = get_fitnum_and_index_in_pf( peaknum )
    successful_fit_values = df.successful_fit.values[ fitnum : : NUM_FITS_PER_PIXEL ]
    pf_values = df.pf.values [ fitnum : : NUM_FITS_PER_PIXEL ]
    return np.array( [ json.loads(pf_values[i])[ peaknum_in_pf ] if successful_fit_values[i] \
            else np.nan for i in range(pf_values.size) ] ).reshape(32,32) 
    
    
    ## populate 
    #mu_differences_grid = np.empty(32,32)
    #count = 0
    #for x in range(32):
    #    for y in range(32):
    #        if( mu1_successful_fit & mu2_successful_fit ) 
    #        mu_differences_grid[x][y] = 
    #    
#     mu1_mu_values =  np.where( mu1_successful_fit, , nan_array )  
    
    
    
def plot_mu_variations( ax, peaknum, all_data ):
    mu_grid = get_mu_grid_where_valid( peaknum, all_data )
    im = ax.imshow( mu_grid, cmap='plasma', interpolation='none' )
    plt.colorbar(im)
    # ax.xlabel('X')
    # ax.ylabel('Y')


def imshow_3d_with_projections(ax, X, Y, Z):

    ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
    cset = ax.contourf(X, Y, Z, zdir='z', offset=-100,
            levels=np.linspace(-100,100,1200),cmap=plt.cm.jet)
    cset = ax.contourf(X, Y, Z, zdir='x', offset=-40, cmap=plt.cm.jet)
    cset = ax.contourf(X, Y, Z, zdir='y', offset=40, cmap=plt.cm.jet)
    ax.set_xlabel('X')
    # ax.set_xlim(-40, 40)
    ax.set_ylabel('Y')
    # ax.set_ylim(-40, 40)
    ax.set_zlabel('Z')
    # ax.set_zlim(-100, 100)    



def imshow_3d( ax, X, Y, Z ):
    x = X.flatten()
    y = Y.flatten()
    z = Z.flatten()
    ax.bar3d( x, y, np.zeros(len(z)), 1, 1, z )
    ax.set_zlim( min(z), max(z) )
    # ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)



analysis()

    
#    
## input: number between 0 and 5 corresponding to the peak number in one of the plots
## output: contour plot of the mu value throughout the grid
#def get_mu_variations( peaknum, all_data ):
#    
#    # construct grids 
#    xcoords = range(32) 
#    ycoords = range(32)
#    xgrid, ygrid = np.meshgrid( xcoords, ycoords )
#    mu_grid = np.empty( (32,32) )
#
#    # determine the fit_id from specified peaknum (between 0 and 4) 
#    # this hack here works because of the fact that the num peaks in each fit is 
#    # [2,2,1]
#    fit_id = peaknum // 2 
#    
#    #if peaknum >= 0 and peaknum < 2:
#    #    fit_id = 0            
#    #elif peaknum >= 2 and peaknum < 4:
#    #    fit_id = 1
#    #elif peaknum == 4:
#    #    peaknum = 2
#    #else:
#    #    print 'ERROR: invalid peak number in plot_mu_variations.'
#    #    return 0
#    
#    # because of poor programming construction in the past we are stuck with this.
#    mu_index_in_pf = 5 + 2* peaknum % 2
#    
#    # construct appropriate subset of the data    
#    subset = all_data[ all_data.fit_id == fit_id ]
#    
#    
#    # populate mu_grid
#    for i in xcoords:
#        for j in ycoords:
#            
#            # extract the element with the right coords 
#            extracted = subset[ (subset.x==i) & (subset.y==j) ]
#
#            if extracted.successful_fit.iloc[0]:            
#                
#                # load the extracted data as an array, hack necessary because of our 
#                # choice to store the arrays as strings in the DB.
#                pf = json.loads( (extracted.pf).iloc[0] ) 
#                
#                # finally take out the right value and add to mu_grid.
#                mu_grid[i][j] = pf[ mu_index_in_pf ]
#                
#            else:
#                mu_grid[i][j] = np.nan
#   
#    #cp = plt.contourf( xgrid, ygrid, mu_grid )
#    #plt.colorbar(cp)
#    #plt.title('Filled Contours Plot')
#    #plt.xlabel('x (cm)')
#    #plt.ylabel('y (cm)')
#    return mu_grid
   
   
   


#
#   
#   
## function that allows you to extract either amplitude or mu from pf or pf_err, 
## use this instead of playing with the indexing. data_idx is the index of either 
## pf or pf_err and start_idx_in_arr is the starting index of the data to be 
## accessed within pf or pf_err
#def retrieve_values_from_fit_array( arr, results, data_idx, start_idx_in_arr ):
#    arr = []
#    for i in range(NUM_FITS_PER_PIXEL):
#        successful_fit = results[i][0]
#        if successful_fit:
#            arr.extend( (results[data_idx])[start_idx_in_arr : start_idx_in_arr + NUM_PEAKS_PER_FIT[i]*2 : 2 ] )
#        else: 
#            arr.extend( [np.nan] * NUM_PEAKS_PER_FIT[i] )
#
#
#
#def retrieve_mu_values( results, db_num, x, y, table ):
#    retrieve_values_from_fit_array( table[db_num, x, y], results, 3, 4 )
#    
#def retrieve_mu_err_values( results, db_num, x, y, table ):
#    retrieve_values_from_fit_array( table[db_num, x, y], results, 4, 4 )
#
#
#
## print relevant stats 
#def stats( all_data ):
#    pass 
#
##
##
### this function populates table, which must have dims (2,32,32,x), with data from the dbs.
### retrieval function must take args (sql_conn, x, y, table) and must construct an array 
### which will be put in the table. 
##def populate_table( table, retrieval_function ):
##
##    if table.shape != (2,32,32):
##        print 'ERROR: shape of table must be (2,32,32)'
##        return 0
##    
##    databases = [ sql_db_manager.rotated_db, sql_db_manager.centered_db ]
##    
##    # populate mu_table and mu_err_table
##    for db_num in range(2):
##        with sqlite3.connect( databases[db_num] ) as sql_conn:
##            for x in range(32):
##                for y in range(32):
##                    results = [ sql_db_manager.read_data_from_db( sql_conn, (x,y), i )
##                            for i in range(NUM_FITS_PER_PIXEL)  ]
##                    retrieval_function( results, db_num, x, y, table )
##                        # successful_fit, fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data = result
##    return 1        
##
##
##def analysis():
##    # 4d arrays to store all the mu values 
##    mu_table = np.empty( (2,32,32) )
##    mu_err_table = np.empty( (2,32,32) )
##    
##    populate_table( mu_table, retrieve_mu_values )
##    populate_table( mu_err_table, retrieve_mu_err_values )
##    


