# this script is made for managing the DB, including creating the initial schema,
# adding columns later for data analysis, wrappers for reading / writing to the 
# schema, and 
# 
# source: https://pymotw.com/3/sqlite3/

import os
import sqlite3
import json
import pandas as pd
import numpy as np
import libjacob.jmeas as meas

from libjacob.jutils import isint
import sys
import jspectroscopy as spec 

import dill # allow lambda function pickle
# import _pickle



# config 
DEBUG_DB = 0









################################
# begin module #################
################################


_current_abs_path = os.path.dirname( __file__ ) + '/'



schema_filename = _current_abs_path + 'initial_schema.sql'
_tablename = 'fits'


schema_cols = [ 'x', 'y', 'fit_id', 'successful_fit',
                'npeaks', 'last_attempt', 
                'params_guess', 'fit_bounds',
                'peak_guesses', 'model']



# class providing complete set of operations for writing
# and specialized querying for the DBs.

class db( object ):

    
    def __init__( self, name, dimensions = None, features_shape = None ):

        # state the number of x and y dimensions of the detector
        # to be stored in the database. for example, if you passed
        # (32, 32) it will generate a sql database storing the
        # info of a 32 x 32 grid 

        if dimensions is None:
            self.dimensions = ( 1, 1 )

        self.xdim = dimensions[0]
        self.ydim = dimensions[1] 
        self.dimensions = dimensions
        
        
        # must be 1D list, tuple, or ndarray such that the
        # number of entries is equal to the number of fits
        # that will be performed on each pixel, and each entry
        # is the number of peaks that will be searched for on
        # the corresponding feature.
        #
        #TODO: only works for alpha peaks. add feature to specify
        # different types of peaks at each feature.
        
        if features_shape is None:
            self.features_shape = ( 1, )

        self.features_shape = features_shape
        
        self.name = name

        self.path = ( _current_abs_path + '../../../fit_results/databases/'
                      + name + '_fits_data.db' ) 

        self.mu_vals_dir = _current_abs_path + '../../../fit_results/mu_values/' 
        

        self.peak_vals_dir = _current_abs_path + '../../../fit_results/peak_values/'

        for dir_ in [ self.mu_vals_dir, self.peak_vals_dir ] :
            if not os.path.exists( dir_ ) :
                os.mkdir( dir_ )
                
                             
        self.conn = None
        
        self.num_peaks_per_feature = ( 2, 2, 2 )

        self.num_features = len( self.num_peaks_per_feature ) 

        # self.alternate_shape = [ 0, 0, [0,1] ]
        
        # self.num_sources = len( num_peaks_per_fit ) 
        

        
        
    # check if this db has been created. 
    def exists( self ):
        return os.path.exists( self.path ) 

    

    
    def connect( self, empty = 0 ):

        if self.conn is not None:
            print( 'ERROR: db is already open.' )
            return 0

        # throw an error for an attempt to connect to a nonexistent
        # database, unless we insist that it is supposed to be empty.
        
        if not empty:
            if not self.exists():
                print( 'ERROR: db has not been created.' )
                return 0 
            
        self.conn = sqlite3.connect( self.path ) 

        # this is set to allow you to read into an sqlite3.Row,
        # which is a dict-like object. 
        self.conn.row_factory = sqlite3.Row

        return 1

    

    
    def disconnect( self ):

        if self.conn is None:
            return 0 
        
        self.conn.close()
        self.conn = None
        
        return 1

    
    
    
    # call this before doing any reading or writing.
    def assert_open( self ):
        if self.conn is None:
            raise ValueError( 'db is not open.' )

    



    def flatten( self ):
        return _meas_no_checks( self.x.flatten, self.dx.flatten() )
        
        

    # this function shall only write the fit parametetrs and things
    # that can be obtained only from the histograms, especially the
    # fits. the fwhms etc. can be extracted later much more rapidly by
    # reconstructing the fit function.  p0 and fit_bounds also saved
    # because they could be used to reconstruct the fit.  the default
    # behavior is to overwrite data, this is to keep the function as
    # efficient as possible. there several instances i can think of in
    # which you would want to put data in without bothering to check
    # the current state of the DB. as such you have to do such a query
    # beforehand if necessary.
    
    def insert_fit_data( self, x, y, fit_id,
                         successful_fit = 0,
                         npeaks = -1,
                         last_attempt = -1, 
                         params_guess = None,
                         fit_bounds=None,
                         peak_guesses = None,
                         model = None,
                         db_is_empty = 0 ):        

        
        if self.conn is None:
            raise ValueError( 'Cannot insert data, sqlite3 connection is not open. Call db.connect().' )

        
        # if db is empty, then we are doing an insert.
        if db_is_empty:
            query = ( 'insert into ' + _tablename + ' ('   
                      + ', '.join(schema_cols)    
                      + ') values ('         
                      + ':' + ', :'.join(schema_cols) 
                      + ');'
            )
            
        # otherwise we are doing an update.
        else:
            query = ( 'update ' + _tablename   
                      + ' set ' + ', '.join( [ col + '=:' + col
                                               for col in schema_cols[3:] ] )   
                      + ' where ' + ' and '.join( [ col + '=:' + col
                                                    for col in schema_cols[:3] ] ) )  
                
        query_dict =  {
            'x' : x,
            'y' : y,
            'fit_id' : fit_id,
            'successful_fit' : successful_fit, 
            'npeaks' : npeaks,
            'last_attempt' : last_attempt,
            'params_guess' : _to_bin( params_guess ),
            'fit_bounds' : _to_bin( fit_bounds ),
            'peak_guesses' : _to_bin( peak_guesses ),
            'model' : _to_bin( model )
        }
        
        self.conn.cursor().execute( query, query_dict )
        self.conn.commit()


        


    def read_fit_data( self, x, y, fit_id ):

        # print( (x,y,fit_id) )

        if self.conn is None:
            raise ValueError( 'Cannot read db, sqlite3 connection is not open. Call db.connect().' )
                
        query = ( 'select ' + ', '.join( schema_cols[3:] )   
                  + ' from ' + _tablename    
                  + ' where ' + ' and '.join( [ col + '=:' + col
                                                for col in schema_cols[:3] ] ) )
            
        query_dict = {
            'x' : x,
            'y' : y,
            'fit_id' : fit_id
        }

        cursor = self.conn.cursor()

        cursor.execute( query, query_dict )
        
        result = dict( cursor.fetchone() )

        # unpickle the pickled objects 
        for key in [ 'params_guess', 'fit_bounds', 'peak_guesses', 'model' ]:
            result[ key ] = _from_bin( result[ key ] )

        # print( 'result: ' + str( result ) )

        return result

    



    # only needs to be called once. has guard against future calls (unless the 
    # location of the db is changed, then it won't work)
    def create( self ):
        
        db_is_new = not os.path.exists( self.path )

        if db_is_new:
            print( 'INFO: creating DB and schema for ' + self.path + '...' )
            
        else:
            print( 'ERROR: database already exists, returning' )
            return 0
        
        # with sqlite3.connect( self.path ) as conn:

        self.connect( empty = 1 )
        
        with open( schema_filename, 'rt' ) as f:
            schema = f.read()
            self.conn.executescript( schema )    

        print( 'INFO: success, now populating...' )
            
        self.populate( 32, 32, 3 )

        self.disconnect()
            


            
    # fill the db with each fit id that we need, giving the needs update flag for
    # everything. this is meant to only be called when the table is empty 
    def populate( self, numx, numy, numfits ):

            for y in range(numx):
                for x in range(numy):
                    for i in range(numfits):
                        self.insert_fit_data( x, y, i, db_is_empty = 1 ) 
                        




    def delete( self ):
        
        if os.path.exists( self.path ):
            
            ans = raw_input( 'PROMPT: delete ' + filename + ', are you sure (y/n) ?  ' )
            
            if( ans == 'y' ):

                os.remove( self.path ) 
                print( 'INFO: deleted db.' )         
                return 1

            print( 'INFO: did not delete db.' )

            return 0

        return 1
    






    
    def get_mu_for_x_strip( self, x ):

        if not isint( x ):
            raise( ValueError( 'x is not an int.' ) )

        # print( 'INFO: loading mu values for x strip ' + str( x ) + ' ...' )

        disconnect_conn_when_done = 0
        if self.conn is None:
            self.connect()
            disconnect_conn_when_done = 1

        # populate the grid 
        strip = np.empty( (3, 2), dtype = 'object' )
        for i in range(3):
            for j in range(2):
                strip[i][j] =  meas.meas.empty( 32 )
        
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename + ' WHERE x = (?)', (x,) )

        for row in cursor:

            # if fit didn't converge then all the mu values for that
            # feature are assigned np.nan
            
            if not row[ 'successful_fit' ]:
                for i in range( 2 ):
                    # grid[ row[ 'fit_id' ][ i ] ].x[ row['y'] ] = np.nan
                    # grid[ row[ 'fit_id' ][ i ] ].dx[ row['y'] ] = np.nan
                    strip[ row[ 'fit_id' ], i ][ row['y'] ] = meas.nan
                    
                continue

            mu = spec.get_alpha_params_dict( _from_bin( row['model'] ) )[ 'mu' ] 
            if mu.size() == 1:
                # mu = meas.meas( [ np.nan, mu.x[0] ], [ np.nan, mu.dx[0] ] )
                mu = meas.append( meas.nan, mu )
                
            for i in range(2):
                strip[ row[ 'fit_id' ], i ][ row['y'] ] = mu[i]

                
        if disconnect_conn_when_done:
            self.disconnect()
                
        return strip 


    

                               

    def get_mu_for_y_strip( self, y ):
        raise( NotImplementedError( 'mu is not const in y direction' ) )


    



    
    # get the array of mu values for each pixel and for each peak.
    
    def get_all_mu_grids( self, read_from_file = 0 ) :

        # construct appropriately sized container for all the mu grids
        mu_grids = [ [ 0 ] * self.num_peaks_per_feature[i]
                     for i in range( self.num_features ) ]

        for i in range( self.num_features ):
            mu_grids[i] = [ meas.meas.empty( self.dimensions ) ] * self.num_peaks_per_feature[i] 


        # meas.meas.empty( self.num_peaks_per_fit + self.dimensions )

        if read_from_file :

            mu_paths = [ [ [  self.mu_vals_dir + self.name + '_%d_%d_%s.bin' % (i,j,s)
                              for s in [ 'x', 'dx' ] ]
                           for j in range( self.num_peaks_per_feature[i] ) ]
                         for i in range( self.num_features) ]
                        
            if all( [ os.path.exists( path ) for path in np.array( mu_paths ).flatten() ] ) :

                for i in range( self.num_features ) :
                    for j in range( self.num_peaks_per_feature[i] ) :

                        mu = np.fromfile( mu_paths[i][j][0] ).reshape( self.dimensions )
                        mu_delta = np.fromfile( mu_paths[i][j][1] ).reshape( self.dimensions )
                        mu_grids[i][j] = meas.meas( mu, mu_delta )
                        
                return mu_grids

            else:
                print( 'INFO: requested to read mu and mu_delta from files, but they aren\'t there. constructing them now...' )


                

        # construct the array from the db if the files don't exist
        # yet, or we requested a fresh fetch.

        disconnect_conn_when_done = 0
        if self.conn is None:
            self.connect()
            disconnect_conn_when_done = 1

            
        # populate the grid 
        for i in range(3):
            for j in range(2):
                mu_grids[i][j] =  meas.meas.empty( self.dimensions )
        
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        for row in cursor:

            x = row['x']
            y = row['y'] 
            feature = row[ 'fit_id' ]
            
            # if fit didn't converge then all the mu values for that
            # feature are assigned np.nan            
            if not row[ 'successful_fit' ]: 

                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    mu_grids[ feature ][ i ][ x,y ] = meas.nan
                continue

            mu = spec.get_alpha_params_dict( _from_bin( row['model'] ) )[ 'mu' ] 

            # check if the fit used the alternate shape 
            if len( mu ) != self.num_peaks_per_feature[ feature ] :

                # todo: this doesn't work in the general case.
                # this should depend on the specified alternate shape.
                mu = meas.append( meas.nan, mu ) 


            # populate the corresponding entry of the grid.
            for i in range( self.num_peaks_per_feature[ feature ] ):
                mu_grids[ feature ][ i ][ x,y ] = mu[i]

                
        if disconnect_conn_when_done:
            self.disconnect()


        # write to a file if requested. note that we have already constructed
        # all the file names.
        if read_from_file :
            
            for i in range( self.num_features ) :
                for j in range( self.num_peaks_per_feature[i] ) :
                    
                    mu_grids[i][j].x.tofile(  mu_paths[i][j][0] )
                    mu_grids[i][j].dx.tofile( mu_paths[i][j][1] )
                    
        
        return mu_grids
            

        
#         for x in range( 32 ) :


            
            # mu = self.get_mu_for_x_strip( x )


            







                    
    # # if populating an entire array, we can get a bit more efficiency by not calling
    # # get_mu_values
    # def get_mu_grids_where_valid( self ):

    #     print( 'INFO: loading mu grid...' )


    #     disconnect_conn_when_done = 0
    #     if self.conn is None:
    #         self.connect()
    #         disconnect_conn_when_done = 1
            

    #     # to be combined in a meas.
    #     grid = np.empty( (3, 2), dtype = 'object' )
    #     for i in range(3):
    #         for j in range(2):

    #             grid[i][j] = meas.meas( np.empty( ( 32, 32 ) ),
    #                                np.empty( ( 32, 32 ) ) )
        
    #     cursor = self.conn.cursor()
    #     cursor.execute( 'SELECT * FROM ' + _tablename )

    #     for row in cursor:

    #         # if fit didn't converge then all the mu values for that
    #         # feature are assigned np.nan
            
    #         if not row[ 'successful_fit' ]:
    #             for i in range( 2 ):
    #                 grid[ row[ 'fit_id' ] ][ i ][ row['x'], row['y'] ] = meas.nan
    #             continue
            
    #             # vals[ row[ 'fit_id' ], :, row['x'], row['y'] ] = np.nan
    #             # deltas[ row[ 'fit_id' ], :, row['x'], row['y'] ] = np.nan
    #             continue

    #         mu = spec.get_alpha_params_dict( _from_bin( row['model'] ) )[ 'mu' ] 
    #         if mu.size() == 1:
    #             mu = meas.meas( [ np.nan, mu.x[0] ], [ np.nan, mu.dx[0] ] )

    #         for i in range(2):
    #             # print( ( row['x'], row['y'], row['fit_id'] ), i )
    #             # print( mu[i] ) 
    #             grid[ row[ 'fit_id' ] ][ i ][ row['x'], row['y'] ] = mu[i]

                
    #     if disconnect_conn_when_done:
    #         self.disconnect()
                
    #     return grid 

        

    



    
# return binary dump that can be
# inserted in the DB as a BLOB
def _to_bin( data ): 
    return sqlite3.Binary( dill.dumps( data, protocol = 2 ) )




# read binary as written in the _to_bin format
def _from_bin( bin ):
    return dill.loads( bin )
    
    

###########################################
## DECLARE GLOBALS FOR USE ELSEWHERE ######
###########################################

all_dbs = [ db( name, (32,32), (2,2,2) )
            for name in [ 'centered', 'moved', 'flat', 'angled' ] ]
centered, moved, flat, angled = all_dbs

normal_dbs = [ centered, moved, flat ]


# # fitnum is the fit that this peak belongs to and index_in_pf is 
# # the index of the peak within pf. a dedicated function is necessary because of 
# # the convoluted way in which the data is stored.

# def get_fitnum_and_mu_index_in_pf( peaknum ):
#     return ( peaknum // 2, 5 + 2* (peaknum % 2) )


    
