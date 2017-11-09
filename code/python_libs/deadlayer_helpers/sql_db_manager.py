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





# constants, do not touch
NUM_FITS_PER_PIXEL = 3
NUM_PEAKS_PER_FIT = [ 2, 2, 1 ]
NUM_SOURCES = 3



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

        self.path = ( _current_abs_path + '../../../databases/'
                      + name + '_fits_data.db' ) 

        self.conn = None
        

        
        
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
    


    ############################################
    ## FUNCTIONS FOR SPECIALIZED QUERYING ######
    ############################################
    


    # # read all data for a particular pixel into a df.
    # def read_pixel_into_df( self, x, y ):

    #     self.assert_open()
        
    #     # check that x and y are ints to avoid injection attack
    #     if not ( isint(x) and isint(y) ):

    #         # print( type(x ) )
    #         # print( type( y ) ) 
    #         raise ValueError( 'x or y is not an int: x=%s, y=%s' 
    #                           % ( str(x), str(y) ) )
        
        
    #     return pd.read_sql_query( 'SELECT * from ' + _tablename
    #                               + ' WHERE x=%d AND y=%d ORDER BY fit_id'
    #                               % ( x, y ), self.conn )
    

    

    # # this function returns a list of 5 entries. each entry gives the
    # # fit parameters for a SINGLE alpha peak. note that this is not
    # # the same as the peak parameters obtained since we actually only
    # # do 3 fits. this is meant to be used for processing data from the
    # # single peak parameters, such as FWHM coords or peak positions.
    # def get_single_peak_fit_parameters( self, x, y ):

    #     self.assert_open()

    #     # to be returned
    #     pf_arr = []
    #     pf_delta_arr = []
        
    #     df = self.read_pixel_into_df( x, y )
        
    #     # print( df.shape )
        
    #     # print df 
        
    #     for peaknum in np.arange( 6 ):
            
    #         # these will accumulate the fit params and be appended to
    #         # pf_arr and pf_delta_arr
    #         pf_current = []
    #         pf_delta_current = []
            
    #         fitnum, mu_index_in_pf = get_fitnum_and_mu_index_in_pf( peaknum )
            
    #         row =  fitnum 
            
    #         successful_fit = df.successful_fit.values[ row ]
            
    #         if successful_fit:
                
    #             # read from db
    #             pf_from_db = json.loads( df['pf'].loc[row] )
    #             pf_delta_from_db = json.loads( df['pferr'].loc[row] ) 
                
    #             # add in the 2 tau values, sigma, and eta value
    #             pf_current.extend( pf_from_db[ 0:4 ] )
    #             pf_delta_current.extend( pf_delta_from_db[ 0:4 ] )
                
    #             # add in the A and mu values
    #             A = pf_from_db[ mu_index_in_pf - 1 ]
    #             mu = pf_from_db[ mu_index_in_pf ]
    #             A_delta = pf_delta_from_db[ mu_index_in_pf - 1 ]
    #             mu_delta = pf_delta_from_db[ mu_index_in_pf ]
                
    #             # extend current arrays 
    #             pf_current.extend( [A, mu] )
    #             pf_delta_current.extend( [A_delta, mu_delta ] )
                
    #             # add to pf_arr
    #             pf_arr.append( pf_current )
    #             pf_delta_arr.append( pf_delta_current )
                
                
    #             # otherwise we can't extract single-fit params for both
    #             # peaks.
                
    #         else:
    #             pf_arr.append( [ np.nan ] * 6  )
    #             pf_delta_arr.append( [ np.nan ] * 6 ) 
                
    #     # print( pf_arr )
    #     # print( pf_delta_arr ) 
                
    #     # return the values and deltas as a measurement.
    #     return meas.meas( pf_arr, pf_delta_arr ) 
            




    # # input filename of database, return DataFrame containing the DB  
    # def read_db_into_df( self ):
    #     self.assert_open()
    #     return pd.read_sql_query( 'SELECT * from ' + _tablename, self.conn )
    


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
    

               
        
    # if populating an entire array, we can get a bit more efficiency by not calling
    # get_mu_values
    def get_mu_grids_where_valid( self ):

        print( 'INFO: loading mu grid...' )


        disconnect_conn_when_done = 0
        if self.conn is None:
            self.connect()
            disconnect_conn_when_done = 1
            

        # to be combined in a meas.
        grid = np.empty( (3, 2), dtype = 'object' )
        for i in range(3):
            for j in range(2):

                grid[i][j] = meas.meas( np.empty( ( 32, 32 ) ),
                                   np.empty( ( 32, 32 ) ) )
        
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        for row in cursor:

            # if fit didn't converge then all the mu values for that
            # feature are assigned np.nan
            
            if not row[ 'successful_fit' ]:
                for i in range( 2 ):
                    grid[ row[ 'fit_id' ] ][ i ][ row['x'], row['y'] ] = meas.nan
                continue
            
                # vals[ row[ 'fit_id' ], :, row['x'], row['y'] ] = np.nan
                # deltas[ row[ 'fit_id' ], :, row['x'], row['y'] ] = np.nan
                continue

            mu = spec.get_alpha_params_dict( _from_bin( row['model'] ) )[ 'mu' ] 
            if mu.size() == 1:
                mu = meas.meas( [ np.nan, mu.x[0] ], [ np.nan, mu.dx[0] ] )

            for i in range(2):
                # print( ( row['x'], row['y'], row['fit_id'] ), i )
                # print( mu[i] ) 
                grid[ row[ 'fit_id' ] ][ i ][ row['x'], row['y'] ] = mu[i]

                
        if disconnect_conn_when_done:
            self.disconnect()
                
        return grid 

        

    



    
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



# # fitnum is the fit that this peak belongs to and index_in_pf is 
# # the index of the peak within pf. a dedicated function is necessary because of 
# # the convoluted way in which the data is stored.

# def get_fitnum_and_mu_index_in_pf( peaknum ):
#     return ( peaknum // 2, 5 + 2* (peaknum % 2) )


    
