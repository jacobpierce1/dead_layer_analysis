# this script is made for managing the DB, including creating the initial schema,
# adding columns later for data analysis, wrappers for reading / writing to the 
# schema, and 
# 
# source: https://pymotw.com/3/sqlite3/

import os
import sqlite3
import json
import numpy as np

from libjacob.jutils import isint, estimate_time_left
import sys
import jspectroscopy as spec 

import datetime

import dill # allow lambda function pickle




# config 
DEBUG_DB = 0



# _current_abs_path = os.path.dirname( __file__ ) + '/'



schema_filename = 'initial_schema.sql'
_tablename = 'fits'


schema_cols = [ 'x', 'y', 'fit_id', 'successful_fit',
                'npeaks', 'last_attempt', 
                'params_guess', 'fit_bounds',
                'peak_guesses', 'model']




# class providing complete set of operations for writing
# and specialized querying for the DBs.
# det_number is optional number of the detector, in case
# multiple detectors were used for the same dataset.
# sources is a list of the names of the sources used.

class spectrum_db( object ):

    
    def __init__( self, path, peak_types = None, dimensions = None ) :

        create_new_db = ( peak_types is not None ) and ( dimensions is not None )
        
        self.path = path
        self.conn = None 

        if not os.path.exists( path ) :
            if create_new_db :
                self.is_empty = 1 
                self.create( peak_types, dimensions )
            else :
                print( 'ERROR: attempted to load db, but the path does not exist' ) 

        else : 
            self.is_empty = 0
            self.connect()
            self.read_metadata()
            
        
        
        return

    
        
    def create( self, peak_types, dimensions ):
        
        self.peak_types = peak_types
        self.dimensions = dimensions

        if not os.path.exists( self.path ) :
            print( 'INFO: creating DB and schema for ' + self.path + '...' )
            
        else:
            print( 'ERROR: database already exists, returning' )
            return 0

        self.is_empty = 1
        self.conn = None
        
        self.xdim = dimensions[0]
        self.ydim = dimensions[1] 
        self.dimensions = dimensions
                                     
        self.num_peaks_per_group = [ len(x) for x in peak_types ] 

        self.num_groups = len( self.num_peaks_per_group )
                
        # with sqlite3.connect( self.path ) as conn:

        self.connect()
        
        with open( os.path.dirname( __file__ ) + '/'
                   + schema_filename, 'rt' ) as f:

            schema = f.read()
            self.conn.executescript( schema )    

        print( 'INFO: success, now populating...' )
            
        self.populate( *dimensions, len( peak_types ) )        
        self.write_metadata()

        self.is_empty = 0
        

    

            
    # fill the db with each fit id that we need, giving the needs update flag for
    # everything. this is meant to only be called when the table is empty 
    def populate( self, numx, numy, numfits ):

            for x in range( self.xdim ):
                for y in range( self.ydim ):
                    for i in range( self.num_groups ):
                        self.insert_fit_data( x, y, i ) 



                        

    def write_metadata( self ) :
        
        time = str( datetime.datetime.now() )

        if self.is_empty :
            self.conn.cursor().execute( 'INSERT INTO metadata VALUES ( ?, ?, ?, ?, ? ) ',
                                        ( self.xdim, self.ydim,
                                          _to_bin( self.peak_types ),
                                          _to_bin( None ),
                                          time ) )

            self.conn.commit() 

            

            
    def read_metadata( self ) :

        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM metadata' )
        # xdim, ydim, peak_types, yield_data, timestamp = cursor.fetchall()
        metadata  = dict( cursor.fetchone() )
                
        self.xdim = metadata[ 'xdim' ] 
        self.ydim = metadata[ 'ydim' ] 
        self.peak_types = metadata[ 'peak_types' ]
        self.yield_data = metadata[ 'yield_data' ]
        self.timestamp = metadata[ 'timestamp' ]

        

        
        
    # check if this db has been created. 
    def exists( self ):
        return os.path.exists( self.path ) 

    

    
    def connect( self ):

        if self.conn is not None:
            print( 'ERROR: db is already open.' )
            return 0

        # throw an error for an attempt to connect to a nonexistent
        # database, unless we insist that it is supposed to be empty.
        
        if not self.is_empty :
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

    


    # TODO: remove
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
                         model = None ):        

        
        if self.conn is None:
            raise ValueError( 'Cannot insert data, sqlite3 connection is not open. Call db.connect().' )

        
        # if db is empty, then we are doing an insert.
        if self.is_empty :
            
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

            params = spec.get_alpha_params_dict( _from_bin( row['model'] ) ) 

            if not spec.alpha_model_valid_params_predicate( params ) :
            
                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    mu_grids[ feature ][ i ][ x,y ] = meas.nan
                continue

            mu = params[ 'mu' ]

            if len( mu ) > 1 :
                if ( np.abs( mu[1].x - mu[0].x - 20 ) > 10
                     or any( mu.dx < 0.1 ) ):
                    
                    for i in range( self.num_peaks_per_feature[ feature ]  ):
                        mu_grids[ feature ][ i ][ x,y ] = meas.nan
                    continue
            
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



    




    
    

    def get_all_peak_grids( self, read_from_file = 1 ) :
        
        # construct appropriately sized container for all the peak grids
        peak_grids = [ [ 0 ] * self.num_peaks_per_feature[i]
                     for i in range( self.num_features ) ]

        for i in range( self.num_features ):
            peak_grids[i] = [ meas.meas.empty( self.dimensions ) ] * self.num_peaks_per_feature[i] 


        # meas.meas.empty( self.num_peaks_per_fit + self.dimensions )

        if read_from_file :

            peak_paths = [ [ [  self.peak_vals_dir + self.name + '_%d_%d_%s.bin' % (i,j,s)
                              for s in [ 'x', 'dx' ] ]
                           for j in range( self.num_peaks_per_feature[i] ) ]
                         for i in range( self.num_features) ]
                        
            if all( [ os.path.exists( path ) for path in np.array( peak_paths ).flatten() ] ) :

                for i in range( self.num_features ) :
                    for j in range( self.num_peaks_per_feature[i] ) :

                        peak = np.fromfile( peak_paths[i][j][0] ).reshape( self.dimensions )
                        peak_delta = np.fromfile( peak_paths[i][j][1] ).reshape( self.dimensions )
                        peak_grids[i][j] = meas.meas( peak, peak_delta )
                        
                return peak_grids

            else:
                print( 'INFO: requested to read peak and peak_delta from files, but they aren\'t there. constructing them now...' )


                

        # construct the array from the db if the files don't exist
        # yet, or we requested a fresh fetch.

        disconnect_conn_when_done = 0
        if self.conn is None:
            self.connect()
            disconnect_conn_when_done = 1

            
        # populate the grid 
        for i in range(3):
            for j in range(2):
                peak_grids[i][j] =  meas.meas.empty( self.dimensions )
        
        cursor = self.conn.cursor()
        cursor.execute( 'SELECT * FROM ' + _tablename )

        rownum = 0
        total_iterations = ( self.dimensions[0] * self.dimensions[1]
                             * self.num_features ) 
        
        for row in cursor:

            rownum += 1
            
            estimate_time_left( rownum, total_iterations, num_updates = 100 ) 

            x = row['x']
            y = row['y'] 
            feature = row[ 'fit_id' ]

            print( (x,y,feature)) 
            
            # if fit didn't converge then all the peak values for that
            # feature are assigned np.nan
            
            if not row[ 'successful_fit' ]: 
                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    peak_grids[ feature ][ i ][ x,y ] = meas.nan
                continue

            
            model = _from_bin( row['model'] )


            # throw out the data if we didn't estimate errorbars.
            if not model.errorbars :
                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    peak_grids[ feature ][ i ][ x,y ] = meas.nan
                continue

            # do check on chi2. TODO: this should be removed and incorporated
            # in the fit itself.
            if 1 - chi2.cdf( model.chisqr, model.nfree ) < 0.05 :
                for i in range( self.num_peaks_per_feature[ feature ]  ):
                    peak_grids[ feature ][ i ][ x,y ] = meas.nan
                continue
                
            # do a monte carlo simulation for each peak in this feature.
                
            peaks = spec.estimate_alpha_peakpos( model, plot = 0 )

            # check if the fit used the alternate shape
            if len( peaks ) != self.num_peaks_per_feature[ feature ] :
                    
                # todo: this doesn't work in the general case.
                # this should depend on the specified alternate shape.
                peaks = meas.append( meas.nan, peaks ) 
                    
                    
            # populate the corresponding entry of the grid.
            for i in range( self.num_peaks_per_feature[ feature ] ):
                peak_grids[ feature ][ i ][ x,y ] = peaks[i]

        # reset the counter 
        estimate_time_left( 0, 0, reset = 1 ) 
                
        if disconnect_conn_when_done:
            self.disconnect()


        # write to a file if requested. note that we have already constructed
        # all the file names.
        if read_from_file :
            
            for i in range( self.num_features ) :
                for j in range( self.num_peaks_per_feature[i] ) :
                    
                    peak_grids[i][j].x.tofile(  peak_paths[i][j][0] )
                    peak_grids[i][j].dx.tofile( peak_paths[i][j][1] )
                    
        
        return peak_grids    
    





    



    
# return binary dump that can be
# inserted in the DB as a BLOB
def _to_bin( data ): 
    return sqlite3.Binary( dill.dumps( data, protocol = 2 ) )




# read binary as written in the _to_bin format
def _from_bin( bin ):
    return dill.loads( bin )
    

