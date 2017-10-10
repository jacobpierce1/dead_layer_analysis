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



# all_db_ids = [ 'centered', 'moved', 'flat', 'angled' ]


# # construct the database names. these can be called
# # from any module that includes this to access the dbs.

# all_dbs_list = [ _current_abs_path + '../../../databases/'
#                  + _name + '_fits_data.db'
#                  for _name in all_db_ids ]


# all_dbs_dict = dict( zip( all_db_ids,
#                      all_dbs_list ) ) 


# centered_db, moved_db, flat_db, angled_db = all_dbs_list
# centered_id, moved_id, flat_id, angled_id = all_db_ids




schema_filename = _current_abs_path + 'initial_schema.sql'
tablename = 'fits_and_extracted_data'

schema_cols = [ 'x', 'y', 'fit_id', 'successful_fit', 'fit_attempt', 'reduc_chisq', 
                    'pf', 'pferr', 'p0', 'fit_bounds', 'peak_detect']



# class providing complete set of operations for writing
# and specialized querying for the DBs.

class db( object ):

    def __init__( self, name ):

        self.name = name

        self.path = ( _current_abs_path + '../../../databases/'
                      + name + '_fits_data.db' ) 

        self.conn = None
        

    def connect( self ):

        if self.conn is not None:
            return 0
        
        self.conn = sqlite3.connect( self.path ) 
        return 1

    
    def disconnect( self ):

        if self.conn is not None:
            return 0 
        
        self.conn.close()
        self.conn = None
        return 1

    
    # call this before doing any reading or writing.
    def assert_open( self ):
        if self.conn is None:
            raise ValueError( 'db is not open.' )

        

    # this function shall only write the fit parametetrs and things that
    # can be obtained only from the histograms, especially the fits. the
    # fwhms etc. can be extracted later much more rapidly by
    # reconstructing the fit function.  p0 and fit_bounds also saved
    # because they could be used to reconstruct the fit.  the default
    # behavior is to overwrite data, this is to keep the function as
    # efficient as possible. there several instances i can think of in
    # which you would want to put data in without bothering to check the
    # current state of the DB. as such you have to do such a query
    # beforehand if necessary.
    
    def insert_fit_data( self, x, y, fit_id, successful_fit=0,
                                 fit_attempt=-1, reduc_chisq=-1, pf=None, pferr=None,
                                 p0=None, fit_bounds=None, peak_detect=None,
                                 db_is_empty=0 ):        

        if self.conn is None:
            raise ValueError( 'Cannot insert data, sqlite3 connection is not open. Call db.connect().' )

        
        # if db is empty, then we are doing an insert.
        if db_is_empty:
            query = 'insert into ' + tablename + ' ('   \
                    + ', '.join(schema_cols)    \
                    + ') values ('         \
                    + ':' + ', :'.join(schema_cols)   \
                    + ');'
            
        # otherwise we are doing an update.
        else:
            query = 'update ' + tablename   \
                    + ' set ' + ', '.join( [ col + '=:' + col for col in schema_cols[3:] ] )   \
                    + ' where ' + ' and '.join( [ col + '=:' + col for col in schema_cols[:3] ] )   
            
                            
        cursor = sql_conn.cursor()
                
        query_dict =    \
                        {
                            schema_cols[0]:pixel_coords[0],
                            schema_cols[1]:pixel_coords[1],
                            schema_cols[2]:fit_id,
                            schema_cols[3]:successful_fit, # by default say that the fit was unsuccessful
                            schema_cols[4]:fit_attempt,
                            schema_cols[5]:reduc_chisq,
                            schema_cols[6]:json.dumps(pf),
                            schema_cols[7]:json.dumps(pferr),
                            schema_cols[8]:json.dumps(p0),
                            schema_cols[9]:json.dumps(fit_bounds),
                            schema_cols[10]:json.dumps(peak_detect) #[0,[0,0],[0,0]]
                        }
        
        self.conn.cursor().execute( query, query_dict )
        self.conn.commit()




    def read_fit_data( self, pixel_coords, fit_id ):

        if self.conn is None:
            raise ValueError( 'Cannot read db, sqlite3 connection is not open. Call db.connect().' )
        
        cursor = self.conn.cursor()
        
        query = 'select ' + ', '.join( schema_cols[3:] )   \
                + ' from ' + tablename    \
                + ' where ' + ' and '.join( [ col + '=:' + col for col in schema_cols[:3] ] )   
        
            
        query_dict = {
            schema_cols[0]:pixel_coords[0],
            schema_cols[1]:pixel_coords[1],
            schema_cols[2]:fit_id
        }
        
        cursor.execute( query, query_dict )
        
        result = cursor.fetchone()
        
        successful_fit, fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, peak_detect = result
        
        pf = json.loads(pf)
        pferr = json.loads(pferr)
        p0 = json.loads(p0)
        fit_bounds = json.loads(fit_bounds)
        peak_detect = json.loads(peak_detect)
        
        return ( successful_fit, fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, peak_detect )
    




    # only needs to be called once. has guard against future calls (unless the 
    # location of the db is changed, then it won't work)
    def create( self ):
        
        db_is_new = not os.path.exists( self.path )

        if db_is_new:
            print('INFO: creating DB and schema for ' + self.path + '...')
            
        else:
            print('ERROR: database already exists, returning')
            return 0
        
        with sqlite3.connect( self.path ) as conn:

            with open( schema_filename, 'rt' ) as f:
                schema = f.read()
                conn.executescript( schema )    

            print( 'INFO: success, now populating...' )
                
            populate_db( self.conn, 32, 32, 3 )
            


            
    # fill the db with each fit id that we need, giving the needs update flag for
    # everything. this is meant to only be called when the table is empty 
    def populate( self, numx, numy, numfits ):

        with sqlite3.connect( self.path ) as conn:
            for i in range(numx):
                for j in range(numy):
                    for k in range(numfits):
                        insert_fit_data_into_db( conn, (i,j), k, db_is_empty = 1 ) 
                        




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
    


    # read all data for a particular pixel into a df.
    def read_pixel_into_df( self, x, y ):

        self.assert_open()
        
        # check that x and y are ints to avoid injection attack
        if not ( isint(x) and isint(y) ):

            # print( type(x ) )
            # print( type( y ) ) 
            raise ValueError( 'x or y is not an int: x=%s, y=%s' 
                              % ( str(x), str(y) ) )
        
        
        return pd.read_sql_query( 'SELECT * from ' + tablename
                                  + ' WHERE x=%d AND y=%d ORDER BY fit_id'
                                  % ( x, y ), self.conn )
    


    # this function returns a list of 5 entries. each entry gives the fit parameters for a
    # SINGLE alpha peak. note that this is not the same as the peak parameters obtained since
    # we actually only do 3 fits. this is meant to be used for processing data from the single
    # peak parameters, such as FWHM coords or peak positions.
    def get_single_peak_fit_parameters( self, x, y ):

        self.assert_open()

        # to be returned
        pf_arr = []
        pf_delta_arr = []
        
        df = self.read_pixel_into_df( x, y )
        
        # print( df.shape )
        
        # print df 
        
        for peaknum in np.arange( 6 ):
            
            # these will accumulate the fit params and be appended to pf_arr and pf_delta_arr
            pf_current = []
            pf_delta_current = []
            
            fitnum, mu_index_in_pf = get_fitnum_and_mu_index_in_pf( peaknum )
            
            row =  fitnum 
            
            successful_fit = df.successful_fit.values[ row ]
            
            if successful_fit:
                
                # read from db
                pf_from_db = json.loads( df['pf'].loc[row] )
                pf_delta_from_db = json.loads( df['pferr'].loc[row] ) 
                
                # add in the 2 tau values, sigma, and eta value
                pf_current.extend( pf_from_db[ 0:4 ] )
                pf_delta_current.extend( pf_delta_from_db[ 0:4 ] )
                
                # add in the A and mu values
                A = pf_from_db[ mu_index_in_pf - 1 ]
                mu = pf_from_db[ mu_index_in_pf ]
                A_delta = pf_delta_from_db[ mu_index_in_pf - 1 ]
                mu_delta = pf_delta_from_db[ mu_index_in_pf ]
                
                # extend current arrays 
                pf_current.extend( [A, mu] )
                pf_delta_current.extend( [A_delta, mu_delta ] )
                
                # add to pf_arr
                pf_arr.append( pf_current )
                pf_delta_arr.append( pf_delta_current )
                
                
                # otherwise we can't extract single-fit params for both
                # peaks.
                
            else:
                pf_arr.append( [ np.nan ] * 6  )
                pf_delta_arr.append( [ np.nan ] * 6 ) 
                
        # print( pf_arr )
        # print( pf_delta_arr ) 
                
        # return the values and deltas as a measurement.
        return meas.meas( pf_arr, pf_delta_arr ) 
            




    # input filename of database, return DataFrame containing the DB  
    def read_db_into_df( self ):
        self.assert_open()
        return pd.read_sql_query( 'SELECT * from ' + tablename, self.conn )
    
    

        
    # if populating an entire array, we can get a bit more efficiency by not calling
    # get_mu_values
    def get_mu_grid_where_valid( self, peaknum ):
        
        df = self.read_db_into_df()
        
        fitnum, peaknum_in_pf = get_fitnum_and_mu_index_in_pf( peaknum )
        
        successful_fit_values = df.successful_fit.values[ fitnum : : NUM_FITS_PER_PIXEL ]
        
        pfs = df.pf.values [ fitnum : : NUM_FITS_PER_PIXEL ]
        
        pf_deltas = df.pferr.values[ fitnum : : NUM_FITS_PER_PIXEL ]
        
        # return meas constructed from the mu and mu_delta values
        # return np.asarray( [ json.loads(pfs[i])[ peaknum_in_pf ] if successful_fit_values[i]
        #                 else np.nan for i in range(pfs.size) ] ).reshape(32,32)
        
        vals = np.asarray( [ json.loads(pfs[i])[ peaknum_in_pf ] if successful_fit_values[i]
                             else np.nan for i in range(pfs.size) ] ).reshape(32,32)
        
        deltas = np.asarray( [ json.loads( pf_deltas[i][ peaknum_in_pf ] ) if successful_fit_values[i]
                               else np.nan for i in range(pf_deltas.size) ] ).reshape(32,32) 
        
        return meas.meas( vals, deltas )
    



###########################################
## DECLARE GLOBALS FOR USE ELSEWHERE ######
###########################################

all_dbs = [ db( name ) for name in [ 'centered', 'moved', 'flat', 'angled' ] ]
centered, moved, flat, angeled = all_dbs



# fitnum is the fit that this peak belongs to and index_in_pf is 
# the index of the peak within pf. a dedicated function is necessary because of 
# the convoluted way in which the data is stored.

def get_fitnum_and_mu_index_in_pf( peaknum ):
    return ( peaknum // 2, 5 + 2* (peaknum % 2) )


    
