# this script is made for managing the DB, including creating the initial schema,
# adding columns later for data analysis, wrappers for reading / writing to the 
# schema, and 
# 
# source: https://pymotw.com/3/sqlite3/

import os
import sqlite3
import json


DEBUG_DB = 0


db_filename = '../fits_and_spectral_data.db'
schema_filename = 'initial_schema.sql'
tablename = 'fits_and_extracted_data'

schema_cols = [ 'x', 'y', 'fit_id', 'successful_fit', 'fit_attempt', 'reduc_chisq', 
                    'pf', 'pferr', 'p0', 'fit_bounds', 'fwhm_data'  ]





             
                                       
def get_schema_cols():
    pass





def refresh_needs_processing():
    pass



# todo: add a col to the db which was not there at time of creation
def add_col( colname, datatype ):
    pass





# this function shall only write the fit parametetrs and things that can be obtained only from the histograms,
# especially the fits. the fwhms etc. can be extracted later much more rapidly by reconstructing the fit function.
# p0 and fit_bounds also saved because they could be used to reconstruct the fit.
# the default behavior is to overwrite data, this is to keep the function as efficient as 
# possible. there several instances i can think of in which you would want to put data
# in without bothering to check the current state of the DB. as such you have to do such a 
# query beforehand if necessary. 

def insert_fit_data_into_db( sql_conn, pixel_coords, fit_id, successful_fit, \
        fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data=[0,[0,0],[0,0]], \
        db_is_empty=0 ):        
    
    # the c programmer within.
    query = ''
    
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
                
            
    if DEBUG_DB:
        print query 
    
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
        schema_cols[10]:json.dumps(fwhm_data) #[0,[0,0],[0,0]]
    }

    cursor.execute( query, query_dict )
    sql_conn.commit()





def read_data_from_db( sql_conn, pixel_coords, fit_id ):
        
    cursor = sql_conn.cursor()
    
    query = 'select ' + ', '.join( schema_cols[3:] )   \
        + ' from ' + tablename    \
        + ' where ' + ' and '.join( [ col + '=:' + col for col in schema_cols[:3] ] )   
             
        
    if DEBUG_DB:
        print query
    
    query_dict = \
    {
        schema_cols[0]:pixel_coords[0],
        schema_cols[1]:pixel_coords[1],
        schema_cols[2]:fit_id
    }
    
    cursor.execute( query, query_dict )
    
    result = cursor.fetchone()
    
    if DEBUG_DB:
        print result    
            
    successful_fit, fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data = result
        
    pf = json.loads(pf)
    pferr = json.loads(pferr)
    p0 = json.loads(p0)
    fit_bounds = json.loads(fit_bounds)
    fwhm_data = json.loads(fwhm_data)
    
    return ( successful_fit, fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data )
    




# only needs to be called once. has guard against future calls (unless the 
# location of the db is changed, then it won't work)
def create_db():

    db_is_new = not os.path.exists(db_filename)
    
    with sqlite3.connect(db_filename) as conn:
        if db_is_new:
            print('INFO: creating DB and schema...')
            with open(schema_filename, 'rt') as f:
                schema = f.read()
            conn.executescript(schema)    

        else:
            print('ERROR: database already exists, returning')



# fill the db with each fit id that we need, giving the needs update flag for
# everything. this is meant to only be called when the table is empty 
def populate_db(conn, numx, numy, numfits):
    
    for i in range(numx):
        for j in range(numy):
            for k in range(numfits):
                insert_fit_data_into_db( conn, (i,j), k, 0, -1, -1, [-1], [-1], [-1], \
                        [-1], fwhm_data=[-1,[-1,-1],[-1,-1]], db_is_empty=1 ) 





def delete_db():
    ans = raw_input( 'PROMPT: delete DB, are you sure (y/n) ?  ' )
    if( ans == 'y' ):
        os.remove( db_filename ) 
        print 'INFO: deleted db.' 
        return 1
    print 'INFO: did not delete db.'
    return 0
    


# todo: when making a new db with some matching col names, fill the matching col
# entries. tricky to implement...
def repopulate_db():
    pass



def reset_fit_attempt_number():
    pass


# this should be called after a col is added to the schema or a change is made to
# the way the db is organized. delete the old db, create a new one, and then populate it.
# todo: option to populate relevant entries of the new db with ones from the old db.
def recreate_db():
    
    if not delete_db():
        return 0
    
    create_db()
    
    with sqlite3.connect( db_filename ) as sql_conn:
        populate_db(sql_conn, 32, 32, 3, )
    
    print 'INFO: successfully populated the db with defaults.'
    return 1
    
    
    
#with sqlite3.connect( db_filename ) as sql_conn:
#
#    # insert_fit_data_into_db( sql_conn, [0,0], 3, 1, [0,0], [0,0], [0,0], [0,0], [0,[0,0],[0,0]], 1 )
#    #successful_fit, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data  \
#    #                = read_data_from_db( sql_conn, [16,16], 3 )
#    #
#    result = read_data_from_db( sql_conn, [16,16], 2 )
    
## add new col: 0
# https://stackoverflow.com/questions/4253804/insert-new-column-into-table-in-sqlite
