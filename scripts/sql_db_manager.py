# this script is made for managing the DB, including creating the initial schema,
# adding columns later for data analysis, wrappers for reading / writing to the 
# schema, and 
# 
# source: https://pymotw.com/3/sqlite3/

import os
import sqlite3
import json


db_filename = '../fits_and_spectral_data.db'
schema_filename = 'initial_schema.sql'
tablename = 'fits_and_extracted_data'

schema_cols = [ 'x', 'y', 'fit_id', 'needs_proccessing', 'reduc_chisq', 
                    'pf', 'pferr', 'p0', 'fit_bounds', 'fwhm_data'  ]

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



def get_schema_cols():
    pass




# todo: add a col to the db which was not there at time of creation
def add_col( colname, datatype ):
    pass





# this function shall only write the fit parametetrs and things that can be obtained only from the histograms,
# especially the fits. the fwhms etc. can be extracted later much more rapidly by reconstructing the fit function.
# p0 and fit_bounds also saved because they could be used to reconstruct the fit.
#schema_cols = [ 'x', 'y', 'needs_proccessing', 'reduc_chisq', 
#                    'pf', 'pferr', 'p0', 'fit_bounds', 'fwhm_data'  ]
def insert_fit_data_into_db( sql_conn, pixel_coords, fit_id, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data=[0,[0,0],[0,0]], overwrite=0 ):        
    
    query = ''
    
    if not overwrite:
        query = 'insert into ' + tablename + ' ('   \
                + ', '.join(schema_cols)    \
                + ') values ('         \
                + ':' + ', :'.join(schema_cols)   \
                + ');'
    
    else:
        query = 'update ' + tablename   \
                + ' set ' + ', '.join( [ col + '=:' + col for col in schema_cols[3:] ] )   \
                + ' where ' + ' and '.join( [ col + '=:' + col for col in schema_cols[:3] ] )   
                
    print query 
    
    cursor = sql_conn.cursor()
    
    query_dict =    \
    {
        schema_cols[0]:pixel_coords[0],
        schema_cols[1]:pixel_coords[1],
        schema_cols[2]:fit_id,
        schema_cols[3]:1, # by default say that it hasnt been processed
        schema_cols[4]:reduc_chisq,
        schema_cols[5]:json.dumps(pf),
        schema_cols[6]:json.dumps(pferr),
        schema_cols[7]:json.dumps(p0),
        schema_cols[8]:json.dumps(fit_bounds),
        schema_cols[9]:json.dumps(fwhm_data) #[0,[0,0],[0,0]]
    }

    cursor.execute( query, query_dict )
    sql_conn.commit()





def read_all_data_from_db( sql_conn, pixel_coords, fit_id ):
        
    cursor = sql_conn.cursor()
    
    query = 'select ' + ', '.join( schema_cols[3:] )   \
        + ' from ' + tablename    \
        + ' where ' + ' and '.join( [ col + '=:' + col for col in schema_cols[:3] ] )   
                
    print query
    
    query_dict = \
    {
        schema_cols[0]:pixel_coords[0],
        schema_cols[1]:pixel_coords[1],
        schema_cols[2]:fit_id
    }
    
    cursor.execute( query, query_dict )
    
    needs_processing, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data = cursor.fetchone()
        
    pf = json.loads(pf)
    pferr = json.loads(pferr)
    p0 = json.loads(p0)
    fit_bounds = json.loads(fit_bounds)
    fwhm_data = json.loads(fwhm_data)
    
    return ( needs_processing, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data )
    



def delete_db():
    ans = raw_input( 'PROMPT: delete DB, are you sure (y/n) ?' )
    if( ans == 'y' ):
        os.remove( db_filename ) 
        print 'INFO: deleted db.' 
        return 1
    print 'INFO: did not delete db.'
    return 0
    



with sqlite3.connect( db_filename, detect_types=sqlite3.PARSE_DECLTYPES, ) as sql_conn:

    insert_fit_data_into_db( sql_conn, [0,0], 3, 1, [0,0], [0,0], [0,0], [0,0], [0,[0,0],[0,0]], 1 )
    needs_processing, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data  \
                    = read_all_data_from_db( sql_conn, [0,0], 3 )
    
    
    
## add new col: 
# https://stackoverflow.com/questions/4253804/insert-new-column-into-table-in-sqlite