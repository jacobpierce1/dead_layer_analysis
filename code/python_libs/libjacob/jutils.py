import time
from math import log10, floor
import numpy as np 


# check if x is a python int or a numpy int.
def isint( x ):
    return isinstance( x, (int, np.integer ) )


class time_estimator( object ) :

    def __init__( self, num_iterations, num_updates ) :
        self.iteration = 0 
        self.counter = 1
        self.num_updates = num_updates 
        self.num_iterations = num_iterations
        self.start_time = time.time()

        
    def update( self ) :
        self.iteration += 1

        if( 1. * self.iteration / self.num_iterations 
            >= 1. * self.counter / self.num_updates ): 

            current_time = time.time()
            self.counter += 1
            print( "%d/%d complete, %f mins remaining" \
                   % ( self.counter, self.num_updates,
                       (current_time - self.start_time) / 60.0
                       * ( self.num_updates - self.counter )
                       / self.counter ) )

            
    def reset( self, num_updates = None ) :
        self.iteration = 0 
        self.counter = 1
        self.start_time = time.time()
        if num_updates is not None :
            self.num_updates = num_updates 
        
      




# Flattne a list of lists.
def flatten_list( list_ ):
    return [ x for sublist in list_ for x in sublist ]




# input: x and uncertainty dx. round dx to nearest sig fig, figure out what it is,
# and round x to the same thing. they are returned as strings.
def sigfig( x, dx, sci=0 ):
    
    precision = -int(floor(log10(abs(dx))))
    
    ## edge case: add an extra digit of precision if we are starting with the digit 1.
    #if( ('%+e' % dx)[1] == '1' ):
    #    precision -= 1
    
    # perform the rounding
    dx_rounded =  round(dx,precision)
    x_rounded = round(x,precision)
        
    ret = [""] * 2
    
    if( dx_rounded >=1 ):
        ret = [ str(int(var)) for var in [x_rounded, dx_rounded] ]
    else:
        ret = [ str(float(var)) for var in [x_rounded, dx_rounded] ]
        
        #handle edge case: x_rounded has lower precision than dx_rounded, the convention
        # is then to append 0s until the precision matches.
        dx_rounded_precision = len(ret[1]) - ret[1].find('.')
        x_rounded_precision = len(ret[0]) - ret[0].find('.')
        for i in range( dx_rounded_precision - x_rounded_precision ):
            ret[0] += "0"
    
    
    return ret 
    
    
    

# take an array containing floats, another array of the same size containing floats,
# and create a string of the form "$ (A1, .., An) = (p[0] \pm perr[0], ... ) $"
def format_measurement_vector( variable_name, p, perr ):
    
    # if we are not using a list, thn this is the precision for the measurement.
    use_uncert = isinstance( perr, list )  # whether or not to include uncertainties in the string.
    
    if use_uncert:
        if( len(p) != len(perr) ):
            print( "ERROR: size of p and perr are different" )
            return ""
    
    use_paren = len(p) > 1 
    
    variable_str = ""
    if(use_paren):    variable_str += "("
    variable_str += ", ".join( [variable_name + "_" + str(i) for i in range(len(p)) ] )
    if(use_paren):    variable_str += ")"
    
    values_str = ""
    if(use_paren):    values_str += "("
    if use_uncert:
        values_str += ", \;".join( [ "%s \\pm %s" % tuple(sigfig(p[i], perr[i])) for i in range(len(p)) ] ) 
    else:
        values_str += ", \;".join( [ "%s" % tuple(sigfig(p[i], perr))[0] for i in range(len(p)) ] ) 
    if(use_paren):    values_str += ")"

    output_str = "$ " + variable_str + " = " + values_str + " $"
    return output_str
    


                            
    
def round_to_1(x):
    x = round(x, -int(floor(log10(abs(x)))))
    if x >= 1:
        x = int(x)
    return x




