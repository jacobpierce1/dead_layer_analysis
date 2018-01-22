import time
from math import log10, floor
import numpy as np 


# check if x is a python int or a numpy int.
def isint( x ):
    return isinstance( x, (int, np.integer ) )




# this function estimates the remaining time for a for loop in which the 
# output is expected to take roughly the same amount of time per run.
# thiis is the equivaletn of a static variable in c 
# currently implemented to only be used once per script, could be changed by adding bool reset.
def estimate_time_left( current_iteration, total_iterations, num_updates=10, reset=0 ):

    # reset if called
    
    if reset == 1 : 
        estimate_time_left.counter = 0
        
        
    # set the start time if it's the first call to this function.
    
    if estimate_time_left.start_time == 0 :
        estimate_time_left.start_time += 1
        estimate_time_left.start_time = time.time()
        return 
        
    if( current_iteration * 1.0 / total_iterations >= estimate_time_left.counter * 1.0 / num_updates ): 
        current_time = time.time()
        estimate_time_left.counter += 1
        print( "%d/%d complete, %f mins remaining" \
               % ( estimate_time_left.counter, num_updates,
                   (current_time - estimate_time_left.start_time) / 60.0
                   * (num_updates - estimate_time_left.counter )
                   / estimate_time_left.counter ) )

estimate_time_left.counter = 0
estimate_time_left.start_time = 0





# flattne a list of lists.
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




