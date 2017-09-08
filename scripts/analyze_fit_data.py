## includes 
import numpy as np
import matplotlib.pyplot as plt
import time 
import sys
# import simplejson as json 
import json


# this script is meant to be run after parse_all_data.py. that script writes all the data
# about all the fits (including failed ones) to a json file. that file is then read 
# here and analyzed. there is no data processing performed in the previous script 
# other than fitting the functions. in this script we reconstruct the fit functinos 
# using the fit parameters in the json file and pick up the analysis from there.
# this is because the curve fitting process takes a while but is final, the data
# analysis phase can take much longer. 

