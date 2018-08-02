import matplotlib.pyplot as plt
import numpy as np


x = np.linspace( 0, 100, 100 ) 

f = plt.figure( figsize = (5,5) )

ax = plt.axes()

ax.plot( x, x ) 

plt.show()
