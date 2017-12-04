import matplotlib.pyplot as plt

import numpy as np 

x = np.array( [1,3,2,5] )
y = x ** 2
dy = x / 2

ax = plt.axes()
ax.errorbar( x, y, dy, ls='none' )

plt.show() 
