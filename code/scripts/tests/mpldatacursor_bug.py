import numpy as np
import matplotlib.pyplot as plt
import mpldatacursor

x, y, z = np.random.random((3, 10))
fig, ax = plt.subplots()

# works as expected with this uncommented:
# ax.scatter(x, y, label = 'test' )

# doesn't work with an errorbar:
ax.errorbar(x, y, xerr=1, yerr=1, ls='none', label = 'test' )

ax.legend()

mpldatacursor.datacursor( formatter = '{label}'.format )
plt.show( )
