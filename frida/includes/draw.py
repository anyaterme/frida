import matplotlib.pyplot as plt
import numpy as np

f = np.loadtxt('skyradiance_mean.dat')
x = f[:,0]
y = f[:,1]

plt.plot(x,y)
plt.show()
