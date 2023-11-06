#===================================================================
# Copyright (c) 2023 Chris Balzer
#     Simple plotting function for the example results
#======================================================================
import os
import numpy as np
import matplotlib.pyplot as plt

# Load data
if os.path.exists('example/data.dat'):
    data = np.loadtxt('example/data.dat')
else:
    ValueError('No datafile found!')

# Plot the excess adsorption as a function of the bulk density
plt.figure()
plt.semilogx(data[data[:,-1] == 0,12],data[data[:,-1] == 0,4],'black',label='Stable (dilute)')
plt.semilogx(data[data[:,-1] == 1,12],data[data[:,-1] == 1,4],'red',label='Unstable')
plt.semilogx(data[data[:,-1] == 2,12],data[data[:,-1] == 2,4],'blue',label='Stable (dense)')
plt.xlabel('Polymer Bulk Volume Fraction')
plt.ylabel('Excess Polymer Adsorption')
plt.legend()
plt.savefig('example/ExcessAdsorption.png', format='png',dpi=300)

# Plot the surface tension (surface excess free energy) as a function of the bulk density
plt.figure()
plt.semilogx(data[data[:,-1] == 0,12],data[data[:,-1] == 0,11],'black',label='Stable (dilute)')
plt.semilogx(data[data[:,-1] == 1,12],data[data[:,-1] == 1,11],'red',label='Unstable')
plt.semilogx(data[data[:,-1] == 2,12],data[data[:,-1] == 2,11],'blue',label='Stable (dense)')
plt.ylim([-0.001, 0.001])
plt.xlim([1e-4, 5e-4])
plt.xlabel('Polymer Bulk Volume Fraction')
plt.ylabel('Surface Tension')
plt.legend()
plt.savefig('example/SurfaceTension.png', format='png',dpi=300)