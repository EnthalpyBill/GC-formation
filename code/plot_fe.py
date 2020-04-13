'''
Plot stellar metallicity distribution

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('niceplot2')

fig, ax1 = plt.subplots()

fe, n_fe = np.loadtxt("./output/GC_formation_25_07_3e-5.txt")

bins_fe = np.linspace(-3,0,100+1)

f_fe = np.histogram(fe, weights=n_fe, bins=bins_fe)[0]

ax1.step(bins_fe[1:], f_fe, c='k')

ax1.set_xlim([-2,0])
ax1.set_ylim([0, 0.004])
plt.xticks([-2, -1.5, -1, -0.5, 0],
           [r'$\rm -2$', r'$\rm 1.5$', r'$\rm -1$',
            r'$\rm -0.5$', r'$\rm 0$'])
plt.yticks([], [])
ax1.set_xlabel(r'$\rm [Fe/H]$')
ax1.set_ylabel(r'$\rm Number Density$')

plt.minorticks_off()

plt.savefig('./figure/fe')
plt.clf()
