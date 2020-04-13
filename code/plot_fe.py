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

x, y = np.loadtxt("./output/fe_distrib_25_07.txt")

ax1.plot(x, y, c='k')

ax1.set_xlim([-2.2,0])
ax1.set_ylim([0, 1])
plt.xticks([-2, -1.5, -1, -0.5, 0],
           [r'$\rm -2$', r'$\rm -1.5$', r'$\rm -1$',
            r'$\rm -0.5$', r'$\rm 0$'])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1],
           [r'$\rm 0$', r'$\rm 0.2$', r'$\rm 0.4$',
            r'$\rm 0.6$', r'$\rm 0.8$', r'$\rm 1$'])
ax1.set_xlabel(r'$\rm [Fe/H]$')
ax1.set_ylabel(r'$\rm Number\ Density$')

plt.minorticks_off()

plt.savefig('./figure/fe')
plt.clf()
