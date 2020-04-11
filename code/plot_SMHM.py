'''
Plot stellar mass-halo mass relation

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np
from SMHM import log_smhm_b13
import matplotlib.pyplot as plt

log_mh = np.linspace(10, 15, 100+1)
z = [0,1,2,3,5]

for i in range(len(z)):
    plt.plot(log_mh, log_smhm_b13(log_mh, z[i])-log_mh)

plt.show()
