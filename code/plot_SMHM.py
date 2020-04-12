'''
Plot stellar mass-halo mass relation

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np
from SMHM import log_smhm_b13
import BSG_ratio as bsg
import matplotlib.pyplot as plt

log_mh = np.linspace(9.5, 14, 100+1)
z = np.array([0, 1, 3, 5])
c = ["k", "r", "b", "g"]

for i in range(len(z)):
    plt.plot(log_mh, np.log10(bsg.star_ratio_lg14(10**log_mh, z[i])), c=c[i], ls="--")
    plt.plot(log_mh, np.log10(bsg.gas_ratio(10**log_mh, z[i])), c=c[i], ls="-")

plt.show()
