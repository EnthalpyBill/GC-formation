'''
Plot stellar mass-halo mass relation

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np
from SMHM import log_smhm_b13
import BSG_ratio as bsg
import MMR as mmr
import matplotlib.pyplot as plt

log_mh = np.linspace(9.5, 14, 100+1)
mh = 10**log_mh
z = np.array([0, 1, 3, 5])
c = ["k", "r", "b", "g"]

for i in range(len(z)):
    plt.plot(mh, bsg.star_ratio_lg14(mh, z[i]), c=c[i], ls="--")
    # plt.plot(mh, bsg.gas_ratio_lg14(mh, z[i]), c=c[i], ls=":")
    plt.plot(mh, bsg.gas_ratio_cgl18(mh, z[i]), c=c[i], ls="-")

plt.xscale('log')
plt.yscale('log')
plt.show()
