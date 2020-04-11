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

log_m_star = log_smhm_b13(log_mh, 0)

plt.plot(log_mh, log_m_star-log_mh)
plt.show()
