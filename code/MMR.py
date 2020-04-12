'''
Mass-metallicity relation

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np

# ***** MMR in Li & Gnedin (2014) *****

def mmr_lg14(m_star, z):
    n = 2.8
    mmr = 0.4*np.log10(m_star / (10**10.5)) - (0.216*n*np.log10(1+z))
    return mmr + (0.2-mmr)*np.array(mmr>0.2, dtype=np.int)


# ***** MMR in Choksi et al. (2018) *****

def mmr_cgl18(m_star, z):
    mmr = 0.35*np.log10(m_star / (10**10.5)) - (0.9*np.log10(1+z))
    return mmr + (0.3-mmr)*np.array(mmr>0.3, dtype=np.int)
