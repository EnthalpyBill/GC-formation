'''
Mass evolution of halo

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np

# Note: all times in [Gyr]

# ***** Mass evolution of halo in Wechsler et al. (2001) *****

def alpha_w01(m0):
    # Simple version!
    return 1

def halo_mass_evolv_w01(m0, z):
    # Mass evolution of halo in Wechsler et al. (2001)
    return m0*np.exp(-alpha_w01(m0)*z)
