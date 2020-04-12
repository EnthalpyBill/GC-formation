'''
Mass evolution of GC

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np

# Note: all times in [Gyr]

# ***** Dynamic evolution of GC in Choksi & Gnedin (2018) *****

def t_tid_cg18(m):
    # Tidally-limited disruption timescale in Choksi & Gnedin (2018)
    p = 0.5
    return 5 * ((m/2e5)**(2/3)) * (p/0.5)

def GC_dyn_evolv_cg18(m0, t):
    # Dynamic evolution of GC in Choksi & Gnedin (2018)
    return (1 - (2/3)*(t/t_tid_cg18(m0)))**(2/3) # m(t)/m0


# ***** Stellar evolution of GC in Prieto & Gnedin (2008) *****

def GC_star_evolv_pg08(m0, t):
    # Simple version!
    return 0.6 # m(t)/m0


# ***** Mass evolution of GC in Choksi & Gnedin (2018) *****

def GC_mass_evolv_cg18(m0, t):
    # Mass evolution of GC in Choksi & Gnedin (2018)
    return m0 * GC_dyn_evolv_cg18(m0, t) * GC_star_evolv_pg08(m0, t)
