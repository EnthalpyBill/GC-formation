'''
Baryon, star, and gas ratios to dark matter

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np
from SMHM import log_smhm_b13
from astropy.cosmology import FlatLambdaCDM

# ***** Baryon ratio in Muratov & Gnedin (2010) *****

def cutoff_mass_mg10(z, h, H02Hz):
    m1 = 3.6 * (10**9) * np.exp(-0.6*(1+z)) / h
    m2 = 1.5 * (10**10) * (180**(-0.5)) * (H02Hz) / h
    return m2 + (m1-m2)*np.array(m1>m2, dtype=np.int) # max

def baryon_ratio_mg10(mh, z, H0=70, Om0=0.3):
    # Baryon ratio in Muratov & Gnedin (2010)

    h = H0 / 100
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    Hz = cosmo.H(z).value
    H02Hz = H0 / Hz

    mc = cutoff_mass_mg10(z, h, H02Hz)
    return (1 + mc/mh)**(-3)

# ***** Baryon ratio in Muratov & Gnedin (2010) *****

# ***** Star ratio in Li & Gnedin (2014) *****

def star_ratio_lg14(mh, z, fb=0.16):
    # Star ratio in Li & Gnedin (2014)

    log_mh = np.log10(mh)
    log_m_star = log_smhm_b13(log_mh, z)
    return 10**(log_m_star - log_mh) / fb

# ***** Star ratio in Li & Gnedin (2014) *****

# ***** Gas ratio in Li & Gnedin (2014) *****

def alpha_lg14(m_star):
    # power-law index of eta in Li & Gnedin (2014)
    return 0.19 + (0.68-0.19)*np.array(m_star>1e9, dtype=np.int)

def eta_lg14(z, m_star, alpha, n=2.8):
    # Gas mass to star mass ratio
    return 1.8 * (m_star/1e9)**(-alpha(m_star)) * ((1+z)**n)

def gas_ratio(mh, z, fb=0.16, alpha=alpha_lg14, eta=eta_lg14, H0=70, Om0=0.3):
    # Universal gas ratio
    star_ratio = star_ratio_lg14(mh, z, fb)
    baryon_ratio = baryon_ratio_mg10(mh, z, H0, Om0)
    m_star = star_ratio * mh * fb

    fg1 = star_ratio * eta(z, m_star, alpha_lg14)
    fg2 = baryon_ratio - star_ratio
    return fg1 + (fg2-fg1)*np.array(fg1>fg2, dtype=np.int) # min

def gas_ratio_lg14(mh, z, fb=0.16):
    # Gas ratio in Li & Gnedin (2014)
    return

# ***** Gas ratio in Li & Gnedin (2014) *****
