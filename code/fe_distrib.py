'''
Calculate stellar metallicity distribution

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np
import SMHM as smhm
import BSG_ratio as bsg
import MMR as mmr
import GC_mass_evolv as gme
import halo_mass_evolv as hme
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt

cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # cosmology
fb = 0.16

m_halo_0 = 1e12
m_halo_c = 1e9

z_start = np.log(m_halo_0/m_halo_c) # Simple version of w01!
z_end = 0

GC_pop_m_tot = 9214517.654711826 # calculated by dN/dM ~ M^-2

def dz_dt(z):
    dz = 0.0001
    t0 = cosmo.lookback_time(z-dz).value
    t1 = cosmo.lookback_time(z+dz).value
    return 2*dz / (t1-t0)

def dN_dz(z):
    # Simple version of w01!
    mh = hme.halo_mass_evolv_w01(m_halo_0, z)
    return (np.tanh(2.5*(1*dz_dt(z)-0.7))+1) * \
           bsg.gas_ratio_cgl18(mh, z) * mh / dz_dt(z)

def fe(z):
    mh = hme.halo_mass_evolv_w01(m_halo_0, z)
    return mmr.mmr_lg14(bsg.star_ratio_lg14(mh, z)*fb*mh, z)

def dfe_dz(z):
    dz = 0.0001
    fe0 = fe(z-dz)
    fe1 = fe(z+dz)
    return (fe1-fe0) / (2*dz)

def fe_distrib(z):
    return dN_dz(z) / dfe_dz(z)

z_slice = 1000
z_list = np.linspace(z_start, z_end, z_slice+1)

x = np.zeros(z_slice)
y = np.zeros(z_slice)

for i in range(z_slice):
    x[i] = fe(z_list[i])
    y[i] = fe_distrib(z_list[i])

y = y / np.trapz(y, x)
np.savetxt("./output/fe_distrib_25_07.txt",
           [x, y])
