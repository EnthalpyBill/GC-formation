'''
GC formation

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
fb=0.16

m_halo_0 = 1e12
m_halo_c = 1e9

z_start = np.log(m_halo_0/m_halo_c) # Simple version of w01!
z_end = 0

z_slice = 10000
z_list = np.linspace(z_start, z_end, z_slice+1)

GC_m0_list = np.zeros(z_slice)
GC_mass_list = np.zeros(z_slice)
GC_time_list = np.zeros(z_slice)
GC_fe_list = np.zeros(z_slice)

pop_slice = 100
GC_pop_mat = np.zeros([z_slice, pop_slice]) # population matrix
log_GC_pop_m = np.linspace(6, 10, pop_slice+1)
log_GC_pop_m_mean = (log_GC_pop_m[1:] + log_GC_pop_m[:-1]) / 2
log_GC_pop_f = log_GC_pop_m * (-2) + 20
log_GC_pop_f = (log_GC_pop_f[1:] + log_GC_pop_f[:-1]) / 2
GC_pop_m = 10**log_GC_pop_m
GC_pop_dm = GC_pop_m[1:] - GC_pop_m[:-1]
GC_pop_m_mean = 10**log_GC_pop_m_mean
GC_pop_f = 10**log_GC_pop_f
GC_pop_p = GC_pop_f*GC_pop_dm / np.sum(GC_pop_f*GC_pop_dm)
GC_pop_m_tot = np.sum(GC_pop_p * GC_pop_m_mean)

for i in range(z_slice):
    z0 = z_list[i]
    z1 = z_list[i+1]
    t0 = -cosmo.lookback_time(z0).value # current time is zero
    t1 = -cosmo.lookback_time(z1).value
    GC_time_list[i] = t1

    mh0 = hme.halo_mass_evolv_w01(m_halo_0, z0)
    mh1 = hme.halo_mass_evolv_w01(m_halo_0, z1)
    halo_dm = mh1 - mh0
    halo_dm_rate = halo_dm/(t1-t0)/mh1
    # print(halo_dm/(t1-t0)/mh1)

    m_star = bsg.star_ratio_lg14(mh1, z1) * fb * mh1
    m_gas = bsg.gas_ratio_cgl18(mh1, z1) * fb * mh1

    GC_dm = 0.5*(np.tanh(2.5*(halo_dm_rate-0.7))+1) * (t1-t0) * \
            (3e-5) * m_gas / fb
    GC_m0_list[i] = GC_dm

    GC_pop = GC_pop_p * GC_dm / GC_pop_m_tot
    GC_pop_mat[i] = GC_pop

    metallicity = mmr.mmr_lg14(m_star, z1)
    GC_fe_list[i] = metallicity

GC_pop_mat_mass = np.zeros([z_slice, pop_slice])
# update mass
# for j in range(z_slice):
#     for k in range(pop_slice):
#         GC_pop_mat_mass[j][k] = gme.GC_mass_evolv_cg18(GC_pop_m_mean[k], t1 - GC_time_list[j])

print(GC_pop_m_tot)

# np.savetxt("./output/GC_formation_25_07_3e-5.txt",
#            [GC_fe_list, np.sum(GC_pop_mat, axis=1)])
