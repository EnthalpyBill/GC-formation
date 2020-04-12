'''
Stellar mass-halo mass relation

Created Apr. 2020
Last Edit Apr. 2020

By Bill Chen
'''

import numpy as np

# ***** SMHM relation in Behroozi et al. (2013) *****

def scale_factor(z):
    return 1 / (1+z)

def cutoff(a):
    # Exponential cutoff of evolution of SMHM relation with scale factor
    return np.exp(-4*a*a)

def fun1(x, alpha, delta, gamma):
    return -np.log10(10**(alpha*x) + 1) + \
                     (delta * np.power(np.log10(1+np.exp(x)), gamma) / \
                     (1 + np.exp(10**(-x))))

def log_smhm_b13(log_mh, z):
    # SMHM relation in Behroozi et al. (2013)

    a = scale_factor(z)
    # Exponential cutoff of evolution of SMHM relation with scale factor
    nu = cutoff(a)
    # Characteristic stellar mass to halo mass ratio
    log_epsilon = -1.777 + (-0.006*(a-1) + \
                  0.000*z)*nu - 0.119*(a-1)
    # Characteristic halo mass
    log_m1      = 11.514 + (-1.793*(a-1) + \
                  (-0.251)*z)*nu
    # Faint-end slope of SMHM relation
    alpha       = -1.412 + 0.731*(a-1)*nu
    # Faint-end slope of SMHM relation
    delta       = 3.508 + (2.608*(a-1) + \
                  (-0.043)*z)*nu
    # Index of subpower law at massive end of SMHM relation
    gamma       = 0.316 + (1.319*(a-1) + \
                  0.279*z)*nu

    return log_epsilon + log_m1 + \
           fun1(log_mh-log_m1, alpha, delta, gamma) - \
           fun1(0, alpha, delta, gamma)

# ***** Eq. (3) in Behroozi et al. (2013) *****
