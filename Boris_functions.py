# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:04:56 2020

@author: Olly
"""
import numpy as np

#constants
mass_e = 9.11*(10**-31)
charge_e = -1.6*(10**-19)
c = 3*(10**8)
epp0 = 8.85*(10**-12)

def modulus(lst):
    mod = np.sqrt((lst[0]**2)+(lst[1]**2)+(lst[2]**2))
    return mod

def E_boris(p, q, E):
    p_new = p + q*E
    return p_new

def B_boris(p, q, B):
    mod_p = (p[0]**2)+(p[1]**2)+(p[2]**2)
    gamma = (1 + (mod_p)/((mass_e*c)**2))**0.5
    t = (q/(mass_e*gamma))*B
    mod_t = (t[0]**2)+(t[1]**2)+(t[2]**2)
    s = (2/(1+mod_t))*t
    p_prime = p + np.cross(p, t)
    p_plus = p + np.cross(p_prime, s)
    return p_plus

def damp_boris(p, q, E, B):
    mod_p = (p[0]**2)+(p[1]**2)+(p[2]**2)
    gamma = (1 + (mod_p)/((mass_e*c)**2))**0.5
    #v = c*p/((((mass_e*c)**2)+mod_p)**0.5)
    v = (1/(gamma*mass_e))*p
    factor_0 = np.dot(v, v)
    if factor_0 == 0:
        p_fin = p
    else:
        factor_1 = (-1)*((charge_e)**3)/(3*(np.pi)*(epp0)*(mass_e**3)*(c**5))
        dot_factor = ((np.dot(E, v))/(np.dot(v, v)))
        E_perp = E - dot_factor*v
        B_term = np.cross(v, B)
        factor_2 = E_perp + B_term
        factor_3 = np.dot(factor_2, factor_2)
        p_fin = p + q*(factor_1)*(gamma)*(factor_3)*p
    return p_fin


