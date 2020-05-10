# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 11:22:19 2020

@author: Olly
"""

import numpy as np
#constants
c = 3*(10**8)

def E_field(z, t, k, w, k0, E0, phase, angle):
    phi = (k*z - w*t - phase)
    phi0 = (k0/k)*(k*z-w*t)
    if ((-1)*np.pi <= phi0 <= 0):
        E_x = (np.cos(angle))*E0*(((np.sin(phi))*(np.sin(phi0))**2)-2*(k0/k)*((np.cos(phi))*(np.sin(phi0))*(np.cos(phi0))))
        E_y = (np.sin(angle))*E0*(((np.sin(phi))*(np.sin(phi0))**2)-2*(k0/k)*((np.cos(phi))*(np.sin(phi0))*(np.cos(phi0))))
        E = np.array([E_x, E_y, 0])
    else:
        E = np.array([0, 0, 0])
    return E

def B_field(z, t, k, w, k0, E0, phase, angle):
    phi = (k*z - w*t - phase)
    phi0 = (k0/k)*(k*z-w*t)
    if ((-1)*np.pi <= phi0 <= 0):
        B_x = (np.cos(angle+(np.pi/2)))*(E0/w)*(((k)*(np.sin(phi))*(np.sin(phi0))**2)-2*(k0)*((np.cos(phi))*(np.sin(phi0))*(np.cos(phi0))))
        B_y = (np.sin(angle+(np.pi/2)))*(E0/w)*(((k)*(np.sin(phi))*(np.sin(phi0))**2)-2*(k0)*((np.cos(phi))*(np.sin(phi0))*(np.cos(phi0))))
        B = np.array([B_x, B_y, 0])
    else:
        B = np.array([0, 0, 0])
    return B

def a_field(z, t, k, w, k0, E0, phase, angle):
    phi = (k*z - w*t - phase)
    phi0 = (k0/k)*(k*z-w*t)
    if ((-1)*np.pi <= phi0 <=0):
        a_x = (np.cos(angle))*(-E0/(w))*(np.cos(phi))*((np.sin(phi0))**2)
        a_y = (np.sin(angle))*(-E0/(w))*(np.cos(phi))*((np.sin(phi0))**2)
    else:
        a_x = 0
        a_y = 0
    a =np.array([a_x,a_y,0])
    return a
