# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 13:17:54 2020

@author: Olly
"""

import numpy as np
#constants
c = 3*(10**8)

def E_field(z, t, k, w, E0, phase, angle):
    phi = (k*z - w*t - phase)
    E_x = (np.cos(angle))*E0*(np.sin(phi))
    E_y = (np.sin(angle))*E0*(np.sin(phi))
    E = np.array([E_x, E_y, 0])
    return E

def B_field(z, t, k, w, E0, phase, angle):
    phi = (k*z - w*t - phase)
    B_x = (np.cos(angle+(np.pi/2)))*(E0*k/w)*(np.sin(phi))
    B_y = (np.sin(angle+(np.pi/2)))*(E0*k/w)*(np.sin(phi))
    B = np.array([B_x, B_y, 0])
    return B

def a_field(z, t, k, w, E0, phase, angle):
    phi = (k*z - w*t - phase)
    a_x = (np.cos(angle))*(-E0/(w))*(np.cos(phi))
    a_y = (np.sin(angle))*(-E0/(w))*(np.cos(phi))
    a =np.array([a_x,a_y,0])
    return a
