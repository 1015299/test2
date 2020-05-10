# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 14:31:09 2020

@author: Olly
"""

#imports 
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from Boris_functions import E_boris
from Boris_functions import B_boris
from Boris_functions import damp_boris
from field_functions4 import E_field 
from field_functions4 import B_field
from field_functions4 import a_field

#constants
mass_e = 9.11e-31
charge_e = -1.6e-19
c = 2.998e8
epp0 = 8.85e-12
Thom_cross = 6.652e-29
red_planck = 1.0546e-34

#functions
def modulus(lst):
    mod2 = 0
    for i in range(len(lst)):
        mod2 += lst[i]**2
    mod = np.sqrt(mod2)
    return mod

def Gamma_p(p):
    modp = modulus(p)
    gamma = np.sqrt(1 + (modp**2)/((mass_e*c)**2))
    return gamma

def eta(p, E, B):
    factor_1 = (charge_e*red_planck)/((mass_e**2)*(c**3))
    gam = Gamma_p(p)
    v = p/(gam*mass_e)
    factor_0 = np.dot(v, v)
    if factor_0 == 0:
        factor_2 = E
        Eta = abs(factor_1*modulus(factor_2) )
    else:
        dot_factor = ((np.dot(E, v))/(np.dot(v, v)))
        E_perp = E - dot_factor*v
        B_term = np.cross(v, B)
        factor_2 = E_perp + B_term
        factor_3 = dot_factor*(modulus(v))/gam
        factor_4 = np.dot(factor_2, factor_2) + np.dot(factor_3, factor_3)
        Eta = gam*abs(factor_1*np.sqrt(factor_4))
    return Eta

#beam characheristics
pulse_length = 3.4e-14
lamda = 1e-6
k = 2*np.pi/lamda
w = c*k
rc = 7

e = (1/6)*(charge_e**2)*w/((np.pi)*(epp0)*mass_e*(c**3))
I = ((mass_e**2)*(w**2)*(c**3)*(epp0)/(4*(charge_e**2)))*((rc/e)**(2/3))
print(I)
ratio = 0.2
#m = 10
#pulse_length = (2*m)*(np.pi/w)
E0 = np.sqrt((2*I)/(c*epp0))
alpha = ((np.pi)/(w*pulse_length))
k0 = alpha*k

#initial conditions and setup
part = 1500
dt = 2*(np.pi)/(part*w)
n = 1
iter_num = int(n*pulse_length/dt + part/2)
num_pos = 50

#position, momentum and E_field
x_ = []
y_ = []
z_ = []
px_ = []
py_ = []
pz_ = []
t_ = []
gam_ = []
eta_ = []
x_1 = []
y_1 = []
z_1 = []
px_1 = []
py_1 = []
pz_1 = []
t_1 = []
gam_1 = []
eta_1 = []

#field functions
def e_field(z, t):
    e1 = E_field(z, t, k, w, k0, E0, 0, 0)
    e2 = E_field(z, t, k, w, k0, E0, -(np.pi/2), (np.pi/2))
    e3 = ratio*E_field(z, t, -k, w, -k0, E0, 0, 0)
    e4 = ratio*E_field(z, t, -k, w, -k0, E0, -(np.pi/2), (np.pi/2))
    E1 = (np.sqrt(0.5))*(e1 + e2 + e3 + e4)
    return E1

def b_field(z, t):
    b1 = B_field(z, t, k, w, k0, E0, 0, 0)
    b2 = B_field(z, t, k, w, k0, E0, -(np.pi/2), (np.pi/2))
    b3 = ratio*B_field(z, t, -k, w, -k0, E0, 0, 0)
    b4 = ratio*B_field(z, t, -k, w, -k0, E0, -(np.pi/2), (np.pi/2))
    E1 = (np.sqrt(0.5))*(b1 + b2 + b3 + b4)
    return E1

def A_field(z, t):
    A1 = a_field(z, t, k, w, k0, E0, 0, 0)
    A2 = a_field(z, t, k, w, k0, E0, -(np.pi/2), (np.pi/2))
    A3 = ratio*a_field(z, t, -k, w, -k0, E0, 0, 0)
    A4 = ratio*a_field(z, t, -k, w, -k0, E0, -(np.pi/2), (np.pi/2))
    a1 = (np.sqrt(0.5))*(A1 + A2 + A3 + A4)
    return a1

for j in range(num_pos):
    print(j)
    x=[0]
    y=[0]
    z=[(j-(num_pos/2))*lamda/num_pos]
    px= [0]
    py= [0]
    pz= [0]
    l = 0
    l1 = 0
    Ein = e_field(z[0], (-1)*part*dt/(2))
    Bin = b_field(z[0], (-1)*part*dt/(2))
    p_in = np.array([px[0], py[0], pz[0]])
    ETa = [eta(p_in, Ein, Bin)]
    for i in range(iter_num):
        r_old = np.array([x[i], y[i], z[i]])
        p_old = np.array([px[i], py[i], pz[i]])
        T = (i-part/2)*dt
        E = e_field(z[i], T)
        B = b_field(z[i], T)
        A = A_field(z[i], T)
        q = (charge_e*dt)/2
        p1 = E_boris(p_old, q, E)
        p2 = B_boris(p1, q/2, B)
        p3 = damp_boris(p2, q, E, B)
        p4 = B_boris(p3, q/2, B)
        p5 = E_boris(p4, q, E)
        px.append(p5[0])
        py.append(p5[1])
        pz.append(p5[2])
        #checks
        gamma = Gamma_p(p5)
        #positions
        modp2 = modulus(p5)**2
        v_new = c*p5/(np.sqrt(((mass_e*c)**2)+modp2))
        r_new = r_old + dt*v_new
        x.append(r_new[0])
        y.append(r_new[1])
        z.append(r_new[2])
        E2 = e_field(z[i+1], T)
        B2 = b_field(z[i+1], T)
        ETa.append(eta(p5, E2, B2))
        if i%(part) == 0:
            l += 1
            x_.append((1e6)*x[i])
            y_.append((1e6)*y[i])
            z_.append((1e6)*z[i])
            px_.append(px[i]/(mass_e*c))
            py_.append(py[i]/(mass_e*c))
            pz_.append(pz[i]/(mass_e*c))
            t_.append(T)
            gam_.append(gamma)
            eta_.append(ETa[i])
        if i%(part/15) == 0:
            l1 += 1
            x_1.append((1e6)*x[i])
            y_1.append((1e6)*y[i])
            z_1.append((1e6)*z[i])
            px_1.append(px[i]/(mass_e*c))
            py_1.append(py[i]/(mass_e*c))
            pz_1.append(pz[i]/(mass_e*c))
            t_1.append(T)
            gam_1.append(gamma)
            eta_1.append(ETa[i])
            

#colour map
colors = cm.rainbow(np.linspace(0, 1, l))
colors2 = cm.rainbow(np.linspace(0, 1, l))
for k in range((num_pos)-1):
    colors = np.append(colors, colors2, axis = 0)

colors3 = cm.rainbow(np.linspace(0, 1, l1))
colors4 = cm.rainbow(np.linspace(0, 1, l1))
for k in range((num_pos)-1):
    colors3 = np.append(colors3, colors4, axis = 0)
        
#graphs
plt.figure()
plt.xlim(1.1*min(t_), 1.1*max(t_))
plt.ylim(1.1*min(z_), 1.1*max(z_))
plt.scatter(t_, z_, s = 0.6 ,c = colors)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$z/\mu m$', fontsize = 16)
plt.savefig('tz_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(px_), 1.1*max(px_))
plt.ylim(1.1*min(py_), 1.1*max(py_))
plt.scatter(px_, py_, s = 0.6, c = colors)
plt.xlabel('$p_x$/$m_e c$', fontsize = 16)
plt.ylabel('$p_y$/$m_e c$', fontsize = 16)
plt.savefig('pxpy_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(z_), 1.1*max(z_))
plt.ylim(1.1*min(pz_), 1.1*max(pz_))
plt.scatter(z_, pz_, s=0.6, c = colors)
plt.xlabel('$z/\mu m$', fontsize = 16)
plt.ylabel('$p_z$/$m_e c$', fontsize = 16)
plt.savefig('zpz_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(t_), 1.1*max(t_))
plt.ylim(1.1*min(pz_), 1.1*max(pz_))
plt.scatter(t_, pz_, s=0.6, c = colors)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$p_z$/$m_e c$', fontsize = 16)
plt.savefig('tpz_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(t_), 1.1*max(t_))
plt.ylim(0, 1.1*max(gam_))
plt.scatter(t_, gam_, s = 0.6 ,c = colors)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$\gamma$', fontsize = 16)
plt.savefig('tgamma_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(z_), 1.1*max(z_))
plt.ylim(0, 1.1*max(gam_))
plt.scatter(z_, gam_, s = 0.6 ,c = colors)
plt.xlabel('$z/\mu m$', fontsize = 16)
plt.ylabel('$\gamma$', fontsize = 16)
plt.savefig('zgamma_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

#graphs
plt.figure()
plt.xlim(1.1*min(t_1), 1.1*max(t_1))
plt.ylim(1.1*min(z_1), 1.1*max(z_1))
plt.scatter(t_1, z_1, s = 0.6 ,c = colors3)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$z/\mu m$', fontsize = 16)
plt.savefig('tz_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(px_1), 1.1*max(px_1))
plt.ylim(1.1*min(py_1), 1.1*max(py_1))
plt.scatter(px_1, py_1, s = 0.6, c = colors3)
plt.xlabel('$p_x$/$m_e c$', fontsize = 16)
plt.ylabel('$p_y$/$m_e c$', fontsize = 16)
plt.savefig('pxpy_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(z_1), 1.1*max(z_1))
plt.ylim(1.1*min(pz_1), 1.1*max(pz_1))
plt.scatter(z_1, pz_1, s=0.6, c = colors3)
plt.xlabel('$z/\mu m$', fontsize = 16)
plt.ylabel('$p_z$/$m_e c$', fontsize = 16)
plt.savefig('zpz_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env__norm.pdf')

plt.figure()
plt.xlim(1.1*min(t_1), 1.1*max(t_1))
plt.ylim(1.1*min(pz_1), 1.1*max(pz_1))
plt.scatter(t_1, pz_1, s=0.6, c = colors3)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$p_z$/$m_e c$', fontsize = 16)
plt.savefig('tpz_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env__norm.pdf')

plt.figure()
plt.xlim(1.1*min(t_1), 1.1*max(t_1))
plt.ylim(0, 1.1*max(gam_1))
plt.scatter(t_1, gam_1, s = 0.6 ,c = colors3)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$\gamma$', fontsize = 16)
plt.savefig('tgamma_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(z_1), 1.1*max(z_1))
plt.ylim(0, 1.1*max(gam_1))
plt.scatter(z_1, gam_1, s = 0.6 ,c = colors3)
plt.xlabel('$z/\mu m$', fontsize = 16)
plt.ylabel('$\gamma$', fontsize = 16)
plt.savefig('zgamma_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env_norm.pdf')

#other graphs

plt.figure()
plt.xlim(1.1*min(t_), 1.1*max(t_))
plt.ylim(1.1*min(eta_), 1.1*max(eta_))
plt.scatter(t_, eta_, s = 0.6 ,c = colors)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$\eta$', fontsize = 16)
plt.savefig('etat_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(z_), 1.1*max(z_))
plt.ylim(1.1*min(eta_), 1.1*max(eta_))
plt.scatter(z_, eta_, s=0.6, c = colors)
plt.xlabel('$z/\mu m$', fontsize = 16)
plt.ylabel('$\eta$', fontsize = 16)
plt.savefig('etaz_rc_'+str(rc)+'_rat_'+str(ratio)+'_final_env_norm.pdf')

#graphs
plt.figure()
plt.xlim(1.1*min(t_1), 1.1*max(t_1))
plt.ylim(1.1*min(eta_1), 1.1*max(eta_1))
plt.scatter(t_1, eta_1, s = 0.6 ,c = colors3)
plt.xlabel('$t/s$', fontsize = 16)
plt.ylabel('$\eta$', fontsize = 16)
plt.savefig('etat_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env_norm.pdf')

plt.figure()
plt.xlim(1.1*min(z_1), 1.1*max(z_1))
plt.ylim(1.1*min(eta_1), 1.1*max(eta_1))
plt.scatter(z_1, eta_1, s=0.6, c = colors3)
plt.xlabel('$z/\mu m$', fontsize = 16)
plt.ylabel('$\eta$', fontsize = 16)
plt.savefig('etaz_rc_'+str(rc)+'_rat_'+str(ratio)+'_detailed_final_env_norm.pdf')
