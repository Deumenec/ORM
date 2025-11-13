# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 10:30:02 2025

@author: dhuerta
In this file, I test the formulas for all partial derivatives of optical functions
with respect to quadrupole strength provided in the report, thus I can quantify
how bad with respect to each component the linear aproximation is
"""

import at
import os
import numpy as np
import tqdm
import matplotlib.pyplot as plt

def vectorPlot(vectors, vecNames, savename=False):
    """ Creates a plot commparing the components of vectors of equal length"""
    colors = plt.cm.rainbow(np.linspace(0, 1, len(vectors)))
    for i, vector in enumerate(vectors):
        plt.plot(vector, color = colors[i], label= vecNames[i])
    plt.legend()
    plt.show()
    if (savename != False): plt.savefig("savename")
    return

def Cijn(phi_i, phi_j, nu, n):
    """Special Cosine function"""
    return(np.cos(n*abs(phi_i-phi_j)-n*np.pi*nu))
def Sijn(phi_i, phi_j, nu, n):
    """Special Sine function"""
    return(np.sign(phi_i-phi_j)*np.sin(n*abs(phi_i-phi_j)-n*np.pi*nu))


def dnu_dq(ring, ind_quad, step):
    """ Returns the derivative of the tune with respect to change along each quadrupole strength by computing the change in optics"""
    optics0    = at.get_optics(ring, refpts=ind_quad)
    quadLen    = [element.Length for element in ring[ind_quad]]
    result_ana = np.zeros([len(ind_quad)])
    result_num = np.zeros([len(ind_quad)])
    for i, ind in enumerate(ind_quad):
        ring[ind].PolynomB[1]+=step
        optics1 = at.get_optics(ring,refpts= [ind] )
        ring[ind].PolynomB[1]-=step
        #We compute numerically the derivative
        result_num[i] = (optics1[1]["tune"][1]-optics0[1]["tune"][1])/step 
        #Also we use the analytical formula
        result_ana[i] = -optics0[2]["beta"][i][1]*quadLen[i]/(4*np.pi) 
    return result_ana, result_num

def dbi_dq(ring, ind_el, ind_quad, step):
    """ Returns the derivative of the beta in elements with respect to change along each quadrupole strength by computing the change in optics"""
    optics0el= at.get_optics(ring, refpts=ind_el)
    optics0quad= at.get_optics(ring, refpts=ind_quad)
    tune       = optics0el[1]["tune"][1]
    el_beta    = [i[1] for i in optics0el[2]["beta"]]
    el_psi     = [i[1] for i in optics0el[2]["mu"]]
    quad_beta  = [i[1] for i in optics0quad[2]["beta"]]
    quad_psi   = [i[1] for i in optics0quad[2]["mu"]]
    quadLen    = [element.Length for element in ring[ind_quad]]
    
    result_ana = np.zeros([len(ind_quad), len(ind_el)])
    result_num = np.zeros([len(ind_quad), len(ind_el)])
    
    for i, ind in enumerate(ind_quad):
        #First I compute numerically the partial derivatives
        ring[ind].PolynomB[1]+=step
        optics1el = at.get_optics(ring,refpts= ind_el )
        ring[ind].PolynomB[1]-=step
        for j in range(len(ind_el)):
            #Here I calculate the derivative numerically
            result_num[i][j] = (optics1el[2]["beta"][j][1] -el_beta[j] )/step
            #Here I calculate it with the formula
            result_ana[i][j] = el_beta[i]*quad_beta[i]*quadLen[i]/(2*np.sin(2*np.sin(2*np.pi*tune)))*Cijn(el_psi[j], quad_psi[i], tune, 2)
    return result_ana, result_num

    
    
os.chdir('Z:\Projectes\ORM\ORM') #Set my working directory!
ring = at.load_mat('THERING.mat')
step = 0.000001


ind_bpm = at.get_refpts(ring, 'BPM')   # family name, adjust to your lattice
ind_cor = at.get_refpts(ring, 'COR')  # horizontal correctors
ind_Vquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QV'))
ind_Hquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QH'))


#ana_dnu, num_dnu = dnu_dq(ring, ind_Vquad, step)
#vectorPlot([ana_dnu, num_dnu], ["Analytical dtune/dq", "Numerical dtune/dq"] ,savename = "dtune.pdf")

ana_db, num_db =  dbi_dq(ring, ind_bpm, ind_Vquad, step)  