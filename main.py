# -*- coding: utf-8 -*-
"""
Created on Wen Nov 5 22:26:41 2025

@author: deumenec
"""

import at
import os
import numpy as np
import tqdm

calc_num_dORM = False
calc_analytic_ORM = True

def dORM_dqV (ring, dipoleIndex, ResponseClass, step, originalORM):
    ring[dipoleIndex].PolynomB += step
    ResponseClass.build_tracking()
    ring[dipoleIndex].PolynomB += -step  
    return (ResponseClass.response-originalORM)/step


def Cijn(phi_i, phi_j, nu, n):
    """
    shorthand for the cosine defined in the report
    """
    return(np.cos(n*abs(phi_i-phi_j)-n*np.pi*nu))
def Sijn(phi_i, phi_j, nu, n):
    """
    shorthand for the sine defined in the report
    """
    return(np.sign(phi_i-phi_j)*np.sin(n*abs(phi_i-phi_j)-n*np.pi*nu))

def dRij_dk(bi, bj, bk, Lk, nu, phi_i, phi_j, phi_k):
    """

    Parameters
    ----------
    bi : np.float 
        beta function in i 
    bj : np.float 
        beta function in j
    bk : TYPE
        beta function in k
    Lk : np.float 
        beta 
    nu : np.float 
        betatron tune of the ring
    phi_i : np.float 
        phase at the i BPM
    phi_j : np.float 
        phase at the j corrector
    phi_k : np.float 
        phase at the kth quadrupole

    Returns
    -------
    The orbit response matrix (i j) (BPM, corrector) coeficient in thin lens 
    approximation derivative with respect to the kth quadrupole

    """
    return np.sqrt(bi*bj)*bk*Lk/(8*np.sin(np.pi*nu)*np.sin(2*np.pi*nu))*(Cijn(phi_i, phi_j, nu, 1)*(Cijn(phi_i, phi_k,nu, 2)+Cijn(phi_j, phi_k,nu, 2)+2*np.cos(np.pi*nu)^2)+Sijn(phi_i, phi_j, nu, 1)*(Sijn(phi_i, phi_k,nu, 2 )-Sijn(phi_j, phi_k, 2)+np.sin(2*np.pi*nu)*(2*np.heaviside(phi_i-phi_k)-2*np.heaviside(phi_j-phi_k)-np.sign(phi_i-phi_j))))

def analyticdORM_dqV(ring, ind_bpm, ind_cor, ind_quad):
    """
    Parameters
    ----------
    ring : at.ring
        
    ind_bpm : list of indices for bpms

    ind_cor : list of indices for correctors
    
    ind_quad : list of indices for quadrupoles
    
    Returns
    -------
    the derivative of the ORM matrix computed analytically

    """
    bpmOptics = at.get_optics(ring, refpts=ind_bpm)
    corOptics = at.get_optics(ring, refpts=ind_cor)
    quadOptics = at.get_optics(ring, refpts=ind_quad)
    tune = bpmOptics[0]["mu"][2]
    bpmBeta = bpmOptics[2]["beta"][:][1]
    corBeta = corOptics[2]["beta"][:][1]
    quadBeta= quadOptics[2]["beta"][:][1]
    quadlength = [element.Length for element in ring[ind_quad]]
    bpmTunes= bpmOptics[2]["mu"][:][1]/(2*np.pi)
    corTunes = corOptics[2]["mu"][:][1]/(2*np.pi)
    quadTunes= quadOptics[2]["mu"][:][1]/(2*np.pi)
    
    dORM_dqV = []
    for quad in range(len())
    return dORM_dqV
    
os.chdir('/Users/deumenec/Documents/Uni/9é semestre/ALBA/Teoria/ORM_compare/') #Set my working directory!
ring = at.load_mat('THERING.mat')

ind_bpm = at.get_refpts(ring, 'BPM')   # family name, adjust to your lattice
ind_cor = at.get_refpts(ring, 'COR')  # horizontal correctors
ind_Vquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QV'))
ind_Hquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QH'))


if(calc_num_dORM == True):
    Rv = at.latticetools.OrbitResponseMatrix(ring, "v", ind_bpm, ind_cor) #class for computing the original ORM
    Rv.build_tracking()
    ORMv = Rv.response
    
    dORMijk = []
    for ind in tqdm.tqdm(ind_Vquad):
        dORMijk.append(dORM_dqV(ring, ind, Rv, 0.000001, ORMv))
    np.save("numdORM", np.array(dORMijk))
    
if(calc_num_dORM ==False):
    num_dORMijk = np.load("numdORM.npy")
    
if(calc_analytic_ORM==True):
    ana_dORMijk = analyticdORM_dqV(ring, ind_bpm, ind_cor, ind_Vquad)
    

