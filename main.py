# -*- coding: utf-8 -*-
"""
Created on Wen Nov 5 22:26:41 2025

@author: deumenec
"""

import at
import os
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import comparisons

calc_num_dORM = False
calc_analytic_dORM = True

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
    Shorthand for the sine defined in the report
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
    return np.sqrt(bi*bj)*bk*Lk/(8*np.sin(np.pi*nu)*np.sin(2*np.pi*nu))*(Cijn(phi_i, phi_j, nu, 1)*(Cijn(phi_i, phi_k,nu, 2)+Cijn(phi_j, phi_k,nu, 2)+2*np.cos(np.pi*nu)**2)+Sijn(phi_i, phi_j, nu, 1)*(Sijn(phi_i, phi_k,nu, 2 )-Sijn(phi_j, phi_k,nu, 2)+np.sin(2*np.pi*nu)*(2*np.heaviside(phi_i-phi_k, 1)-2*np.heaviside(phi_j-phi_k, 1)-np.sign(phi_i-phi_j))))

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
    tune = bpmOptics[1]["tune"][1]
    bpmBeta = [i[1] for i in bpmOptics[2]["beta"]]
    corBeta = [i[1] for i in corOptics[2]["beta"]]
    quadBeta= [i[1] for i in quadOptics[2]["beta"]]
    quadLen = [element.Length for element in ring[ind_quad]]
    bpmTune = [i[1] for i in bpmOptics[2]["mu"]] #Important, mu doesn't have the /2pi factor in atcollab!
    corTune = [i[1] for i in corOptics[2]["mu"]]
    quadTune= [i[1] for i in quadOptics[2]["mu"]]
    dORM_dqV = []
    for k in range(len(ind_quad)):
        dORM_dqV.append([])
        for i in range(len(ind_bpm)):
            dORM_dqV[k].append([])
            for j in range(len(ind_cor)):
                dORM_dqV[k][i].append(dRij_dk(bpmBeta[i], corBeta[j], quadBeta[k], quadLen[k], tune, bpmTune[i], corTune[j], quadTune[k]))
    return np.array(dORM_dqV)
    
def analyticORM(ring, ind_bpm, ind_cor):
    """
    Parameters
    ----------
    ring : TYPE
        DESCRIPTION.
    ind_bpm : TYPE
        DESCRIPTION.
    ind_cor : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    bpmOptics = at.get_optics(ring, refpts=ind_bpm)
    corOptics = at.get_optics(ring, refpts=ind_cor)
    tune = bpmOptics[1]["tune"][1]
    bpmBeta = [i[1] for i in bpmOptics[2]["beta"]]
    corBeta = [i[1] for i in corOptics[2]["beta"]]
    bpmTune = [i[1] for i in bpmOptics[2]["mu"]]
    corTune = [i[1] for i in corOptics[2]["mu"]]
    ORM = []
    for i in range(len(ind_bpm)):
        ORM.append([])
        for j in range(len(ind_cor)):
            ORM[i].append(np.sqrt(bpmBeta[i]*corBeta[j])/(2*np.sin(np.pi*tune))*np.cos(abs(bpmTune[i]-corTune[j])-np.pi*tune))
    return np.array(ORM)
def dimSSD(matrix1, matrix2, dim):
    """Sum of square differences along que specified dimension"""
    return np.sum((matrix1-matrix2)**2, axis =dim)

    
os.chdir('Z:\Projectes\ORM\ORM') #Set my working directory!
ring = at.load_mat('THERING.mat')

ind_bpm = at.get_refpts(ring, 'BPM')   # family name, adjust to your lattice
ind_cor = at.get_refpts(ring, 'COR')  # horizontal correctors
ind_Vquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QV'))
ind_Hquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QH'))


Rv = at.latticetools.OrbitResponseMatrix(ring, "v", ind_bpm, ind_cor) #class for computing the original ORM
Rv.build_tracking()
ORMv = Rv.response

if(calc_num_dORM == True):    
    dORMijk = []
    for ind in tqdm.tqdm(ind_Vquad):
        dORMijk.append(dORM_dqV(ring, ind, Rv, 0.000001, ORMv))
    np.save("numdORM", np.array(dORMijk))
    
if(calc_num_dORM ==False):
    num_dORMijk = np.load("numdORM.npy")
    
if(calc_analytic_dORM==True):
    ana_dORMijk = analyticdORM_dqV(ring, ind_bpm, ind_cor, ind_Vquad)
    np.save("anadORM", ana_dORMijk)
    
if(calc_num_dORM ==False):
    ana_dORMijk = np.load("anadORM.npy")

ana_ORMv = analyticORM(ring, ind_bpm, ind_cor)

dades1= 100000*dimSSD(ORMv, ana_ORMv, 1)
dades2= dimSSD(ORMv, np.zeros([120, 88]), 1)
dades3= dimSSD(ana_ORMv, np.zeros([120, 88]), 1)

plt.plot(dades1, color = "red", label = "10^5 ORM difference")
plt.plot(dades2, color = "blue", label = "analytical ORM")
plt.plot(dades3,color ="green", linestyle = "--", label = "numerical ORM")
plt.xlabel("BPM index", loc = "right")
plt.legend()

#plt.savefig("ORMSSDBPMS.pdf")


