 # -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 10:30:53 2025

@author: dhuerta

In this code I implement advices by Zeus to better compute numerically the CS alpha from beta

"""
import numpy as np
import at
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

corrector_optics = at.get_optics(ring, refpts=range(0,len(ring)))
    