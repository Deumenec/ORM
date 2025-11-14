# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 10:30:02 2025

@author: dhuerta
In this file, I test the formulas for all partial derivatives of optical functions
with respect to quadrupole strength provided in the report, thus I can quantify
how bad with respect to each component the linear aproximation is.

In derivatives, the derivative with respect to que quadrupolar strength is the
first component.
Comprovem que el mètode, encara que comporta cert error, és correcte!
El problema abans sembla que estava en la manera de la qual calculava la 
"""

import at
import os
import numpy as np
import tqdm
import matplotlib.pyplot as plt

def vectorPlot(vectors, vecNames, savename):
    """ Creates a plot commparing the components of vectors of equal length"""
    colors = plt.cm.rainbow(np.linspace(0, 1, len(vectors)))
    for i, vector in enumerate(vectors):
        plt.plot(vector, color = colors[i], label= vecNames[i])
    plt.legend()
    plt.savefig("dORMs_comparison.pdf", format = "pdf")
    plt.show() #Això ha de ser l'últim perqué es peta la figura
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
            result_ana[i][j] = el_beta[j]*quad_beta[i]*quadLen[i]/(2*np.sin(2*np.pi*tune))*Cijn(el_psi[j], quad_psi[i], tune, 2)
    return result_ana, result_num

def dpsi_dq(ring, ind_el, ind_quad, step):
    """ Returns the derivative of the phase in elements with respect to change along each quadrupole strength by computing the change in optics"""
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
            result_num[i][j] = (optics1el[2]["mu"][j][1] -el_psi[j] )/step
            #Here I calculate it with the formula
            result_ana[i][j] = -quad_beta[i]*quadLen[i]/(4*np.sin(2*np.pi*tune))*(Sijn(el_psi[j], quad_psi[i], tune, 2)+np.sin(2*quad_psi[i]-2*np.pi*tune)+2*np.heaviside(el_psi[j]-quad_psi[i],1)*np.sin(2*np.pi*tune))
    return result_ana, result_num

def analyticORM(ring, ind_bpm, ind_cor):
    """Calculates analytically the ORM (this one has been checked to work)"""
    bpmOptics = at.get_optics(ring, refpts=ind_bpm)
    corOptics = at.get_optics(ring, refpts=ind_cor)
    tune = bpmOptics[1]["tune"][1]
    bpmBeta = [i[1] for i in bpmOptics[2]["beta"]]
    corBeta = [i[1] for i in corOptics[2]["beta"]]
    bpmTune = [i[1] for i in bpmOptics[2]["mu"]]
    corTune = [i[1] for i in corOptics[2]["mu"]]
    ORM = np.zeros([len(ind_bpm), len(ind_cor)])
    for i in range(len(ind_bpm)):
        for j in range(len(ind_cor)):
            ORM[i][j] = np.sqrt(bpmBeta[i]*corBeta[j])/(2*np.sin(np.pi*tune))*np.cos(abs(bpmTune[i]-corTune[j])-np.pi*tune)
    return np.array(ORM)

def dORMdqi(ring, ind_bpm, ind_cor, dbbpm_dqk, dbcor_dqk, dnu_dqk, dpsibpm_dqk, dpsicor_dqk):
    """Calculates the derivative of the ORM with respect to dqi Jacobian using partial derivatives comuted before"""
    #We begin initializing the required variables
    bpmOptics = at.get_optics(ring, refpts=ind_bpm)
    corOptics = at.get_optics(ring, refpts=ind_cor)
    tune = bpmOptics[1]["tune"][1]
    bpmBeta = [i[1] for i in bpmOptics[2]["beta"]]
    corBeta = [i[1] for i in corOptics[2]["beta"]]
    bpmTune = [i[1] for i in bpmOptics[2]["mu"]]
    corTune = [i[1] for i in corOptics[2]["mu"]]
    
    result = np.zeros([len(dbcor_dqk), len(ind_bpm), len(ind_cor)])
    for k in range(len(dbcor_dqk)):
        for i in range(len(ind_bpm)):
            for j in range(len(ind_cor)):
                #I calculate separatelly the 5 terms given by the chain rule
                dbi_term   = np.sqrt(corBeta[j]/bpmBeta[i])/(4*np.sin(np.pi*tune))*Cijn(bpmTune[i], corTune[j], tune, 1)*dbbpm_dqk[k][i]
                dbj_term   = np.sqrt(bpmBeta[i]/corBeta[j])/(4*np.sin(np.pi*tune))*Cijn(bpmTune[i], corTune[j], tune, 1)*dbcor_dqk[k][j]
                dtune_term = -np.pi*np.sqrt(bpmBeta[i]*corBeta[j])/(2*np.sin(np.pi*tune)**2)*(Cijn(bpmTune[i], corTune[j], tune, 1)-np.sign(bpmTune[i]-corTune[j])*Sijn(bpmTune[i], corTune[j], tune, 1)*np.sin(np.pi*tune))*dnu_dqk[k]
                dpsii_term = -np.sqrt(bpmBeta[i]*corBeta[j])/(2*np.sin(np.pi*tune))*Sijn(bpmTune[i], corTune[j], tune, 1)*dpsibpm_dqk[k][i]
                dpsij_term = np.sqrt(bpmBeta[i]*corBeta[j])/(2*np.sin(np.pi*tune))*Sijn(bpmTune[i], corTune[j], tune, 1)*dpsicor_dqk[k][j]
                result[k][i][j]= dbi_term + dbj_term + dtune_term + dpsii_term + dpsij_term
    return result

#os.chdir('Z:\Projectes\ORM\ORM') #Set my working directory!
os.chdir('/Users/deumenec/Documents/Uni/9é semestre/ALBA/Teoria/ORM_compare') #Set my working directory!
ring = at.load_mat('THERING.mat')
step = 0.0000001


ind_bpm = at.get_refpts(ring, 'BPM')   # family name, adjust to your lattice
ind_cor = at.get_refpts(ring, 'COR')  # horizontal correctors
ind_Vquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QV'))
ind_Hquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QH'))

#ana_dnu, num_dnu = dnu_dq(ring, ind_Vquad, step)
#vectorPlot([ana_dnu, num_dnu], ["Analytical dtune/dq", "Numerical dtune/dq"] ,savename = "dtune.pdf")

#ana_db, num_db =  dbi_dq(ring, ind_bpm, ind_Vquad, step)
#vectorPlot([ana_db[10], num_db[10]], ["Analytic db/dq", "Numerical db/dq"], savename = "db for k=10 in BPMs.pdf")

#ana_dpsi, num_dpsi = dpsi_dq(ring, ind_bpm, ind_Vquad, step)
#vectorPlot([ana_dpsi[20], num_dpsi[20]], ["Analytic dpsi/dq", "Numerical dpsi/dq"], savename = "dpsi for k=20 in BPMs.pdf")

ORM = analyticORM(ring, ind_bpm, ind_cor)
ring[ind_Vquad[0]].PolynomB[1] +=step
ORM2 = analyticORM(ring, ind_bpm, ind_cor)
ring[ind_Vquad[0]].PolynomB[1] -=step

a, dbbpm_dqk  = dbi_dq(ring, ind_bpm, ind_Vquad, step)
a, dbcor_dqk  = dbi_dq(ring, ind_cor, ind_Vquad, step)
a, dnu_dqk    = dnu_dq(ring, ind_Vquad, step)
a, dpsibpm_dqk= dpsi_dq(ring, ind_bpm, ind_Vquad, step)
a, dpsicor_dqk= dpsi_dq(ring, ind_cor, ind_Vquad, step)

dORM = dORMdqi(ring, ind_bpm, ind_cor, dbbpm_dqk, dbcor_dqk, dnu_dqk, dpsibpm_dqk, dpsicor_dqk)
dORM_num = (ORM2-ORM)/step

Rv = at.latticetools.OrbitResponseMatrix(ring, "v", ind_bpm, ind_cor) #class for computing the original ORM
Rv.build_tracking()
ORM = Rv.response
ring[ind_Vquad[0]].PolynomB[1] +=step
Rv.build_tracking()
ORM2 = Rv.response
ring[ind_Vquad[0]].PolynomB[1] -=step
dORM_num_num = (ORM2-ORM)/step

vectorPlot([dORM[0][10], dORM_num[10], dORM_num_num[10]], ["Analytic ORM", "Numerical ORM", "Numerical Numerical ORM"], savename = "ORM_0_10_comparison.pdf")
