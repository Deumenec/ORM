# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 10:28:03 2025

@author: dhuerta
"""
def partialCheck(ring, ind_bpm, ind_cor, ind_quad, step):
    """
    numerical check of the partial derivative formulas deduced to try and spot any possible mistake.
    """
    bpmOptics = at.get_optics(ring, refpts=ind_bpm)
    corOptics = at.get_optics(ring, refpts=ind_cor)
    quadOptics = at.get_optics(ring, refpts=ind_quad)
    tune = [bpmOptics[1]["tune"][1]]
    bpmBeta = [[i[1] for i in bpmOptics[2]["beta"]]]
    corBeta = [[i[1] for i in corOptics[2]["beta"]]]
    quadBeta= [[i[1] for i in quadOptics[2]["beta"]]]
    quadLen = [[element.Length for element in ring[ind_quad]]]
    bpmTune = [[i[1] for i in bpmOptics[2]["mu"]]] #Important, mu doesn't have the /2pi factor in atcollab!
    corTune = [[i[1] for i in corOptics[2]["mu"]]]
    quadTune= [[i[1] for i in quadOptics[2]["mu"]]]
    ring[ind_quad[0]].PolynomB += step
    bpmOptics = at.get_optics(ring, refpts=ind_bpm)
    corOptics = at.get_optics(ring, refpts=ind_cor)
    quadOptics = at.get_optics(ring, refpts=ind_quad)
    tune += [bpmOptics[1]["tune"][1]]
    bpmBeta += [[i[1] for i in bpmOptics[2]["beta"]]]
    corBeta += [[i[1] for i in corOptics[2]["beta"]]]
    quadBeta+= [[i[1] for i in quadOptics[2]["beta"]]]
    quadLen += [[element.Length for element in ring[ind_quad]]]
    bpmTune += [[i[1] for i in bpmOptics[2]["mu"]]] #Important, mu doesn't have the /2pi factor in atcollab!
    corTune += [[i[1] for i in corOptics[2]["mu"]]]
    quadTune+= [[i[1] for i in quadOptics[2]["mu"]]]
    tune, bpmBeta, corBeta, bpmTune, corTune = np.array(tune), np.array(bpmBeta), np.array(corBeta), np.array(bpmTune), np.array(corTune)
    num_d_dq = {"tune": (tune[1]-tune[0])/step,
                "bpmBeta": (bpmBeta[1]- bpmBeta[0])/step,
                "corBeta": (corBeta[1]-corBeta[0])/step,
                "bpmTune": (bpmTune[1]-bpmTune[0])/step,
                "corTune": (corTune[1]-corTune[0])/step}
    ana_d_dq = {"tune":     -quadBeta[0][0]*quadLen[0][0]/(4*np.pi),
                "bpmBeta": (bpmBeta[0]*quadBeta[0][0]*quadLen[0][0]/2/np.sin(2*np.pi*tune[0])),
                "corBeta": (bpmBeta[0]*quadBeta[0][0]*quadLen[0][0]/2/np.sin(2*np.pi*tune[0])),
                #"bpmTune": quadBeta[0][0]*quadLen[0][0]/(4*np.sin(2*np.pi*tune[0]))*(Sijn(bpmTune[0], quadTune[0], tune, 2)+np.sin(2*quadTune[0]-2*np.pi*tune[0])+2*np.heaviside(corTune[0]-quadTune[0])*np.sin(2*np.pi*tune)),
                #"corTune": quadBeta[0][0]*quadLen[0][0]/(4*np.sin(2*np.pi*tune[0]))*(Sijn(corTune[0], quadTune[0], tune, 2)+np.sin(2*quadTune[0]-2*np.pi*tune[0])+2*np.heaviside(corTune[0]-quadTune[0])*np.sin(2*np.pi*tune))
                }
    return

partialCheck(ring, ind_bpm, ind_cor, ind_Vquad, 0.00001)