# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 10:30:02 2025

@author: dhuerta
In this file, I test the formulas for partial derivatives of optical functions
with respect to quadrupole strength provided in the report
"""

import at
import os
import numpy as np
import tqdm
import matplotlib.pyplot as plt

os.chdir('Z:\Projectes\ORM\ORM') #Set my working directory!
ring = at.load_mat('THERING.mat')
step = 0.00001
index = 9

ind_bpm = at.get_refpts(ring, 'BPM')   # family name, adjust to your lattice
ind_cor = at.get_refpts(ring, 'COR')  # horizontal correctors
ind_Vquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QV'))
ind_Hquad = at.get_refpts(ring, lambda el: el.FamName.startswith('QH'))

optics1 = at.get_optics(ring, refpts= [ind_bpm[0], ind_cor[index]])
ring[ind_Vquad[index]].PolynomB[1]+=step
optics2 = at.get_optics(ring, refpts= [ind_bpm[0], ind_cor[index]])

dTunedq_ana= optics1[2]["beta"][1][1]*ring[ind_Vquad[index]].Length/(4*np.pi)
dTunedq_num= (optics1[1]["tune"][1]-optics2[1]["tune"][1])/step
 
dTuneerror = abs((dTunedq_ana-dTunedq_num)/dTunedq_num)
print(f"In the thin lens aproximation, the tune derivative displays an error of around {dTuneerror*100}%")

