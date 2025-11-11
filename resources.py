 # -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 10:30:53 2025

@author: dhuerta

In this code I implement advices by Zeus to better compute numerically the CS alpha from beta

"""


import os
import at
import matplotlib.pyplot as plt
import inspect
import numpy as np
import copy

os.chdir('/Users/deumenec/Documents/Uni/9é semestre/ALBA/Teoria/ORM_compare/') #Set my working directory!


def pointPosition (sPoints: np.array,ringGeometry: np.array):
    """
    Parameters
    ----------
    s : np.array
        List of ordered s points
    ringGeometry : np.array
        Geometry of a ring with the start coordinate of each element

    Returns
    -------
    [int, float]
    The element number and length inside of that element where eacch point is
    inside the ring as an array of [elementNum, length]
    """
    sInsidePositions = []
    positionScan = 0
    for s in sPoints:
        if positionScan < len(ringGeometry):
            while (ringGeometry[positionScan]< s):
                positionScan +=1
        sInsidePositions.append([positionScan-1,s-ringGeometry[positionScan-1]])
    return sInsidePositions

def opticalFunctions(ring, ringOptics, point, step):
    """
    Parameters
    ----------
    ring : lattice.lattice_object.Lattice
        ring where optical parameters are calculated
    ringOptics : tuple
        as returned by at.get_optics containing the twiss functions in each element
    point : [int, float]
        point inside an element where the optical parameters are calculated
    step : float
        DESCRIPTION.

    Returns
    -------
    An array and a 1 or zero depending on if the calculation could be performed.
    AT beta, AT alpha, numerical alpha through the at beta and a mask for points were alpha couldn't be calculated.
    For simplicity we only consider the x optical funcitons

    """
    #
    elementLength= ring[point[0]].Length
    if ( elementLength<=point[1]+step):
        print("Fail at element", point[0],"at the lenght" ,point[1])
        return [[], 0]
    #Create a copy of the segment of the ring with the element
    ringSegment = ring[point[0]:point[0]+1]
    newElements= ringSegment[0].divide([point[1]/elementLength, step/elementLength, (elementLength-point[1]-step)/elementLength])
    ringSegment[0] = newElements[0]
    ringSegment.append(newElements[1])
    ringSegment.append(newElements[2])
    #Calculate the optics along that elemment based on the initial ones
    segmentOptics = at.physics.linopt6(ringSegment, twiss_in=ringOptics[point[0]], refpts=[1,2])
    #We calculate the alpha numerically through beta before and after step
    numAlpha = -1/2*(segmentOptics[2].beta[1][0]-segmentOptics[2].beta[0][0])/step
    return [[segmentOptics[2].beta[0][0],segmentOptics[2].alpha[0][0], numAlpha], 1]



def plotOptics (sPoints, opticalFunctions):
    """
    Parameters
    ----------
    sPoints          : points where the opticalFunctions have been atempted to calculate
    opticalFunctions : tuple
        optical functions as returned by opticalFunctions

    Returns
    -------
    None. prints the matplotlib plot with the data

    """
    atBetas      = []
    atAlphas     = []
    numAlphas   = []
    plotPoints  = []

    for i in range(len(sPoints)):
        if (opticalFunctions[i][1]==1):
            atBetas.append(opticalFunctions[i][0][0])
            atAlphas.append(opticalFunctions[i][0][1])
            numAlphas.append(opticalFunctions[i][0][2])
            plotPoints.append(sPoints[i])
    plt.plot(plotPoints, atBetas,label= "atBetas" )
    plt.plot(plotPoints, atAlphas,label= "atAlphas" )
    plt.plot(plotPoints, numAlphas,label= "numAlphas", linestyle ='--')
    plt.legend()
    plt.savefig("optical_functions.pdf")


ring = at.load_lattice('THERING.mat')
ringGeometry = at.lattice.get_geometry(ring)[0] #Start coordinate of each element along the s coordinate

plotInterval    = [0.1,26]        #Interval where the functions are calculated
nPoints         = 2000       #Number of points where alpha and beta are calculated
derivativeStep  = 0.01          #Step used to calculate the derivative
sPoints         = np.linspace(plotInterval[0], plotInterval[1], num = nPoints)
ringOptics      = at.get_optics(ring, refpts=range(0,len(ring))) #Array with the optics at each component

insidePos = pointPosition(sPoints, ringGeometry.x)

computedOpticalFunctions = []

for point in insidePos:
    computedOpticalFunctions.append(opticalFunctions(ring, ringOptics[2], point, derivativeStep))

plotOptics(sPoints, computedOpticalFunctions)

def main():
    return
if __name__ =="__main__":
    main()