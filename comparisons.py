#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 13:13:15 2025

@author: deumenec

This file contains functions used to compare values of ORM matrices 
"""

import numpy as np

def rowComparison(matrix1, matrix2):
    """
    From two imput matrices, returns a vector of size m showing the difference by columns between the two matrices

    Parameters
    ----------
    matrix1 : np.array n*m matrix
        DESCRIPTION.
    matrix2 : np.array 
        DESCRIPTION.

    Returns
    -------
    None.

    """
    n = len(matrix1)
    m = len(matrix1[0])
    return [sum((matrix1[row][column]-matrix2[row][column])**2 for row in range(n) )for column in range(m)]

def columnComparison(matrix1, matrix2):
    """
    From two imput matrices, returns a vector of size n showing the difference by rows between the two matrices

    Parameters
    ----------
    matrix1 : np.array n*m matrix
        DESCRIPTION.
    matrix2 : np.array 
        DESCRIPTION.

    Returns
    -------
    None.

    """
    n = len(matrix1)
    m = len(matrix1[0])
    return [sum((matrix1[row][column]-matrix2[row][column])**2 for column in range(m) )for row in range(n)]


#Tests to see if it works out
matrix1 = np.array([[0,2], [2,3], [5, 6]])
matrix2 = np.array([[0,2], [2,3], [6, 7]])

print(columnComparison(matrix1, matrix2))
