#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 14:01:56 2021

@author: pascal
"""
import numpy as np
import CoolProp.CoolProp as CP

def get_alpha(T_boil, T, tol=1):
    alpha=np.zeros(T.size)
    for i in range(T.size):
        if T[i]>T_boil+0.5*tol:
            alpha[i]=200
        elif T[i]<T_boil-0.5*tol:
            alpha[i]=2000
        else:
            alpha[i]=2300
    return alpha

def add_energy_source(res, T_boil, d, fields, R, T_wall, fluid):
    #calculate source termsource
    #input:
    #output: (4xN) array with added  term, passed by reference
    alpha=get_alpha(T_boil, fields[4,:])
    res[2,:]= -alpha[:]*4/d*(fields[4,:]-T_wall[:])
    
def add_momentum_source(res, fields, d, dx, nu):
    #calculate pipe resistance
    #input:
    #output:
    res[1,:]= -128*nu/(d**3)*dx*fields[1,:]
    
def add_momentum_source_2(res, fields, d, dx, nu, l=0.025):
    #calculate pipe resistance
    #input:
    #output:
    res[1,:]= -0.5*l*fields[1,:]**2/fields[0,:]/d