#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 14:01:56 2021

@author: pascal
"""
import numpy as np
import CoolProp.CoolProp as CP

def get_alpha(T_boil, T, tol=1):
    #calculate heat transfer coefficient based on region (liquid, gaseous, two-phase)
    #inputs: float T_boil (boiling temperature (in K) of fluid at initial pressure)
    #        array T (N, temperature at all locations in the domain)
    #        float tol (tolerance around T_boil to assume two phase region, optional)
    #output: array alpha (N, heat transfer coefficent at all locations)
    alpha=np.zeros(T.size)
    for i in range(T.size):
        if type(T_boil)==type(None) or T[i]>T_boil+0.5*tol:
            alpha[i]=200
        elif T[i]<T_boil-0.5*tol:
            alpha[i]=2000
        else:
            alpha[i]=2500
    return alpha

def add_energy_source(res, T_boil, d, fields, T_wall):
    #calculate source term for energy equation
    #input: array res (3xN, array to store source term)
    #       float T_boil (boling temperature (in K) of fluid at initial pressure)
    #       float d (pipe diameter)
    #       array fields (6xN, containing all fields variables at all locations in the domain)
    #       array T_wall (N, wall temperature at all locations in the domain)
    #output:void, (3xN) array res (3xN, added source term in the last row) passed by reference
    for i in range(res[2,:].size):
        if fields[0,i]==0:
            res[2,i]=0
        else:
           alpha=get_alpha(T_boil, fields[4,:])
           res[2,:]= -alpha[:]*4/d*(fields[4,:]-T_wall[:]) 
    #alpha=get_alpha(T_boil, fields[4,:])
    #res[2,:]= -alpha[:]*4/d*(fields[4,:]-T_wall[:])
    
def add_momentum_source(res, fields, d, dx, nu):
    #calculate source term for momentum equation based on pipe resistance for laminar flows
    #currently not in use, because no change in result was observed
    #input: array res (3xN, array to store source term)
    #       array fields (6xN, containing all fields variables at all locations in the domain)
    #       float d (pipe diameter), dx (grid spacing), nu (kinematic viscosity)
    #output:void, (3xN) array res (3xN, added source term in the second row) passed by reference
    res[1,:]= -128*nu/(d**3)*dx*fields[1,:]
    
def add_momentum_source_2(res, fields, d, dx, nu, l=0.025):
    #calculate source term for momentum flows based on pipe resistance for turbulent flows
    #currently not in use, because no change in result was observed
    #input: array res (3xN, array to store source term)
    #       array fields (6xN, containing all fields variables at all locations in the domain)
    #       float d (pipe diameter), dx (grid spacing), nu (kinematic viscosity)
    #       float l (Darcy friction factor (Rohrreibungszahl) based on Moody Diagram), optional
    #output:void, (3xN) array res (3xN, added source term in the second row) passed by reference
    res[1,:]= -0.5*l*fields[1,:]**2/fields[0,:]/d