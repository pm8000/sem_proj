#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 10:16:43 2021

@author: pascal
"""

import numpy as np
import rs_trial_EOS as eos
import pandas as pd

def minmod(u, u_p, u_m):
    #apply minmod slope limiter
    #input: scalars value of where slope is to be determined u, value of cell right u_p, value of cell left r_m
    #output: scalar slope inside cell
    d1=u_p-u
    d2=u-u_m
    if np.sign(d1)!=np.sign(d2):
        return 0
    return np.sign(d1)*min(abs(d1),abs(d2))

def get_middle_state(fields, i, inlet=[0,0,0], border=0):
    #calculate the values at cell interface assuming linearised problem at boundary and use linear profile within cell limited by minmod
    #input: (4xN) array with stored values of density, velocity and pressure, integer i to denote index,
    #       vector of length 3, containing (in this order) inlet density, velocity and pressure, only required if i==1
    #       parameter border, set to 1 if left cell is boundry cell, -1 if right cell is boundary cell
    #output: tuple of middle state values for density, velocity, pressure and energy
    if (border !=0 and border !=1 and border !=-1):
        print("encountered unknown value for border in function get_middle_state")
        assert(False)
    #get values at cell boundary
    #boundary cases
    if border==1:
        #construct ghost cell to the left to achieve desired flux
        rho_l=fields[0,i-1]+0.5*minmod(fields[0, i-1], fields[0,i], 2*inlet[0]-fields[0,i-1])
        u_l=fields[1,i-1]/fields[0, i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], 2*inlet[1]-fields[1,i-1]/fields[0,i-1])
        p_l=fields[3,i-1]+0.5*minmod(fields[3, i-1], fields[3,i], 2*inlet[2]-fields[3,i-1])
        
    else:
        rho_l=fields[0, i-1]+0.5*minmod(fields[0,i-1], fields[0,i], fields[0,i-2])
        u_l=fields[1,i-1]/fields[0,i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], fields[1,i-2]/fields[0,i-1])
        p_l=fields[3,i-1]+0.5*minmod(fields[3,i-1], fields[3,i], fields[3,i-2])
    if rho_l<0 or pd.isna(rho_l)==True:
        print("density is negative or not a number")
        assert(False)
    cs_l=eos.get_cs(p_l, rho_l)
    
    if border==-1:
        #on right border there is zero gradient BC, hence minmod would return 0 slope
        rho_r=fields[0,i]
        u_r=fields[1,i]/rho_r
        p_r=fields[3,i]
    else:
        rho_r=fields[0,i]+0.5*minmod(fields[0,i], fields[0,i+1], fields[0,i-1])
        u_r=fields[1,i]/fields[0,i]+0.5*minmod(fields[1,i]/fields[0,i], fields[1,i+1]/fields[0,i+1], fields[1,i-1]/fields[0,i-1])
        p_r=fields[3,i]+0.5*minmod(fields[0,i], fields[3,i+1], fields[3,i-1])
    cs_r=eos.get_cs(p_r, rho_r)
        
    #compute middle state with linearised characteristics
    #assume always subsonic
        
    u_m=(rho_l*u_l*cs_l+rho_r*u_r*cs_r+p_l-p_r)/(rho_l*cs_l+rho_r*cs_r)
    p_m=p_l-rho_l*cs_l*(u_m-u_l)
    if u_l+u_r>=0:
        rho_m=(p_m-p_l)/(cs_l*cs_l)+rho_l
    else:
        rho_m=(p_m-p_r)/(cs_r*cs_r)+rho_r
    E_m=eos.get_E(rho_m, u_m, p_m)
    return (rho_m, u_m, p_m, E_m)

def get_outlet_bc(fluxes, fields):
    #calculate outflow flux based on zero gradient
    #input: (4xN) array fields with stored values of density, momentum, energy and pressure, (3xN+1) array fluxes to write fluxes to
    #output:(3xN+1) array with updated fluxes, passed by reference
    
    fluxes[0,-1]=fields[1,-1]
    fluxes[1,-1]=fields[1,-1]**2/fields[0,-1]+fields[3,-1]
    fluxes[2,-1]=fields[1,-1]/fields[0,-1]*(fields[2,-1]+fields[3,-1])

def calc_inlet_bc(rho, u, p):
    #calculate inlet fluxes given by BC. Fluxes are assumed to be constant, hence only one evaluation is required
    #input: scalars inlet density, velocity and pressure
    #output: (1x3) vector with inlet fluxes
    res=np.zeros([3])
    res[0]=rho*u
    res[1]=rho*u*u+p
    res[2]=u*(eos.get_E(rho, u, p)+p)
    return res

def get_inflow_bc(fluxes, inlet):
    #add the inlet fluxes
    #input: (3xN+1) array to write inlet fluxes to, (1x3) vector containing fluxes
    #output: inlet fluxes added to fluxes in the first column, passed by reference
    fluxes[0,0]=inlet[0]
    fluxes[1,0]=inlet[1]
    fluxes[2,0]=inlet[2]
            
def update_fluxes(fluxes, fields, inlet_flux, inlet):
    #calculate fluxes at cell interfaces, assuming Riemann problems at cell boundaries
    #input: (4xN) array fields with stored values of density, momentum, energy and pressure, (3xN+1) array fluxes to write fluxes to
    #       vector of length 3 (inlet_flux) containing inlet fluxes
    #       vector of length 3 (inlet) containing (in this order) inlet density, velocity and pressure
    #output: (3xN+1) array with updated fluxes, passed by reference
    
    get_inflow_bc(fluxes, inlet_flux)
    
    for i in range(1, fluxes.shape[1]-1):
        if i==1:
            rho_m, u_m, p_m, E_m=get_middle_state(fields, i, inlet=inlet, border=1)
        elif i==fluxes.shape[1]-2:
            rho_m, u_m, p_m, E_m=get_middle_state(fields, i, border=-1)
        else:
            rho_m, u_m, p_m, E_m=get_middle_state(fields, i)
        fluxes[0,i]=rho_m*u_m
        fluxes[1,i]=rho_m*u_m*u_m+p_m
        fluxes[2,i]=u_m*(E_m+p_m)
        
    get_outlet_bc(fluxes, fields)
    
def update_fluxes_upwind(fluxes, fields, inlet_flux, inlet):
    #calculate flux at cell interface based on upwind scheme
    get_inflow_bc(fluxes, inlet_flux)
    for i in range(1, fluxes.shape[1]-1):
        if fields[1,i-1]/fields[0,i-1]+fields[1,i]/fields[0,i]<0:
            fluxes[0,i]=fields[1,i]
            fluxes[1,i]=fields[1,i]**2/fields[0,i]+fields[3,i]
            fluxes[2,i]=fields[1,i]/fields[0,i]*(fields[2,i]+fields[3,i])
        else:
            fluxes[0,i]=fields[1,i-1]
            fluxes[1,i]=fields[1,i-1]**2/fields[0,i-1]+fields[3,i-1]
            fluxes[2,i]=fields[1,i-1]/fields[0,i-1]*(fields[2,i-1]+fields[3,i-1])
    get_outlet_bc(fluxes, fields)