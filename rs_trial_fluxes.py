#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 10:16:43 2021

@author: pascal
"""

import numpy as np
import rs_trial_EOS as eos
import CoolProp.CoolProp as CP
import pandas as pd
import rs_trial_parallelisation as para

def minmod(u, u_p, u_m):
    #apply minmod slope limiter
    #input: scalars value of where slope is to be determined u, value of cell right u_p, value of cell left r_m
    #output: scalar slope inside cell
    d1=u_p-u
    d2=u-u_m
    if np.sign(d1)!=np.sign(d2):
        return 0
    return np.sign(d1)*min(abs(d1),abs(d2))

def get_middle_state(u_l, rho_l, p_l, u_r, rho_r, p_r, fluid):
    #calculate middle state based on characteristics of left and right state:
    #inputs:left and right state (all scalar values)
    #output: middle state (tuple of scalar values)
    chi_l=eos.get_chi(p_l, rho_l, fluid)
    chi_r=eos.get_chi(p_r, rho_r, fluid)
    
    if u_l>chi_l:
        print("supersonic velocity encountered")
        u_m=u_l
        p_m=p_l
        rho_m=rho_l
    else:
        u_m=(rho_l*u_l*chi_l+rho_r*u_r*chi_r+p_l-p_r)/(rho_l*chi_l+rho_r*chi_r)
        p_m=p_l-rho_l*chi_l*(u_m-u_l)
        if u_l+u_r>=0:
            rho_m=(p_m-p_l)/(chi_l*chi_l)+rho_l
        else:
            rho_m=(p_m-p_r)/(chi_r*chi_r)+rho_r
    E_m=eos.get_E(rho_m, u_m, p_m, fluid)
    return (rho_m, u_m, p_m, E_m)

def get_border_values(fields, i, fluid, inlet=[0,0,0], border=0, p_atm=None, border_val=None):
    #calculate the values at cell interface assuming linear profile within cell limited by minmod
    #input: (4xN) array with stored values of density, velocity and pressure, integer i to denote index,
    #       vector of length 3, containing (in this order) inlet density, velocity and pressure, only required if i==1
    #       parameter border, set to 1 if left cell is boundry cell, -1 if right cell is boundary cell
    #output: tuple of middle state values for density, velocity, pressure and energy
    if (border <-3 or border >2):
        print("encountered invalid value for border in function get_middle_state")
        assert(False)
    #get values at cell boundary
    #boundary cases
    if border==1:
        #construct ghost cell to the left to achieve desired flux
        rho_l=fields[0,i-1]+0.5*minmod(fields[0, i-1], fields[0,i], 2*inlet[0]-fields[0,i-1])
        u_l=fields[1,i-1]/fields[0, i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], 2*inlet[1]-fields[1,i-1]/fields[0,i-1])
        p_l=fields[3,i-1]
        #p_l=fields[3,i-1]+0.5*minmod(fields[3, i-1], fields[3,i], 2*fields[3,0]-fields[3,i-1])
    
    elif border==2:
        rho_l=fields[0, i-1]+0.5*minmod(fields[0,i-1], fields[0,i], border_val[0])
        u_l=fields[1,i-1]/fields[0,i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], border_val[1])
        p_l=fields[3,i-1]+0.5*minmod(fields[3,i-1], fields[3,i], border_val[2])
        
    elif border==-3:
        rho_l=fields[0, i-1]+0.5*minmod(fields[0,i-1], border_val[0,0], fields[0,i-2])
        u_l=fields[1,i-1]/fields[0,i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], border_val[1,0], fields[1,i-2]/fields[0,i-2])
        p_l=fields[3,i-1]+0.5*minmod(fields[3,i-1], border_val[2,0], fields[3,i-2])
        
    else:
        rho_l=fields[0, i-1]+0.5*minmod(fields[0,i-1], fields[0,i], fields[0,i-2])
        u_l=fields[1,i-1]/fields[0,i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], fields[1,i-2]/fields[0,i-2])
        p_l=fields[3,i-1]+0.5*minmod(fields[3,i-1], fields[3,i], fields[3,i-2])
    if rho_l<0 or pd.isna(rho_l)==True:
        print("density is negative or not a number")
        assert(False)
    
    if border==-1:
        #on right border there is zero gradient BC, hence minmod would return 0 slope
        rho_r=fields[0,i]
        u_r=fields[1,i]/rho_r
        p_r=fields[3,i]-0.5*minmod(fields[3, i-1], fields[3,i], 2*p_atm-fields[3,i-1])
        #p_r=fields[3,i]
        
    elif border==-2:
        #on right border process border shared values are taken
        rho_r=fields[0,i]-0.5*minmod(fields[0,i], border_val[0,0], fields[0,i-1])
        u_r=fields[1,i]/fields[0,i]-0.5*minmod(fields[1,i]/fields[0,i], border_val[1,0], fields[1,i-1]/fields[0,i-1])
        p_r=fields[3,i]-0.5*minmod(fields[0,i], border_val[2,0], fields[3,i-1])

    elif border==-3:
        rho_r=border_val[0,0]-0.5*minmod(border_val[0,0], border_val[0,1], fields[0,i-1])
        u_r=border_val[1,0]-0.5*minmod(border_val[1,0], border_val[1,1], fields[1,i-1]/fields[0,i-1])
        p_r=border_val[2,0]-0.5*minmod(border_val[2,0], border_val[2,1], fields[3,i-1])
    
    else:
        rho_r=fields[0,i]-0.5*minmod(fields[0,i], fields[0,i+1], fields[0,i-1])
        u_r=fields[1,i]/fields[0,i]-0.5*minmod(fields[1,i]/fields[0,i], fields[1,i+1]/fields[0,i+1], fields[1,i-1]/fields[0,i-1])
        p_r=fields[3,i]-0.5*minmod(fields[0,i], fields[3,i+1], fields[3,i-1])
        
    return(get_middle_state(u_l, rho_l, p_l, u_r, rho_r, p_r, fluid))

def zero_gradient_outlet_bc(fluxes, fields, E_correction=0):
    #calculate outflow flux based on zero gradient
    #input: (4xN) array fields with stored values of density, momentum, energy and pressure, (3xN+1) array fluxes to write fluxes to
    #output:(3xN+1) array with updated fluxes, passed by reference
    
    fluxes[0,-1]=fields[1,-1]
    fluxes[1,-1]=fields[1,-1]**2/fields[0,-1]+fields[3,-1]
    fluxes[2,-1]=fields[1,-1]/fields[0,-1]*(fields[2,-1]+fields[3,-1])
    
def no_outlet(fluxes, fields, fluid, E_correction=0):
    #close pipe and restrict outflow
    #input: (4xN) array fields with stored values of density, momentum, energy and pressure, (3xN+1) array fluxes to write fluxes to
    #output:(3xN+1) array with updated fluxes, passed by reference
    fluxes[0,-1]=0
    fluxes[1,-1]=eos.get_p(fields[0,-1], fields[1,-1]/fields[0,-1], fields[2,-1], fluid)
    fluxes[2,-1]=0
    
def get_outlet_bc_p_driven(fluxes, fields, p_atm, fluid, E_correction=0, last_state=[0]):
    p_outlet=p_atm
    #T_outlet=CP.PropsSI('T', 'P', fields[3,-1], 'D', fields[0,-1], fluid)
    #rho_outlet=CP.PropsSI('D', 'T', T_outlet, 'P', p_outlet, fluid)
    rho_outlet=fields[0,-1]
    u_outlet=fields[1,-1]/rho_outlet
    #u_outlet=0
    rho, u, p, E=get_middle_state(fields[1,-1]/fields[0,-1], fields[0,-1], fields[3,-1], u_outlet, rho_outlet, p_outlet, fluid)
    """
    print('left')
    print('rho', fields[0,-1], 'u', fields[1,-1]/fields[0,-1], 'p', fields[3,-1])
    print('right')
    print('rho', rho_outlet, 'u', u_outlet, 'p', p_outlet)
    print('middle')
    print('rho', rho, 'u', u, 'p', p)
    """
    fluxes[0,-1]=rho*u
    fluxes[1,-1]=rho*u*u+p
    fluxes[2,-1]=u*(E+p)
    if len(last_state)>=2:
        last_state[0]=rho
        last_state[1]=u
        last_state[2]=p
        last_state[3]=E
    
def energy_optimised_outlet_bc(fluxes, fields, E_correction):
    fluxes[0,-1]=fields[1,-1]
    fluxes[1,-1]=fields[1,-1]**2/fields[0,-1]+fields[3,-1]
    fluxes[2,-1]=max(0,fields[1,-1]/fields[0,-1]*(fields[2,-1]+fields[3,-1])-E_correction)
    
def calc_inlet_bc(rho, u, p, fluid):
    #calculate inlet fluxes given by BC. Fluxes are assumed to be constant, hence only one evaluation is required
    #input: scalars inlet density, velocity and pressure
    #output: (1x3) vector with inlet fluxes
    res=np.zeros([3])
    res[0]=rho*u
    res[1]=rho*u*u
    res[2]=u*(eos.get_E(rho, u, p, fluid)+p)
    return res

def get_inflow_bc(fluxes, inlet, fields):
    #add the inlet fluxes
    #input: (3xN+1) array to write inlet fluxes to, (1x3) vector containing fluxes
    #output: inlet fluxes added to fluxes in the first column, passed by reference
    fluxes[0,0]=inlet[0]
    fluxes[1,0]=inlet[1]+fields[3,0]
    fluxes[2,0]=inlet[2]
    
def get_inflow_bc_p_driven(fluxes, inlet, fields, fluid):
    #fix inlet pressure (only makes sense for pressure higher than ambient)
    #input: (3xN+1) array to write inlet fluxes to, (1x3) vector containing fluxes
    #output: inlet fluxes added to fluxes in the first column, passed by reference
    rho_l=CP.PropsSI('D', 'P', fields[3,0], 'T', inlet[3], fluid)
    rho, u, p, E=get_middle_state(inlet[1], rho_l, fields[3,0], fields[1,0]/fields[0,0], fields[0,0], fields[3,0], fluid)
    fluxes[0,0]=rho*u
    fluxes[1,0]=rho*u*u+p
    fluxes[2,0]=u*(E+p)

            
def update_fluxes(fluxes, fields, inlet_flux, inlet, p_atm, fluid, par, E_correction=0, last_state=[0]):
    #calculate fluxes at cell interfaces, assuming Riemann problems at cell boundaries
    #input: (4xN) array fields with stored values of density, momentum, energy and pressure, (3xN+1) array fluxes to write fluxes to
    #       vector of length 3 (inlet_flux) containing inlet fluxes
    #       vector of length 3 (inlet) containing (in this order) inlet density, velocity and pressure
    #output: (3xN+1) array with updated fluxes, passed by reference

    #inflow BC
    if par.rank==0:
        get_inflow_bc_p_driven(fluxes, inlet, fields, fluid)
    #exchange field info around border to neighbour cell
    border_upper, border_lower=para.exchange_border(fields, par)
                        
    for i in range(1, fluxes.shape[1]):
        if i>1 and i<fluxes.shape[1]-2:
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid)           
        elif i==1 and par.rank==0:
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, inlet=inlet, border=1)
            
        elif i==1:
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, inlet=inlet, border=2, border_val=border_lower)
            
        elif i==fluxes.shape[1]-2 and par.rank<par.num_procs-1:
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, inlet=inlet, border=-2, border_val=border_upper)
            
        elif i==fluxes.shape[1]-2:
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, border=-1, p_atm=p_atm)
            
        elif i==fluxes.shape[1]-1 and par.rank<par.num_procs-1:
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, inlet=inlet, border=-3, border_val=border_upper)
            
        elif i==fluxes.shape[1]-1:
            #outlet BC
            get_outlet_bc_p_driven(fluxes, fields, p_atm, fluid, E_correction, last_state)
            break
        
        else:
            print('error during flux calculation: not all cases were included')
            assert(False)
            
        fluxes[0,i]=rho_m*u_m
        fluxes[1,i]=rho_m*u_m*u_m+p_m
        fluxes[2,i]=u_m*(E_m+p_m)
    
    #comunicate border flux to neighbouring cell
    fluxes[:,0]=para.exchange_fluxes(fluxes,par)

    
    
def update_fluxes_upwind(fluxes, fields, inlet_flux, inlet, E_correction=0):
    #calculate flux at cell interface based on upwind scheme
    get_inflow_bc(fluxes, inlet_flux, fields)
    for i in range(1, fluxes.shape[1]-1):
        if fields[0,i]<0 or fields[2,i]<0 or fields[3,i]<0:
            assert(False)
        if fields[1,i-1]/fields[0,i-1]+fields[1,i]/fields[0,i]<0:
            fluxes[0,i]=fields[1,i]
            fluxes[1,i]=fields[1,i]**2/fields[0,i]+fields[3,i]
            fluxes[2,i]=fields[1,i]/fields[0,i]*(fields[2,i]+fields[3,i])
        else:
            fluxes[0,i]=fields[1,i-1]
            fluxes[1,i]=fields[1,i-1]**2/fields[0,i-1]+fields[3,i-1]
            fluxes[2,i]=fields[1,i-1]/fields[0,i-1]*(fields[2,i-1]+fields[3,i-1])
    zero_gradient_outlet_bc(fluxes, fields)