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
    #inputs: floats u (function value in current cell), u_p (function value in neigbouring cell to the right)
    #        float u_m (function value in neighbouring cell to the left)
    #output: float slope inside cell
    d1=u_p-u #slope to the right
    d2=u-u_m #slope to the left
    if np.sign(d1)!=np.sign(d2): #local extreme point
        return 0
    return np.sign(d1)*min(abs(d1),abs(d2))

def get_middle_state(u_l, rho_l, p_l, u_r, rho_r, p_r, fluid=None,state='real', ideal_gas=None):
    #calculate middle state based on characteristics of left and right state:
    #inputs: floats u_l (left velocity), rho_l (left density), p_l (left pressure)
    #        floats u_r (right velocity), rho_r (right density), p_r (right pressure)
    #        string state (real or ideal depending on equation of state used), optional (default is real)
    #        string fluid (only required for real gases, must be known to CoolProp)
    #        ideal_gas object ideal_gas (only required for ideal gas)
    #output: tuple of floats for values in middle state
    #        floats rho_m (density), u_m(velocity), p_m(pressure), E_m(total energy per volume)
    if state=='real':
        chi_l=eos.get_chi(p_l, rho_l, fluid)
        chi_r=eos.get_chi(p_r, rho_r, fluid)
    elif state=='ideal':
        chi_l=eos.get_cs(p_l, rho_l, ideal_gas=ideal_gas)
        chi_r=eos.get_cs(p_r, rho_r, ideal_gas=ideal_gas)
    
    if u_l>chi_l:
        print("supersonic velocity encountered")
        u_m=u_l
        p_m=p_l
        rho_m=rho_l
    else: #subsonic case
        u_m=(rho_l*u_l*chi_l+rho_r*u_r*chi_r+p_l-p_r)/(rho_l*chi_l+rho_r*chi_r)
        p_m=p_l-rho_l*chi_l*(u_m-u_l)
        if u_l+u_r>=0: #positive velocity
            rho_m=(p_m-p_l)/(chi_l*chi_l)+rho_l
        else: #negative velocity
            rho_m=(p_m-p_r)/(chi_r*chi_r)+rho_r
    E_m=eos.get_E(rho_m, u_m, p_m, fluid, state=state, ideal_gas=ideal_gas)
    return (rho_m, u_m, p_m, E_m)

def get_border_values(fields, i, fluid=None, inlet=None, border=0, p_out=None, border_val=None, state='real', ideal_gas=None):
    #calculate the values at the cell interface assuming linear profile within cell limited by minmod
    #inputs: array fields(6xN, all field variables at all locations in the domain)
    #        int i (current cell, whose interface to the left to be evaluated)
    #        int border (if cell is close to (whole or thread domain) border, additional information is required:
    #                   set border=1 for interface between first and second cell of the whole domain (use BC)
    #                   set border=2 for interface between first and second cell of the thread domain (use values passed from neighbouring thread)
    #                   set border=-1 for interface between second last and last cell of whole domain (use BC)
    #                   set border=-2 for interface between second last and last cell of thread domain (use values passed from neighbouring thread)
    #                   set border =-3 for interface between last cell of current thread and first cell of next thread (use values passed from neighbouring thread)
    #        string state (real or ideal depending of equation of state used), optional (default is real)
    #        string fluid (only required for real fluid, must beknown to CoolProp)
    #        ideal_gas object ideal_gas (only required for real flow)
    #        array border_val (3 for positive border, 3x2 for negative border, values in the adjacent cells in the next/last threads), only required if abs(border)>1
    #        float p_out (outlet pressure, required if border=1)
    #        array inlet (at least 3 , element 0 denotes inlet density, 1 denotes inlet velocity, 2 denotes inlet pressure)       
    #output: tuple of floats for values in middle state
    #        floats rho_m (density), u_m(velocity), p_m(pressure), E_m(total energy per volume)
    if (border <-3 or border >2):
        print("encountered invalid value for border in function get_middle_state")
        assert(False)
    #get values at cell boundary
    #boundary cases
    if border==1:
        rho_l=fields[0,i-1]+0.5*minmod(fields[0, i-1], fields[0,i], 2*inlet[0]-fields[0,i-1])
        u_l=fields[1,i-1]/fields[0, i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], 2*inlet[1]-fields[1,i-1]/fields[0,i-1])
        #p_l=fields[3,i-1]
        p_l=fields[3,i-1]+0.5*minmod(fields[3, i-1], fields[3,i], 2*(inlet[2])-fields[3,i-1])
    
    elif border==2:
        rho_l=fields[0, i-1]+0.5*minmod(fields[0,i-1], fields[0,i], border_val[0])
        u_l=fields[1,i-1]/fields[0,i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], border_val[1])
        p_l=fields[3,i-1]+0.5*minmod(fields[3,i-1], fields[3,i], border_val[2])
        
    else:
        rho_l=fields[0, i-1]+0.5*minmod(fields[0,i-1], fields[0,i], fields[0,i-2])
        u_l=fields[1,i-1]/fields[0,i-1]+0.5*minmod(fields[1,i-1]/fields[0,i-1], fields[1,i]/fields[0,i], fields[1,i-2]/fields[0,i-2])
        p_l=fields[3,i-1]+0.5*minmod(fields[3,i-1], fields[3,i], fields[3,i-2])
    if rho_l<0 or pd.isna(rho_l)==True:
        print("density is negative or not a number")
        assert(False)
    
    if border==-1:
        rho_r=fields[0,i]
        u_r=fields[1,i]/rho_r
        p_r=fields[3,i]-0.5*minmod(fields[3, i-1], fields[3,i], 2*(p_out)-fields[3,i-1])
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
        
    return(get_middle_state(u_l, rho_l, p_l, u_r, rho_r, p_r, fluid, state=state, ideal_gas=ideal_gas))

def zero_gradient_outlet_bc(fluxes, fields):
    #calculate outflux based on zero gradient assumption (like incompressible flow)
    #currently not in use
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array fields (6xN, all field variables at all locations in domain)
    #output: void, the outflux is added to the last column of fluxes
    
    fluxes[0,-1]=fields[1,-1]
    fluxes[1,-1]=fields[1,-1]**2/fields[0,-1]+fields[3,-1]
    fluxes[2,-1]=fields[1,-1]/fields[0,-1]*(fields[2,-1]+fields[3,-1])
    
def no_outlet(fluxes, fields, fluid=None, state='real', ideal_gas=None):
    #outlet boundary condition if pipe outlet is closed (no flow)
    #currently not in use
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array fields (6xN, all field variables at all locations in domain)
    #        string state (real or ideal depending of equation of state used), optional (default is real)
    #        string fluid (only required for real fluid, must beknown to CoolProp)
    #        ideal_gas object ideal_gas (only required for real flow)
    #output: void, the outflux is added to the last column of fluxes
    fluxes[0,-1]=0
    fluxes[1,-1]=eos.get_p(fields[0,-1], fields[1,-1]/fields[0,-1], fields[2,-1], fluid, state=state, ideal_gas=ideal_gas)
    fluxes[2,-1]=0
    
def get_outlet_bc_p_driven(fluxes, fields, p_out, fluid=None, state='real', ideal_gas=None):
    #outlet boundary condotion with ghost cell, where values are specified and flux is computed through chracteristics at interface
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array fields (6xN, all field variables at all locations in domain)
    #        float p_out (outlet pressure)
    #        string state (real or ideal depending of equation of state used), optional (default is real)
    #        string fluid (only required for real fluid, must beknown to CoolProp)
    #        ideal_gas object ideal_gas (only required for real flow)
    #output: void, the outflux is added to the last column of fluxes

    p_outlet=p_out
    #T_outlet=CP.PropsSI('T', 'P', fields[3,-1], 'D', fields[0,-1], fluid)
    #rho_outlet=CP.PropsSI('D', 'T', T_outlet, 'P', p_outlet, fluid) #calculate outlet density based on outlet pressure, difference is negigeable
    rho_outlet=fields[0,-1]
    u_outlet=fields[1,-1]/rho_outlet
    #u_outlet=0 #alternative boundary condition (reservoir at outflow)
    rho, u, p, E=get_middle_state(fields[1,-1]/fields[0,-1], fields[0,-1], fields[3,-1], u_outlet, rho_outlet, p_outlet, fluid, state=state, ideal_gas=ideal_gas)
    fluxes[0,-1]=rho*u
    fluxes[1,-1]=rho*u*u+p
    fluxes[2,-1]=u*(E+p)
    
def energy_optimised_outlet_bc(fluxes, fields, E_correction):
    #outlet boundary condition with reduced energy outflow or if resulting energy outflow is negative, no energy outflow
    #energy is usually reduced to force constant energy within the domain (E_in-E_out-Q_source=0)
    #not in use
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array fields (6xN, all field variables at all locations in domain)
    #        float E_correction (amount of energy to reduce)
    #output: void, the outflux is added to the last column of fluxes
    
    fluxes[0,-1]=fields[1,-1]
    fluxes[1,-1]=fields[1,-1]**2/fields[0,-1]+fields[3,-1]
    fluxes[2,-1]=max(0,fields[1,-1]/fields[0,-1]*(fields[2,-1]+fields[3,-1])-E_correction)
    
def calc_inlet_bc(rho, u, p, fluid=None, state='real', ideal_gas=None):
    #calculate inlet flux, for directly specified inlet flux (no use of ghost cells, like incompressible flow)
    #pressure in momentum flux is left out to calculate later based on zero gradient
    #can be calculated once and added at every time step, because this influx is constant in time
    #not in use
    #inputs: floats rho (density), u (velocity), p (pressure)
    #        string state (real or ideal depending of equation of state used), optional (default is real)
    #        string fluid (only required for real fluid, must beknown to CoolProp)
    #        ideal_gas object ideal_gas (only required for real flow)
    #output: array res(3, contains inlet flux, pressure in momentum flux is left out, will be added later based on zero gradient method)
    res=np.zeros([3])
    res[0]=rho*u
    res[1]=rho*u*u
    res[2]=u*(eos.get_E(rho, u, p, fluid, state=state, ideal_gas=ideal_gas)+p)
    return res

def get_inflow_bc(fluxes, inlet, fields):
    #add inlet flux calculated by calc_inlet_bc
    #not in use
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array inlet (3, vector with fluxes, momentum flux must not contain pressure, this is added by zero gradient bc)
    #        array fields (6xN, all field variables at all locations in domain)
    #output: void, the influx is added to the first column of fluxes
    fluxes[0,0]=inlet[0]
    fluxes[1,0]=inlet[1]+fields[3,0]
    fluxes[2,0]=inlet[2]
    
def get_inflow_bc_p_driven(fluxes, inlet, fields, fluid=None, ideal_gas=None, state='real'):
    #inlet boundary condotion with ghost cell, where values are specified and flux is computed through chracteristics at interface
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array inlet (at least 3 , element 0 denotes inlet density, 1 denotes inlet velocity, 2 denotes inlet pressure)
    #        array fields (6xN, all field variables at all locations in domain)
    #        string state (real or ideal depending of equation of state used), optional (default is real)
    #        string fluid (only required for real fluid, must beknown to CoolProp)
    #        ideal_gas object ideal_gas (only required for real flow)
    #output: void, the influx is added to the first column of fluxes
    
    if state=='real':
        rho_l=CP.PropsSI('D', 'P', inlet[2], 'T', inlet[3], fluid) #use this line for fix inlet pressure
        #rho_l=CP.PropsSI('D', 'P', fields[3,0], 'T', inlet[3], fluid) #use this line for zero gradient inlet pressure
    else:
        rho_l=inlet[2]/(ideal_gas().R*inlet[3]) #use this line for fix inlet pressure
        #rho_l=fields[3,0]/(ideal_gas.R*inlet[3]) #use this line for zero gradient inlet pressure
        
    rho, u, p, E=get_middle_state(inlet[1], rho_l, inlet[2], fields[1,0]/fields[0,0], fields[0,0], fields[3,0], fluid, state=state, ideal_gas=ideal_gas) #use this line for fix inlet pressure
    #rho, u, p, E=get_middle_state(inlet[1], rho_l, fields[3,0], fields[1,0]/fields[0,0], fields[0,0], fields[3,0], fluid) #use this line for zero gradient inlet pressure
    fluxes[0,0]=rho*u
    fluxes[1,0]=rho*u*u+p
    fluxes[2,0]=u*(E+p)

            
def update_fluxes(fluxes, fields, inlet_flux, inlet, p_out, par, fluid=None, ideal_gas=None, E_correction=0, state='real'):
    #calculate fluxes at cell interfaces, based on characteristics (Godunov method)
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array fields (6xN, all field variables at all locations in domain)
    #        array inlet_flux (3, vector with fluxes for fixed influx bc, momentum flux must not contain pressure, this is added by zero gradient bc)
    #        array inlet (at least 3 , element 0 denotes inlet density, 1 denotes inlet velocity, 2 denotes inlet pressure)
    #        float p_out (outlet pressure)
    #        parallel object par
    #        string state (real or ideal depending of equation of state used), optional (default is real)
    #        string fluid (only required for real fluid, must beknown to CoolProp)
    #        ideal_gas object ideal_gas (only required for real flow)
    #        float E_correction (amount of energy to reduce if energy_optimised_outlet_bc is used), optional (default 0)
    #output: void, fluxes is overwtritten with fluxes at current time step

    #inflow BC
    if par.rank==0:
        get_inflow_bc_p_driven(fluxes, inlet, fields, fluid, state=state, ideal_gas=ideal_gas)
        
    #exchange field info around border to neighbour cell
    border_upper, border_lower=para.exchange_border(fields, par)

    #loop over field array to compute fluxes                   
    for i in range(1, fluxes.shape[1]):
        # get the middle state values
        if i>1 and i<fluxes.shape[1]-2: #regular case, no boundaries involved
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, state=state, ideal_gas=ideal_gas)  
            
        elif i==1 and par.rank==0: #interface between first and second cell of domain
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, inlet=inlet, border=1, state=state, ideal_gas=ideal_gas)
            
        elif i==1: #interface between first and second cell of a thread
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, border=2, border_val=border_lower, state=state, ideal_gas=ideal_gas)
            
        elif i==fluxes.shape[1]-2 and par.rank<par.num_procs-1: #interface between second last and last cell of a thread
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, border=-2, border_val=border_upper, state=state, ideal_gas=ideal_gas)
            
        elif i==fluxes.shape[1]-2: #interface between second last and last cell of domain
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, border=-1, p_out=p_out, state=state, ideal_gas=ideal_gas)
            
        elif i==fluxes.shape[1]-1 and par.rank<par.num_procs-1: #interface between last cell of thread and first of next thread
            rho_m, u_m, p_m, E_m=get_border_values(fields, i, fluid, inlet=inlet, border=-3, border_val=border_upper, state=state, ideal_gas=ideal_gas)
            
        elif i==fluxes.shape[1]-1: #outlet BC (fluxes are added directly, therefore loop is broken out)
            get_outlet_bc_p_driven(fluxes, fields, p_out, fluid, state=state, ideal_gas=ideal_gas)
            break
        
        else:
            print('error during flux calculation: not all cases were included')
            assert(False)
        #caluclate and add fluxes based on middle state    
        fluxes[0,i]=rho_m*u_m
        fluxes[1,i]=rho_m*u_m*u_m+p_m
        fluxes[2,i]=u_m*(E_m+p_m)
    #comunicate border flux to neighbouring cell
    fluxes[:,0]=para.exchange_fluxes(fluxes,par)

    
    
def update_fluxes_upwind(fluxes, fields, inlet_flux, inlet):
    #calculate flux at cell interface based on upwind scheme with incompressible boundary conditions (inlet flux, and zero gradient outflux)
    #not in use
    #inputs: array fluxes (3xN+1, fluxes between cells and in/outflux)
    #        array fields (6xN, all field variables at all locations in domain)
    #        array inlet_flux (3, vector with fluxes for fixed influx bc, momentum flux must not contain pressure, this is added by zero gradient bc)
    #        array inlet (at least 3 , element 0 denotes inlet density, 1 denotes inlet velocity, 2 denotes inlet pressure)
    #output: void, fluxes is overwtritten with fluxes at current time step
    
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