#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 09:58:26 2021

@author: pascal
"""
import numpy as np
import CoolProp.CoolProp as CP
from fits.rs_trial_dedrho_table import get_dedrho
from fits.rs_trial_dedp_table import get_dedp

class ideal_gas:
    #class contains fluid properties for ideal gas
    #simplifies argument transfer to functions
    def __init__(self,gamma,R):
        #inputs: float gamma (specific heat ratio)
        #        float R (specific gas constant)
        self.gamma=gamma
        self.R=R
        self.c_v=self.R/(self.gamma-1) #specific heat constant at constant volume
        self.c_p=self.c_v*self.gamma #specific heat constant at constant pressure

def get_E(rho, u, p, fluid=None, state='real', ideal_gas=None):
    #calculate total energy per volume based on equation of state
    #input: floats rho (density), u (velocity), p (pressure)
    #       string state (real or ideal, choose eq. of state), optional (default is real)
    #       string fluid (must be known to CoolProp), only required for real gas e.o.s.
    #       ideal_gas class ideal_gas, only required for ideal gas e.o.s.
    #output: float E (total energy per volume)
    if state=='real':
        e=CP.PropsSI('Umass', 'P', p, 'D', rho, fluid)
        return rho*(e+0.5*u*u)
    elif state=='ideal':
        return p/(ideal_gas.gamma-1)+0.5*rho*u*u
    else:
        print('unknown equation of state')
        assert(False)

def get_e(rho, u, E):
    #calculate mass specific internal energy, based on total energy per volume
    #inputs: floats rho (density), u (velocity), E (total energy per volume)
    #output: float e (mass specific internal energy)
    return (E/rho-0.5*u*u)

def get_p(rho, u, E, fluid=None, state='real', ideal_gas=None):
    #calculate pressure based on equation of state
    #inputs: floats rho (density), u (velocity), E (total energy per volume)
    #        string state (real or ideal, choose eq. of state), optional (default is real)
    #        string fluid (must be known to CoolProp), only required for real gas e.o.s.
    #        ideal_gas class ideal_gas, only required for ideal gas e.o.s.
    #output float pressure
    if state=='real':    
        e=get_e(rho, u, E)
        return CP.PropsSI('P', 'D', rho, 'Umass', e, fluid)
    elif state=='ideal':
        return (E-0.5*rho*u*u)*(ideal_gas.gamma-1)
    else:
        print('unknown equation of state')
        assert (False)

def get_cs(p, rho, ideal_gas=None):
    #calculate speed of sound for an ideal gas
    #input: floats p (pressure), rho (density)
    #       ideal_gas class ideal_gas
    #output: float speed of sound
    if p/rho<0:
        print("negative speed of sound")
        assert(False)
    return np.sqrt(ideal_gas.gamma*p/rho)

def get_chi(p,rho, fluid, distance=0.001, exact=False):
    #calculate speed of acoustic waves for real gases
    #input: floats p (pressure), rho (density)
    #       string fluid (must be known to CoolProp)
    #       float distance (spacing of the two points to evaluate derivative), optional
    #       bool exact (False uses approximations if possible, True always uses CoolProp), optional (default False)
    #output: float acourstic wave speed
    if rho==0:
        #no mediusm for soundwaves to propagate
        return 0
    if exact==False and fluid=='Water' and p>=99000 and p<=101000:
        dedrho=get_dedrho(rho, fluid)
        dedp=get_dedp(p, rho, fluid)
    else:
        dedrho=(CP.PropsSI('Umass','P',p,'D',rho*(1+distance),fluid)-CP.PropsSI('Umass','P',p,'D',rho*(1-distance),fluid))/(2*rho*distance)
        dedp=(CP.PropsSI('Umass','P',p*(1+distance),'D',rho,fluid)-CP.PropsSI('Umass','P',p*(1-distance),'D',rho,fluid))/(2*p*distance)
    return np.sqrt((p/rho**2-dedrho)/dedp)