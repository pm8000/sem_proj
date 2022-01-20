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
    def __init__(self,gamma,R):
        self.gamma=gamma
        self.R=R
        self.c_v=self.R/(self.gamma-1)
        self.c_p=self.c_v*self.gamma

def get_E(rho, u, p, fluid=None, state='real', ideal_gas=None):
    #calculate energy based on ideal gas EOS
    #input: scalars density rho, velocity u, presure p, isentropic coefficient gamma
    #output: scalar energy E
    if state=='real':
        e=CP.PropsSI('Umass', 'P', p, 'D', rho, fluid)
        return rho*(e+0.5*u*u)
    elif state=='ideal':
        return p/(ideal_gas.gamma-1)+0.5*rho*u*u
    else:
        print('unknown equation of state')
        assert(False)

def get_e(rho, u, E):
    
    return (E/rho-0.5*u*u)

def get_p(rho, u, E, fluid=None, state='real', ideal_gas=None):

    if state=='real':    
        e=get_e(rho, u, E)
        return CP.PropsSI('P', 'D', rho, 'Umass', e, fluid)
    elif state=='ideal':
        return (E-0.5*rho*u*u)*(ideal_gas.gamma-1)
    else:
        print('unknown equation of state')
        assert (False)

def get_cs(p, rho, ideal_gas=None):
    #calculate speed of sound
    #input: scalars isentropic coefficient gamma, pressure p, density rho
    #output: scalar speed of sound
    if p/rho<0:
        print("negative speed of sound")
        assert(False)
    return np.sqrt(ideal_gas.gamma*p/rho)

def get_chi(p,rho, fluid, distance=0.001, exact=False):
    #calculate speed of acoustic waves
    #input: scalars density rho and pressure p, string fluid (must be parsable for CoolProp), distance to adjust derrivative
    #set bool exact true to use CoolProp, otherwise approximations are used
    #output: scalar acoustic waves
    if rho==0:
        #no mediusm for soundwaves to propagate
        return 0
    if exact==False and fluid=='Water' and p>=99500 and p<=100500:
        dedrho=get_dedrho(rho, fluid)
        dedp=get_dedp(p, rho, fluid)
    else:
        dedrho=(CP.PropsSI('Umass','P',p,'D',rho*(1+distance),fluid)-CP.PropsSI('Umass','P',p,'D',rho*(1-distance),fluid))/(2*rho*distance)
        dedp=(CP.PropsSI('Umass','P',p*(1+distance),'D',rho,fluid)-CP.PropsSI('Umass','P',p*(1-distance),'D',rho,fluid))/(2*p*distance)
    return np.sqrt((p/rho**2-dedrho)/dedp)