#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 09:58:26 2021

@author: pascal
"""
import numpy as np
import CoolProp.CoolProp as CP

def get_E(rho, u, p, fluid='Air'):
    #calculate energy based on ideal gas EOS
    #input: scalars density rho, velocity u, presure p, isentropic coefficient gamma
    #output: scalar energy E
    e=CP.PropsSI('Umass', 'P', p, 'D', rho, fluid)
    return rho*(e+0.5*u*u)

def get_e(rho, u, E):
    
    return (E/rho-0.5*u*u)

def get_p(rho, u, E, fluid='Air'):

    e=get_e(rho, u, E)    
    return CP.PropsSI('P', 'D', rho, 'Umass', e, fluid)

def get_cs(p, rho, gamma=1.4):
    #calculate speed of sound
    #input: scalars isentropic coefficient gamma, pressure p, density rho
    #output: scalar speed of sound
    if p/rho<0:
        print("negative speed of sound")
        assert(False)
    return np.sqrt(gamma*p/rho)

def get_chi(p,rho, fluid='Air', distance=0.01):
    #calculate speed of acoustic waves
    #input: scalars density rho and pressure p, string fluid (must be parsable for CoolProp), distance to adjust derrivative
    #output: scalar acoustic waves
    dedrho=(CP.PropsSI('Umass','P',p,'D',rho*(1+distance),fluid)-CP.PropsSI('Umass','P',p,'D',rho*(1-distance),fluid))/(2*rho*distance)
    dedp=(CP.PropsSI('Umass','P',p*(1+distance),'D',rho,fluid)-CP.PropsSI('Umass','P',p*(1-distance),'D',rho,fluid))/(2*p*distance)
    return np.sqrt((p/rho**2-dedrho)/dedp)