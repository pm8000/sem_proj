#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 09:58:26 2021

@author: pascal
"""
import numpy as np

def get_E(rho, u, p, gamma=1.4):
    #calculate energy based on ideal gas EOS
    #input: scalars density rho, velocity u, presure p, isentropic coefficient gamma
    #output: scalar energy E
    return p/(gamma-1)+0.5*rho*u*u

def get_p(rho, u, E, gamma=1.4):
    #calculate pressure based on ideal gas EOS
    #input: scalars or matrices with same dimensions density rho, velocity u, endergy E, isentropic coefficient gamma
    #output: scalar pressure p
    return (gamma - 1) * (E - 0.5*rho*(u*u))

def get_cs(p, rho, gamma=1.4):
    #calculate speed of sound
    #input: scalars isentropic coefficient gamma, pressure p, density rho
    #output: scalar speed of sound
    if p/rho<0:
        assert(False)
    return np.sqrt(gamma*p/rho)