#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 14:01:56 2021

@author: pascal
"""


def add_energy_source(res, alpha, d, fields, R, T_wall):
    #calculate source termsource
    #input:
    #output: (4xN) array with added  term, passed by reference
    res[2,:]= -alpha*4/d*(fields[3,:]/(fields[0,:]*R)-T_wall[:])
    
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