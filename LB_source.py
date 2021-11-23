#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:15:13 2021

@author: pascal
"""

import numpy as np

def get_alpha(U, D, nu, Pr, l_f):
    #calculate alpha by the Dittus-Boelter correletion
    #inputs: scalars inlet velocity U, pipe diameter D, viscosity nu, fluid thermal conductvity l_f
    #output: scalar heat transfer coefficient
    return 0.0023*(U*D/nu)**0.8*Pr**0.4*l_f/D

def get_source(T, T_wall, alpha, d):
    #calculate current macroscopic source term
    #input: (1xN) temperature fields of stream and wall, scalars pipe diameter d and heat transfer coefficient alpha
    #output: (1xN) heat transfer field
    return alpha*4/d*(T[:]-T_wall[:])

def get_dqdt(source, source_last, dt):
    #calculate time derrivative of heat source
    #input: (1xN) source of current and last time step, scalar time step dt
    #output: (1xN) time derrivative
    return (source[:]-source_last[:])/dt

def get_source_populations(source, source_last, T, T_wall, alpha, d, dt, order=2):
    #calculate source term for the individual populations
    #inputs: result matrix source (3xN), temperature fields (1xN) of stream an wall, last source term required for time derrivative for second order,
    #            scalars heat transfer coefficient alpha, pipe diameter d and time step dt, parameter order (1 or 2) to decide on order of accuracy
    #output:source term for populations based on first order approximation, result is passed by reference
    if (order!=1) and (order!=2):
        print("this order for the heat sink is not supported")
        assert(False)
    source_total=get_source(T, T_wall, alpha, d)
    source_correction=np.zeros(T.size)
    if order==2:
        source_correction=dt*get_dqdt(source_total, source_last, dt)
        source_last=source_total
    source[0,:]=-1.0/6*dt*(source_total[:]+source_correction[:])
    source[1,:]=-2.0/3*dt*(source_total[:]+source_correction[:])
    source[2,:]=source[0,:]