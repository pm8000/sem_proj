#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 14:01:56 2021

@author: pascal
"""


def add_source(res, alpha, d, fields, R, T_wall):
    #calculate source term
    #input:
    #output: (4xN) array with added source term, passed by reference
    res[2,:]= -alpha*4/d*(fields[3,:]/(fields[0,:]*R)-T_wall[:])