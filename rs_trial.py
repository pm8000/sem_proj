#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 08:52:35 2021

@author: pascal
"""
import numpy as np
import math
import rs_trial_EOS as eos
import rs_trial_fluxes as flux
import rs_trial_source as src

L=0.1 #pipe length
dx=0.01 #grid density
N=math.ceil(L/dx) #number of cells
d=0.1 #pipe diameter

t_end=1e-3 #end time
cfl=0.5 #cfl number to define time step

alpha=200 #convection heat transfer coefficient
gamma=1.4
R=287.058
T_inlet=200+273
T_amb=273+20
p_amb=1e5
u_inlet=10
p_inlet=p_amb

rho_inlet=p_inlet/(R*T_inlet)
rho_amb=p_amb/(R*T_amb)

inlet=np.array([rho_inlet,u_inlet,p_inlet])
inlet_flux=flux.calc_inlet_bc(rho_inlet, u_inlet, p_inlet)

fields=np.zeros([4,N]) #0 stores density, 1 stores momentum, 2 stores energy, 3 stores pressure
fluxes=np.zeros([3, N+1]) #0 stores density flux, 1 stores momentum flux, 2 stores energy flux
T_wall=np.zeros(N)
T_wall[:]=20+273
T=np.zeros(N)
source=np.zeros([3,N])

#initialise fields

fields[0,:]=p_amb/(R*T_amb)
fields[1,:]=fields[0,:]*u_inlet
fields[2,:]=eos.get_E(rho_amb, 0, p_amb, gamma)
fields[3,:]=p_amb

#determine time step
cs_max=eos.get_cs(p_inlet, rho_inlet)
dt=cfl*dx/cs_max
nt=math.ceil(t_end/dt)

for i in range(50):
    T[:]=fields[3,:]/(fields[0,:]*R)
    print("step",i)
    if (i+1)*dt>t_end:
        dt=t_end-i*dt

    #compute fluxes
    #BC are added in the same step
    flux.update_fluxes(fluxes, fields, inlet_flux, inlet)
    #add source
    #src.add_source(source, alpha, d, fields, R, T_wall)
    #integrate with euler explicit
    fields[:-1,:]+=dt/dx*(fluxes[:,:-1]-fluxes[:,1:]+source[:,:])

