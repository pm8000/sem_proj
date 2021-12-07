#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 08:52:35 2021

@author: pascal
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import rs_trial_EOS as eos
import rs_trial_fluxes as flux
import rs_trial_source as src
import rs_trial_output as out
plt.clf()

L=1 #pipe length
dx=0.01 #grid density
N=math.ceil(L/dx) #number of cells
d=0.1 #pipe diameter

t_end=5e-1 #end time
cfl=0.5 #cfl number to define time step

alpha=200 #convection heat transfer coefficient
gamma=1.4
R=287.058
T_inlet=200+273
T_amb=273+20
p_amb=1e5
u_inlet=10
p_inlet=p_amb
E_loss=10

rho_inlet=p_inlet/(R*T_inlet)
rho_amb=p_amb/(R*T_amb)

inlet=np.array([rho_inlet,u_inlet,p_inlet])
inlet_flux=flux.calc_inlet_bc(rho_inlet, u_inlet, p_inlet)
print(eos.get_E(rho_inlet, u_inlet, p_inlet))

fields=np.zeros([5,N]) #0 stores density, 1 stores momentum, 2 stores energy, 3 stores pressure
fluxes=np.zeros([3, N+1]) #0 stores density flux, 1 stores momentum flux, 2 stores energy flux
T_wall=np.zeros(N)
T_wall[:]=20+273
T=np.zeros(N)
T[:]=T_wall[:]
source=np.zeros([3,N])

#initialise fields

fields[0,:]=p_amb/(R*T[:])
fields[1,:]=fields[0,:]*u_inlet
fields[2,:]=eos.get_E(rho_amb, fields[1,:]/fields[0,:], p_amb, gamma)
fields[3,:]=p_amb
fields[4,:]=fields[1,:]/fields[0,:]

#determine time step
cs_max=eos.get_cs(p_inlet, rho_inlet)
dt=cfl*dx/cs_max
nt=math.ceil(t_end/dt)
E_correction=E_loss*dx/dt
first_point=np.zeros([5,1])
E_tot=np.sum(fields[2,:])
E_tot_hist=np.zeros(nt+1)
E_tot_hist[0]=E_tot
t=np.linspace(0, nt*dt, nt+1)

start=time.time()
for i in range(nt):
    if i%1e4==0:
        plt.plot(np.linspace(0.005,0.995,100),fields[3,:], label='t='+str(i*dt)+' s')
        1+1
    #print(T_list)
    if i%2000==0:
        print("step",i)
    if (i+1)*dt>t_end:
        dt=t_end-i*dt

    #compute fluxes
    #BC are added in the same step
    flux.update_fluxes(fluxes, fields, inlet_flux, inlet, E_correction)
    #add source
    src.add_source(source, alpha, d, fields, R, T_wall)
    plt.plot(np.linspace(0,1,101),fluxes[0,:], label='rho*u')
    plt.plot(np.linspace(0,1,101),fluxes[1,:], label='rho*u*u+p')
    plt.plot(np.linspace(0,1,101),fluxes[2,:], label='u*(E+p)')
    plt.legend()
    plt.show()
    #integrate with euler explicit
    fields[:3,:]+=(dt/dx*(fluxes[:,:-1]-fluxes[:,1:])+dt*source[:,:])
    fields[3,:]=eos.get_p(fields[0,:], fields[1,:]/fields[0,:], fields[2,:])
    T[:]=fields[3,:]/(fields[0,:]*R)
    fields[4,:]=fields[1,:]/fields[0,:]
    #E_tot=np.sum(fields[2,:])
    #out.write_energy_loss(i, dt, fields[2,-1])
    #E_tot_hist[i+1]=E_tot
    """
    first_point[0,i]=fields[0,0]
    first_point[1,i]=fields[1,0]
    first_point[2,i]=fields[2,0]
    first_point[3,i]=fields[3,0]
    first_point[4,i]=T[0]
    """
print("time "+str(time.time()-start))
#plt.plot(t,E_tot_hist)
plt.legend()
plt.show()
"""
plt.plot(first_point[0,:], label='source')
plt.title('density')
plt.savefig('first_cell_density_supers')
plt.clf()
plt.plot(first_point[1,:], label='flux difference')
plt.title('velocity')
plt.savefig('first_cell_velocity_supers')
#plt.legend()
#plt.show()
plt.clf()

print("density influx", inlet_flux[0])
print("momentum influx", inlet_flux[1])
print("energy influx", inlet_flux[2])

plt.plot(first_point[2,:])
plt.title('energy flux')
plt.savefig('first_cell_energy_supers')
plt.clf()
plt.plot(first_point[3,:])
plt.title('pressure')
plt.savefig('first_cell_pressure_supers')
plt.clf()
plt.plot(first_point[4,:])
plt.title('temperature')
plt.savefig('first_cell_temperature_supers')
plt.clf()
"""
