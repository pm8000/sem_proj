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

t_end=0.3  #end time
cfl=0.5 #cfl number to define time step

alpha=200 #convection heat transfer coefficient
nu=15e-6 #kinematic viscosity
gamma=1.4
R=287.058
T_inlet=200+273
T_amb=20+273
p_amb=1e5
u_inlet=10
p_inlet=1*p_amb
E_loss=0
rho_inlet=p_inlet/(R*T_inlet)
rho_amb=p_amb/(R*T_amb)

inlet=np.array([rho_inlet,u_inlet,p_inlet])
inlet_flux=flux.calc_inlet_bc(rho_inlet, u_inlet, p_inlet)
#print(eos.get_E(rho_inlet, u_inlet, p_inlet))

fields=np.zeros([5, N]) #0 stores density, 1 stores momentum, 2 stores energy, 3 stores pressure
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
first_point=np.zeros([6,26160])
#E_tot=np.sum(fields[2,:])
#E_tot_hist=np.zeros(nt+1)
#E_tot_hist[0]=E_tot
#dE_hist=np.zeros(E_tot_hist.size-1)
#E_flow=np.zeros(nt)
#t=np.linspace(0, nt*dt, nt)
#src_sum=np.zeros(22672)
start=time.time()
for i in range(nt):
    #if i%2.6e3==0:
    #    plt.plot(np.linspace(0.005,0.995,100),fields[3,:], label='t='+str(i*dt)+' s')
    # print(T_list)
    if i%2e3==0:
        print("step",i)
    if (i+1)*dt>t_end:
        dt=t_end-i*dt

    #compute fluxes
    #BC are added in the same step
    flux.update_fluxes(fluxes, fields, inlet_flux, inlet, E_correction)
    #add source
    src.add_momentum_source_2(source, fields, d, dx, nu)
    src.add_energy_source(source, alpha, d, fields, R, T_wall)
    #src_sum[i]=np.sum(source[2,:])*dt
    #plt.plot(np.linspace(0,1,101),fluxes[0,:], label='rho*u')
    #plt.plot(np.linspace(0,1,101),fluxes[1,:], label='rho*u*u+p')
    #plt.plot(np.linspace(0.05,0.95,100),fields[3,:], label='pressure')
    #plt.legend()
    #plt.show()
    #integrate with euler explicit
    fields[:3,:]+=(dt/dx*(fluxes[:,:-1]-fluxes[:,1:])+dt*source[:,:])
    fields[3,:]=eos.get_p(fields[0,:], fields[1,:]/fields[0,:], fields[2,:])
    T[:]=fields[3,:]/(fields[0,:]*R)
    fields[4,:]=fields[1,:]/fields[0,:]
    #E_tot=np.sum(fields[2,:])
    #out.write_energy_loss(i, dt, E_tot)
    #E_tot_hist[i+1]=E_tot
    #dE_hist[i]=E_tot_hist[i+1]-E_tot_hist[i]
    #E_flow[i]=(fluxes[2,0]-fluxes[2,-1])*dt/dx
    
    first_point[0,i]=fields[0,20]
    first_point[1,i]=fields[1,20]
    first_point[2,i]=fields[2,20]
    first_point[3,i]=fields[3,20]
    first_point[4,i]=fields[4,20]
    first_point[5,i]=T[20]
    
print("time "+str(time.time()-start)+" sec")
#plt.plot(t,src_sum, label='src')
#plt.plot(t,dE_hist, label='dE')
#plt.plot(t, E_flow[:], label='net flux')
#plt.plot(t, E_flow[:]+src_sum[:], label='sum' )
#plt.plot(first_point[3,:])
plt.title('Pressure')
plt.legend(bbox_to_anchor=(1.01,1),loc='upper left', borderaxespad=0)
plt.show()
#plt.clf()
"""
plt.plot(first_point[0,:])
plt.title('density')
#plt.show()
plt.savefig('first_cell_density_slow_2')
plt.clf()
plt.plot(first_point[1,:], label='momentum')
plt.title('momentum')
plt.savefig('first_cell_momentum_slow_2')
#plt.legend()
#plt.show()
plt.clf()

print("density influx", inlet_flux[0])
print("momentum influx", inlet_flux[1])
print("energy influx", inlet_flux[2])

plt.plot(first_point[2,:])
plt.title('energy')
plt.savefig('first_cell_energy_slow_2')
"""
plt.clf()
plt.plot(first_point[3,:])
plt.title('pressure')
plt.show()
#plt.savefig('first_cell_pressure_slow_2')
"""
plt.clf()
plt.plot(first_point[4,:])
plt.title('velocity')
plt.savefig('first_cell_velocity_slow_2')
plt.clf()
plt.plot(first_point[5,:])
plt.title('temperature')
plt.savefig('first_cell_temperature_slow_2')
plt.clf()
"""
