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
import rs_trial_animation as anim
import CoolProp.CoolProp as CP
import rs_trial_parallelisation as para
from mpi4py import MPI

plt.clf()
par=para.parallel()


L=1 #pipe length
dx=0.01 #grid density
N_tot=math.ceil(L/dx) #number of cells
N=para.get_N(N_tot, par)
d=0.1 #pipe diameter

t_end=0.01 #end time
cfl=0.5 #cfl number to define time step

fluid='Water'
nu=15e-6 #kinematic viscosity
gamma=1.4
R=287.058
T_inlet=120+273
T_amb=20+273
p_amb=1e5
u_inlet=1
p_inlet=1*p_amb
#E_loss=0
rho_inlet=CP.PropsSI('D', 'T',T_inlet, 'P',p_inlet,fluid)
T_boil=CP.PropsSI('T', 'P', p_amb, 'Q', 0.5,fluid)

inlet=np.array([rho_inlet,u_inlet,p_inlet,T_inlet])
inlet_flux=flux.calc_inlet_bc(rho_inlet, u_inlet, p_inlet, fluid)
#print(eos.get_E(rho_inlet, u_inlet, p_inlet))

fields=np.zeros([6, N]) #0 stores density, 1 stores momentum, 2 stores energy, 3 stores pressure, 4 stores temperature(, 5 stores velocity)
fluxes=np.zeros([3, N+1]) #0 stores density flux, 1 stores momentum flux, 2 stores energy flux
#T_wall=np.linspace(473,273,N)
T_wall=np.zeros(N)
T_wall[:]=T_amb
source=np.zeros([3,N])

#initialise fields

fields[0,:]=CP.PropsSI('D', 'P', p_amb, 'T', T_inlet,fluid)
fields[1,:]=fields[0,:]*u_inlet
fields[2,:]=eos.get_E(fields[0,:], fields[1,:]/fields[0,:], p_amb, fluid)
fields[3,:]=p_amb
fields[4,:]=T_inlet
fields[5,:]=fields[1,:]/fields[0,:]

#determine time step
upper_chi_est=eos.get_chi(p_inlet, rho_inlet, fluid) #chi is most likely largest at high temperature
dt=cfl*dx/upper_chi_est #maybe consider adaptive time step
nt=math.ceil(t_end/dt)
#E_correction=E_loss*dx/dt
first_point=np.zeros([6,nt])
last_state=[0,0,0,0]
p_evolution=np.zeros([11,first_point.shape[1]])
animation=np.zeros([int(t_end*4e4/3)+1,N,6])
animation_time=np.zeros(animation.shape[0])
#E_tot=np.sum(fields[2,:])
#E_tot_hist=np.zeros(nt+1)
#E_tot_hist[0]=E_tot
#dE_hist=np.zeros(E_tot_hist.size-1)
#E_flow=np.zeros(nt)
#t=np.linspace(0, nt*dt, nt)
#src_sum=np.zeros(22672)
if par.rank==0:
    start=time.time()
for i in range(nt):
    #if i%2.6e3==0:
    #   plt.plot(np.linspace(0.005,0.995,100),fields[3,:], label='t='+str(i*dt)+' s')
    if i%math.floor(nt/(animation.shape[0]-1))==0 and i/math.floor(nt/(animation.shape[0]-1))<animation.shape[0]:
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,0]=fields[0,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,1]=fields[1,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,2]=fields[2,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,3]=fields[3,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,4]=fields[4,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,5]=fields[5,:]
        animation_time[int(i/math.floor(nt/(animation.shape[0]-1)))]=i*dt
    if par.rank==0 and i%1e2==0:
        print("step",i)
    if (i+1)*dt>t_end:
        dt=t_end-i*dt
    #compute fluxes
    #BC are added in the same step
    #if i==15e3:
    #    flux.update_fluxes(fluxes, fields, inlet_flux, inlet, p_amb, fluid, E_correction, last_state)
    #else:
 
    flux.update_fluxes(fluxes, fields, inlet_flux, inlet, p_amb, fluid, par)
    #add source
    #src.add_momentum_source_2(source, fields, d, dx, nu)
    src.add_energy_source(source, T_boil, d, fields, R, T_wall, fluid)
    #src_sum[i]=np.sum(source[2,:])*dt
    #plt.plot(np.linspace(0,1,101),fluxes[0,:], label='rho*u')
    #plt.plot(np.linspace(0,1,101),fluxes[1,:], label='rho*u*u+p')
    #plt.plot(np.linspace(0.05,0.95,100),fields[3,:], label='pressure')
    #plt.legend()
    #plt.show()
    #integrate with euler explicit
    fields[:3,:]+=(dt/dx*(fluxes[:,:-1]-fluxes[:,1:])+dt*source[:,:])
    fields[3,:]=eos.get_p(fields[0,:], fields[1,:]/fields[0,:], fields[2,:], fluid)
    fields[4,:]=CP.PropsSI('T', 'P', fields[3,:], 'D', fields[0,:], fluid)
    fields[5,:]=fields[1,:]/fields[0,:]
    #E_tot=np.sum(fields[2,:])
    #out.write_energy_loss(i, dt, E_tot)
    #E_tot_hist[i+1]=E_tot
    #dE_hist[i]=E_tot_hist[i+1]-E_tot_hist[i]
    #E_flow[i]=(fluxes[2,0]-fluxes[2,-1])*dt/dx
    """
    for x in range(0,110,10):
        if x==100:
            x-=1
        p_evolution[math.ceil(x/10),i]=fields[5,x]
    """
    #first_point[0,i]=fields[0,0]
    #first_point[1,i]=fields[1,0]
    #first_point[2,i]=fields[2,0]
    #first_point[3,i]=fields[3,0]
    #first_point[4,i]=fields[4,0]
    #first_point[5,i]=fields[5,0]
    
if par.rank==0:
    print("time "+str(time.time()-start)+" sec")
#plt.plot(t,src_sum, label='src')
#plt.plot(t,dE_hist, label='dE')
#plt.plot(t, E_flow[:], label='net flux')
#plt.plot(t, E_flow[:]+src_sum[:], label='sum' )
#plt.plot(first_point[3,:])

#merge animation from different processes
animation=para.merge_results(animation, par)

if par.rank==0:
    #create graphs for output
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,0],'density_evolution_w_120_u_1_test.gif',0.5,0.68,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,1],'momentum_evolution_w_120_u_1_test.gif',0.0,3.25,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,2],'energy_evolution_w_120_u_1_test.gif',1.4e6,1.52e6,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,3],'pressure_evolution_w_120_u_1_test.gif',99600,100400,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,4],'temperature_evolution_w_120__u_1_test.gif',270,400,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,5],'velocity_evolution_w_120_u_1_test.gif',0,1.5,animation_time)
"""
plt.title('Velocity')
for i in range(p_evolution.shape[0]):
    plt.plot(np.linspace(0,t_end,nt),p_evolution[i,:], label='x= '+str(i*0.1)+'m')
plt.legend(bbox_to_anchor=(1.01,1),loc='upper left', borderaxespad=0)
plt.show()

max_pressure=[]
min_pressure=[]
for i in range(1,first_point.shape[1]-1):
    if first_point[3,i]>first_point[3,i-1] and first_point[3,i]>first_point[3,i+1]:
        max_pressure.append(first_point[3,i])
    if first_point[3,i]<first_point[3,i-1] and first_point[3,i]<first_point[3,i+1]:
        min_pressure.append(first_point[3,i])

plt.plot(max_pressure)
plt.savefig('result/max_pressure_w_200')
plt.clf()
plt.plot(min_pressure)
plt.savefig('result/min_pressure_w_200')

plt.clf()
plt.plot(first_point[0,:])
plt.title('density')
#plt.show()
plt.savefig('cell_20_density_fast')
plt.clf()
plt.plot(first_point[1,:], label='momentum')
plt.title('momentum')
plt.savefig('cell_20_momentum_fast')
#plt.legend()
#plt.show()
plt.clf()

print("density influx", inlet_flux[0])
print("momentum influx", inlet_flux[1])
print("energy influx", inlet_flux[2])

plt.plot(first_point[2,:])
plt.title('energy')
plt.savefig('cell_20_energy_fast')

plt.clf()
plt.plot(first_point[3,:])
plt.title('pressure')
#plt.show()
plt.savefig('cell_20_pressure_fast')

plt.clf()
plt.plot(first_point[4,:])
plt.title('velocity')
plt.savefig('cell_20_velocity_fast')
plt.clf()
plt.plot(first_point[5,:])
plt.title('temperature')
plt.savefig('cell_20_temperature_fast')
plt.clf()
"""
