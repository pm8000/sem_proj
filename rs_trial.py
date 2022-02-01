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
import pandas as pd
import rs_trial_EOS as eos
import rs_trial_fluxes as flux
import rs_trial_source as src
import rs_trial_output as out
import rs_trial_animation as anim
import CoolProp.CoolProp as CP
import rs_trial_parallelisation as para


plt.clf()
par=para.parallel()


L=1 #pipe length
dx=0.009 #grid density
N_tot=math.ceil(L/dx) #number of cells
N=para.get_N(N_tot, par)
d=0.1 #pipe diameter

t_end=1.5 #end time
cfl=0.4 #cfl number to define time step

state='real'  #real or ideal
start_region=2 #type 1 if initial condition is single phase, or 2 if initial condition is two phase, only considered for real gases
divisor=1
assert(state=='real' or state=='ideal')
assert(start_region==1 or start_region==2)
fluid='Water' #fluid must be known for CoolProp (real gas only)
gamma=1.4   #ideal gas only
R=287.058   #ideal gas only
if state=='ideal':
    ideal_gas=eos.ideal_gas(gamma, R)
else:
    ideal_gas=None
#nu=15e-6 #kinematic viscosity, required only if wall shear stress is included
T_inlet=120+273 #inlet is always assumed to be in a single phase (gas) region
T_amb=20+273 #initial temperature, considered if start_region==1
q_init=0.2 #initial steam quality, considered if start_region==2
p_amb=1e5
u_inlet=3
p_inlet=1*p_amb
#E_loss=0
if state=='real':
    rho_inlet=CP.PropsSI('D', 'T',T_inlet, 'P',p_inlet,fluid)
    T_boil=CP.PropsSI('T', 'P', p_amb, 'Q', 0.5,fluid) #find boiling temperature (Q value is arbitrary number between (excluding) 0 and 1)
    if start_region==2:
        T_amb=T_boil
else:
    rho_inlet=p_inlet/(ideal_gas.R*T_inlet)
    T_boil=None
appendix='120_to_q_0.2_w_cool_50_u_3_2500'

inlet=np.array([rho_inlet,u_inlet,p_inlet,T_inlet])
inlet_flux=flux.calc_inlet_bc(rho_inlet, u_inlet, p_inlet, fluid,state=state, ideal_gas=ideal_gas)
#print(eos.get_E(rho_inlet, u_inlet, p_inlet))

fields=np.zeros([6, N]) #0 stores density, 1 stores momentum, 2 stores energy, 3 stores pressure, 4 stores temperature(, 5 stores velocity)
fluxes=np.zeros([3, N+1]) #0 stores density flux, 1 stores momentum flux, 2 stores energy flux

T_wall=np.zeros(N)
if start_region==1:
    T_wall[:]=T_amb
elif start_region==2:
    T_wall[:]=T_boil-50
"""
#uncomment and modify if the wall temperature is not uniform
#CAVEAT if code is executed on multiple cores the T_wall array is only the domain of one process
#CAVEAT use if statements and par.rank to assign different profiles to different processes
if par.rank==0:
    T_wall=np.zeros(N)
    T_wall[:]=T_inlet
elif par.rank==1:
    T_wall=np.linspace(473,293,N)
else:
    T_wall=np.zeros(N)
    T_wall[:]=T_amb
"""
source=np.zeros([3,N])

#initialise fields

if state=='real':
    if start_region==1:
        fields[0,:]=CP.PropsSI('D', 'P', p_amb, 'T', T_wall[:],fluid)
    elif start_region==2:
        fields[0,:]=CP.PropsSI('D', 'P', p_amb, 'Q', q_init, fluid)    #use if initial condition is two phase
else:
    fields[0,:]=p_amb/(ideal_gas.R*T_wall[:])
fields[1,:]=fields[0,:]*u_inlet
fields[2,:]=eos.get_E(fields[0,:], fields[1,:]/fields[0,:], p_amb, fluid, state=state, ideal_gas=ideal_gas)
fields[3,:]=p_amb
fields[4,:]=T_wall
fields[5,:]=fields[1,:]/fields[0,:]

#determine time step
if state=='real':
    upper_chi_est=eos.get_chi(p_inlet, rho_inlet, fluid) #chi is most likely largest at high temperature
else:
    upper_chi_est=eos.get_cs(p_inlet, rho_inlet, ideal_gas=ideal_gas)
dt=cfl*dx/upper_chi_est #maybe consider adaptive time step
nt=math.ceil(t_end/dt)
#E_correction=E_loss*dx/dt
first_point=np.zeros([6,nt])
last_state=[0,0,0,0]
p_evolution=np.zeros([3,first_point.shape[1]])
#animation=np.zeros([anim.get_animation_size(nt, 10e3),N,6])
animation=np.zeros([10121+1,N,6])
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
    
    #if i%math.floor(nt/(animation.shape[0]-1))==0 and i/math.floor(nt/(animation.shape[0]-1))<animation.shape[0]:
    #if i%nt/animation.shape[0]==0:    
    if i%divisor==0:
        pass
        """
        animation[int(i/divisor),:,0]=fields[0,:]
        animation[int(i/divisor),:,1]=fields[1,:]
        animation[int(i/divisor),:,2]=fields[2,:]
        animation[int(i/divisor),:,3]=fields[3,:]
        animation[int(i/divisor),:,4]=fields[4,:]
        animation[int(i/divisor),:,5]=fields[5,:]
        animation_time[int(i/divisor)]=i*dt
        
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,0]=fields[0,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,1]=fields[1,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,2]=fields[2,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,3]=fields[3,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,4]=fields[4,:]
        animation[int(i/math.floor(nt/(animation.shape[0]-1))),:,5]=fields[5,:]
        animation_time[int(i/math.floor(nt/(animation.shape[0]-1)))]=i*dt
        """
    if par.rank==0 and i%1e2==0:
        print("step",i)
    if (i+1)*dt>t_end:
        dt=t_end-i*dt
    #compute fluxes
    #BC are added in the same step
    #if i==15e3:
    #    flux.update_fluxes(fluxes, fields, inlet_flux, inlet, p_amb, fluid, E_correction, last_state)
    #else:
    flux.update_fluxes(fluxes, fields, inlet_flux, inlet, p_amb, par, fluid, state=state, ideal_gas=ideal_gas)
    #add source
    #src.add_momentum_source_2(source, fields, d, dx, nu)
    src.add_energy_source(source, T_boil, d, fields, T_wall)
    #src_sum[i]=np.sum(source[2,:])*dt
    #plt.plot(np.linspace(0,1,101),fluxes[0,:], label='rho*u')
    #plt.plot(np.linspace(0,1,101),fluxes[1,:], label='rho*u*u+p')
    #plt.plot(np.linspace(0.05,0.95,100),fields[3,:], label='pressure')
    #plt.legend()
    #plt.show()
    #integrate with euler explicit
    fields[:3,:]+=(dt/dx*(fluxes[:,:-1]-fluxes[:,1:])+dt*source[:,:])
    #print(eos.get_p(fields[0,:], fields[1,:]/fields[0,:], fields[2,:], fluid))
    fields[3,:]=eos.get_p(fields[0,:], fields[1,:]/fields[0,:], fields[2,:], fluid, state=state, ideal_gas=ideal_gas)
    if state=='real':
        fields[4,:]=CP.PropsSI('T', 'P', fields[3,:], 'D', fields[0,:], fluid)
    else:
        fields[4,:]=fields[3,:]/(ideal_gas.R*fields[0,:])
    fields[5,:]=fields[1,:]/fields[0,:]
    """
    for i in range(fields.shape[1]):
        if fields[3,i]==float('inf'):
            fields[3,i]=0
        #if fields[4,i]==float('inf'):
        #    fields[4,i]=T_wall[i]
        if pd.isna(fields[5,i]):
            fields[5,i]=0
    """
    #E_tot=np.sum(fields[2,:])
    #out.write_energy_loss(i, dt, E_tot)
    #E_tot_hist[i+1]=E_tot
    #dE_hist[i]=E_tot_hist[i+1]-E_tot_hist[i]
    #E_flow[i]=(fluxes[2,0]-fluxes[2,-1])*dt/dx
    #CAVEAT: No more than 1 millions steps should be written to output
    #CAVEAT: Only to be used without multithreading
    if i%divisor==0:
        out.write_step(fields, i, dt, appendix)
    """
    for x in range(N):
        if x+(par.rank-1)*N: #only to be used of all processes are of equal length
            if x+(par.rank-1)*N==100:
                x-=1
            p_evolution[math.ceil((x+(par.rank-1))/10),i]=fields[5,x]
    
    if par.rank==par.num_procs-1:
        p_evolution[0,i]=fields[5,79]
        p_evolution[1,i]=fields[5,89]
        p_evolution[2,i]=fields[5,99]
    #first_point[0,i]=fields[0,0]
    #first_point[1,i]=fields[1,0]
    #first_point[2,i]=fields[2,0]
    #first_point[3,i]=fields[3,0]
    #first_point[4,i]=fields[4,0]
    #first_point[5,i]=fields[5,0]
    """  
if par.rank==0:
    runtime=time.time()-start
    out.write_runtime(runtime, appendix)
    print("time "+str(runtime)+" sec")
#plt.plot(t,src_sum, label='src')
#plt.plot(t,dE_hist, label='dE')
#plt.plot(t, E_flow[:], label='net flux')
#plt.plot(t, E_flow[:]+src_sum[:], label='sum' )
#plt.plot(first_point[3,:])

#merge animation from different processes
animation=para.merge_results(animation, par)

if par.rank==0:
    pass
    #output data in csv table
    #out.write_all(animation,animation_time, appendix)
    """
    #create graphs for output
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,0],'result/density'+appendix+'.gif',0.4,2.5,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,1],'result/momentum'+appendix+'.gif',0.0,3.25,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,2],'result/energy'+appendix+'.gif',1.4e6,2.45e6,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,3],'result/pressure'+appendix+'.gif',99000,101000,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,4],'result/temperature'+appendix+'.gif',370,400,animation_time)
    anim.create_gif(np.linspace(0.05,0.95,100),animation.shape[0],animation[:,:,5],'result/velocity'+appendix+'.gif',-4.0,1.5,animation_time)
    """
"""
if par.rank==par.num_procs-1:
    plt.title('Velocity')
    for i in range(p_evolution.shape[0]):
        plt.plot(np.linspace(0,t_end,nt),p_evolution[i,:], label='x= '+str(i*0.1+0.8)+'m')
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
