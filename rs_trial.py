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


#initialise paralellisation
par=para.parallel()

#spatial and temporal dimensions
L=1 #pipe length
dx=0.009 #grid spacing
N_tot=math.ceil(L/dx) #number of cells
N=para.get_N(N_tot, par) #allcoate number of cells to threads
d=0.1 #pipe diameter
t_end=1.5 #end time
cfl=0.4 #cfl number

 #fluid and flow initialisation
state='real'  #real or ideal
start_region=2 #type 1 if initial condition is single phase, or 2 if initial condition is two phase, only considered for real gases
assert(state=='real' or state=='ideal')
assert(start_region==1 or start_region==2)
fluid='Water' #fluid must be known for CoolProp (real gas only)
gamma=1.4   #ideal gas only
R=287.058   #ideal gas only
if state=='ideal':
    ideal_gas=eos.ideal_gas(gamma, R) #create object with ideal gas properties
else:
    ideal_gas=None
#nu=15e-6 #kinematic viscosity, required only if wall shear stress is included
T_inlet=120+273 #inlet temperature (K)
T_amb=20+273 #initial temperature(K), considered if start_region==1
undercool=50 # if start_region=2, the wall needs to be colder than boiling temperature to extract latent heat, specify here how much colder (K) than boiling temperature
q_init=0.2 #initial steam quality, considered if start_region==2
p_amb=1e5 #initial pressure (Pa)
u_inlet=3 #inlet (and initial) velocity (m/s)
p_inlet=1*p_amb #inlet pressure (Pa)
p_outlet=p_amb #outlet pressure (Pa)
#E_loss=0 #energy extracted by the environment, to be determined manually, only used for energy optimised boundary conditions in rs_trial_fluxes.py
if state=='real':
    rho_inlet=CP.PropsSI('D', 'T',T_inlet, 'P',p_inlet,fluid)
    T_boil=CP.PropsSI('T', 'P', p_amb, 'Q', 0.5,fluid) #find boiling temperature (Q value is arbitrary number between (excluding) 0 and 1)
    if start_region==2:
        T_amb=T_boil
else:
    rho_inlet=p_inlet/(ideal_gas.R*T_inlet)
    T_boil=None
    
#output file initialisation
divisor=1 #every n-th step will be written to the output file (if more than ca. 1,000,000 steps are written, the file cannot be opened by excell or similar)
appendix='test' #part of the filename to specify the simulation

#combine inlet conditions to an array
inlet=np.array([rho_inlet,u_inlet,p_inlet,T_inlet])
inlet_flux=flux.calc_inlet_bc(rho_inlet, u_inlet, p_inlet, fluid,state=state, ideal_gas=ideal_gas)


#initialise wall temperature
T_wall=np.zeros(N) #wall temerature
if start_region==1:
    T_wall[:]=T_amb
elif start_region==2:
    T_wall[:]=T_boil-undercool
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



#initialise data arrays
fields=np.zeros([6, N]) #0 stores density, 1 stores momentum, 2 stores energy, 3 stores pressure, 4 stores temperature(, 5 stores velocity)
fluxes=np.zeros([3, N+1]) #0 stores density flux, 1 stores momentum flux, 2 stores energy flux
source=np.zeros([3,N]) #0 stores density source, 1 stores momentum source, 2 stores energy source

if state=='real': #initial density
    if start_region==1:
        fields[0,:]=CP.PropsSI('D', 'P', p_amb, 'T', T_wall[:],fluid)
    elif start_region==2:
        fields[0,:]=CP.PropsSI('D', 'P', p_amb, 'Q', q_init, fluid)    #use if initial condition is two phase
else: #initial density
    fields[0,:]=p_amb/(ideal_gas.R*T_wall[:])
fields[1,:]=fields[0,:]*u_inlet #initial momentum
fields[2,:]=eos.get_E(fields[0,:], fields[1,:]/fields[0,:], p_amb, fluid, state=state, ideal_gas=ideal_gas) #initial energy
fields[3,:]=p_amb #initial pressure
fields[4,:]=T_wall #initial temperature
fields[5,:]=fields[1,:]/fields[0,:] #initial velocity

#determine time step
if state=='real':
    upper_chi_est=eos.get_chi(p_inlet, rho_inlet, fluid) #acoustic wave speed is largest at high temperature
else:
    upper_chi_est=eos.get_cs(p_inlet, rho_inlet, ideal_gas=ideal_gas) #speed of sound is largest at high temperature
dt=cfl*dx/upper_chi_est # time step
nt=math.ceil(t_end/dt)  # number of time steps
#E_correction=E_loss*dx/dt #energy correction, only used for energy optimised boundary conditions in rs_trial_fluxes.py
p_hist=np.zeros(nt) #can be used to detect pressure maxima in time
animation_size=math.ceil(nt/divisor)
animation=np.zeros([animation_size,N,6])
animation_time=np.zeros(animation.shape[0])

if par.rank==0:
    start=time.time()
    
for i in range(nt):
   
    if i%divisor==0:
        
        animation[int(i/divisor),:,0]=fields[0,:]
        animation[int(i/divisor),:,1]=fields[1,:]
        animation[int(i/divisor),:,2]=fields[2,:]
        animation[int(i/divisor),:,3]=fields[3,:]
        animation[int(i/divisor),:,4]=fields[4,:]
        animation[int(i/divisor),:,5]=fields[5,:]
        animation_time[int(i/divisor)]=i*dt
        
    if par.rank==0 and i%1e2==0:
        print("step",i)
    if (i+1)*dt>t_end: #adapt the last time step to meed the final time
        dt=t_end-i*dt
    #compute fluxes and boundary conditions
    flux.update_fluxes(fluxes, fields, inlet_flux, inlet, p_outlet, par, fluid, state=state, ideal_gas=ideal_gas)
    #add source
    #src.add_momentum_source_2(source, fields, d, dx, nu) #wall friction
    src.add_energy_source(source, T_boil, d, fields, T_wall) #energy extraction
    #integrate with euler explicit
    fields[:3,:]+=(dt/dx*(fluxes[:,:-1]-fluxes[:,1:])+dt*source[:,:])
    #calculate other field variables based on equation of state
    fields[3,:]=eos.get_p(fields[0,:], fields[1,:]/fields[0,:], fields[2,:], fluid, state=state, ideal_gas=ideal_gas)
    if state=='real':
        fields[4,:]=CP.PropsSI('T', 'P', fields[3,:], 'D', fields[0,:], fluid)
    else:
        fields[4,:]=fields[3,:]/(ideal_gas.R*fields[0,:])
    fields[5,:]=fields[1,:]/fields[0,:]

    """
    #CAVEAT: No more than 1 millions steps should be written to output
    #CAVEAT: Only to be used without multithreading
    #write step directly to output file
    if i%divisor==0:
        out.write_step(fields, i, dt, appendix)
    """
if par.rank==0:
    runtime=time.time()-start
    out.write_runtime(runtime, appendix)
    print("runtime "+str(runtime)+" sec")
#plt.plot(first_point[3,:])

#merge animation from different processes
animation=para.merge_results(animation, par)

if par.rank==0:
    #output data in csv table
    #CAVEAT: if steps are written individually after step (only for single thread simulations) the next line should be commented out
    out.write_all(animation,animation_time, appendix)
    #pass #if line above is commented out, uncomment this line
    
"""
#extract and plot pressure maxima and minima (x axis has no physical meaning)
max_pressure=[]
min_pressure=[]
for i in range(1,p_hist.size-1):
    if p_hist[i]>p_hist[i-1] and p_hist[i]>p_hist[i+1]:
        max_pressure.append(p_hist[i])
    if p_hist[i]<p_hist[i-1] and p_hist[i]<p_hist[i+1]:
        min_pressure.append(p_hist[i])

plt.plot(max_pressure)
plt.savefig('result/max_pressure')
plt.clf()
plt.plot(min_pressure)
plt.savefig('result/min_pressure')
plt.clf()
"""

