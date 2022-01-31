#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 09:03:17 2022

@author: pascal
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import rs_trial_animation as anim

L=112 #number of data points in space (number of columns in file)
length=1 #pipe length
offset=0.5*length/L
nt=118634 #number of time steps (number of rows in file)
animation=np.zeros([nt,L,6])
animation_time=np.zeros(nt)
path='result/' #folder with .csv files
appendix='_120_to_q_0.2_w_cool_50_u_4_2500' #name of files (without data property)
filetype='.csv'

"""
read all the .csv files of one simulation and store data in numpy array animation (and animation_time)
"""

filename_path=path+'density'+appendix+filetype
with open (filename_path) as csv_file:
            csv_reader = csv.reader(csv_file)
            i=0
            for row in csv_reader:
                for j in range(len(row)):
                    animation[i,j,0]=float(row[j])
                i+=1
                
filename_path=path+'momentum'+appendix+filetype
with open (filename_path) as csv_file:
            csv_reader = csv.reader(csv_file)
            i=0
            for row in csv_reader:
                for j in range(len(row)):
                    animation[i,j,1]=float(row[j])
                i+=1
                
filename_path=path+'energy'+appendix+filetype
with open (filename_path) as csv_file:
            csv_reader = csv.reader(csv_file)
            i=0
            for row in csv_reader:
                for j in range(len(row)):
                    animation[i,j,2]=float(row[j])
                i+=1
                
filename_path=path+'pressure'+appendix+filetype
with open (filename_path) as csv_file:
            csv_reader = csv.reader(csv_file)
            i=0
            for row in csv_reader:
                for j in range(len(row)):
                    animation[i,j,3]=float(row[j])
                i+=1
                
filename_path=path+'temperature'+appendix+filetype
with open (filename_path) as csv_file:
            csv_reader = csv.reader(csv_file)
            i=0
            for row in csv_reader:
                for j in range(len(row)):
                    animation[i,j,4]=float(row[j])-273
                i+=1
                
filename_path=path+'velocity'+appendix+filetype
with open (filename_path) as csv_file:
            csv_reader = csv.reader(csv_file)
            i=0
            for row in csv_reader:
                for j in range(len(row)):
                    animation[i,j,5]=float(row[j])
                i+=1
                
filename_path=path+'time'+appendix+filetype
with open (filename_path) as csv_file:
            csv_reader = csv.reader(csv_file)
            i=0
            for row in csv_reader:
                animation_time[i]=float(row[1])
                i+=1
                
"""
create and save .gif animations
5th and 6th argument of anim.create_gif() are boundries of y-axis in plot
"""
"""
filetype='.gif'                
anim.create_gif(np.linspace(offset,1-offset,L),animation.shape[0],animation[:,:,0],path+'density'+appendix+filetype,0.5,0.65,animation_time)
anim.create_gif(np.linspace(offset,1-offset,L),animation.shape[0],animation[:,:,1],path+'momentum'+appendix+filetype,0,3.2,animation_time)
anim.create_gif(np.linspace(offset,1-offset,L),animation.shape[0],animation[:,:,2],path+'energy'+appendix+filetype,1.35e6,2.6e6,animation_time)
anim.create_gif(np.linspace(offset,1-offset,L),animation.shape[0],animation[:,:,3],path+'pressure'+appendix+filetype,99900,100100,animation_time)
anim.create_gif(np.linspace(offset,1-offset,L),animation.shape[0],animation[:,:,4],path+'temperature'+appendix+filetype,100,103,animation_time)
anim.create_gif(np.linspace(offset,1-offset,L),animation.shape[0],animation[:,:,5],path+'velocity'+appendix+filetype,0.5,1.2,animation_time)
"""

"""
plot profiles (cross sections in the (x,t) plane with consatnat t)
by loop length amount of lines per plot are specified
by changing nplot it can be zoomed in to the beginning of the simulation
e.g. nplot=10 and loop length=6 will plot profiles at 0, 10, 20, 30, 40, 50% of the final time
"""

filetype='_test.png'
nplot=5

X=np.linspace(offset, 1-offset, L)
for i in range(6):
    t=np.floor(nt/nplot)*i
    if t==animation_time.size:
        t-=1
    plt.plot(X,animation[int(t),:,0], label='t='+str(animation_time[int(t)])+' s')
plt.legend()
plt.xlabel('Location (m)')
plt.ylabel('Density (kg/m^3)')
plt.savefig(path+'density'+appendix+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    if t==animation_time.size:
        t-=1
    plt.plot(X,animation[int(t),:,1], label='t='+str(animation_time[int(t)])+' s')
plt.legend()
plt.xlabel('Location (m)')
plt.ylabel('Momentum per volume (Ns/m^3)')
plt.savefig(path+'momentum'+appendix+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    if t==animation_time.size:
        t-=1
    plt.plot(X,animation[int(t),:,2], label='t='+str(animation_time[int(t)])+' s')
plt.legend()
plt.xlabel('Location (m)')
plt.ylabel('Energy per volume (J/m^3)')
plt.savefig(path+'energy'+appendix+filetype)
#plt.show()
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    if t==animation_time.size:
        t-=1
    plt.plot(X,animation[int(t),:,3], label='t='+str(animation_time[int(t)])+' s')
plt.legend()
plt.xlabel('Location (m)')
plt.ylabel('Pressure (Pa)')
plt.savefig(path+'pressure'+appendix+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    if t==animation_time.size:
        t-=1
    plt.plot(X,animation[int(t),:,4], label='t='+str(animation_time[int(t)])+' s')
plt.legend()
plt.xlabel('Location (m)')
plt.ylabel('Temperature (°C)')
plt.savefig(path+'temperature'+appendix+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    if t==animation_time.size:
        t-=1
    plt.plot(X,animation[int(t),:,5], label='t='+str(animation_time[int(t)])+' s')
plt.legend()
plt.xlabel('Location (m)')
plt.ylabel('Velocity (m/s)')
plt.savefig(path+'velocity'+appendix+filetype)

"""
plot evolution (cross sections in the (x,t) plane with consatnat x)
by loop length amount of lines per plot are specified
with end the plot can be cut of to zoom to the transient part of the simulation
by changing nplot it can be zoomed in to the left end of the domain
e.g. nplot=10 and loop length=6 will plot evolutions at 0, 10, 20, 30, 40, 50% of the domain length from the inflow
"""

plt.clf()
end=int(nt/4)
for i in range(6):
    t=np.floor(nt/nplot)*i
    t=int(L/nplot*i)
    if t>111:
        t=111
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Density (kg/m^3)')
plt.savefig(path+'density'+appendix+'_time'+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    t=int(L/nplot*i)
    if t>111:
        t=111
    plt.plot(animation_time[:end],animation[:end,int(t),1], label='x='+str(1/nplot*i)+' m')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Momentum per volume (Ns/m^3)')
plt.savefig(path+'momentum'+appendix+'_time'+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    t=int(L/nplot*i)
    if t>111:
        t=111
    plt.plot(animation_time[:end],animation[:end,int(t),2], label='x='+str(1/nplot*i)+' m')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Energy per volume (J/m^3)')
plt.savefig(path+'energy'+appendix+'_time'+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    t=int(L/nplot*i)
    if t>111:
        t=111
    plt.plot(animation_time[:end],animation[:end,int(t),3], label='x='+str(1/nplot*i)+' m')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.savefig(path+'pressure'+appendix+'_time'+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    t=int(L/nplot*i)
    if t>111:
        t=111
    plt.plot(animation_time[:end],animation[:end,int(t),4], label='x='+str(1/nplot*i)+' m')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.savefig(path+'temperature'+appendix+'_time'+filetype)
plt.clf()
for i in range(6):
    t=np.floor(nt/nplot)*i
    t=int(L/nplot*i)
    if t>111:
        t=111
    plt.plot(animation_time[:end],animation[:end,int(t),5], label='x='+str(1/nplot*i)+' m')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.savefig(path+'velocity'+appendix+'_time'+filetype)
