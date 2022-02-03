#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 15:20:13 2021

@author: pascal
"""
import matplotlib.animation as anim
import matplotlib.pyplot as plt
import math

def get_animation_size(nt,min_size):
    #currently not in use
    #evaluate array size to store all data
    #inputs: floats nt (amount of time steps), min_size (desired amount of rows in output file)
    #output: int array size
    if nt<=min_size:
        return nt
    ratio=float(nt/min_size)
    return int(nt/math.floor(ratio))

def create_animation(i, x, u, y_min,y_max,t):
    #function is called by anim.FuncAnimation to create one frame in animation
    #inputs: int i (counting variable, time_step to be plotted in current frame)
    #       array x(vector contining domain points),u (array containng variable values at all time steps, u.shape[0] must be >= i and u.shape[1]=x.shape)
    #       float y_min,y_max (plot boundaries of y-axis), t (time)
    #NOTE: plots are intended for a domain with length 1 (otherwise change plt.xlim)
    plt.cla()
    plt.xlim(0,1)
    plt.ylim(y_min,y_max)
    plt.plot(x,u[i,:],label=t[i])
    plt.legend()
    
def create_gif(x,nt,u,fname,y_min,y_max,t):
    #create animated plots and save as .gif
    #inputs: arrays x (vector containing domain points), u (array containing variable values at all time steps, u.shape[0] must be >= i and u.shape[1]=x.shape))
    #        array t (vector containing time information, t.size=u.shate[0])
    #        int nt (amount of time steps)
    #        float y_min,y_max(confing y-axis of plot)
    #        string fname (filename of saved plot)
    #output: void (animation is saved on machine)
    fig=plt.figure()
    animation=anim.FuncAnimation(fig,create_animation, fargs=(x,u,y_min,y_max,t),frames=range(0,nt,500),blit=False, interval=1)
    animation.save(fname, dpi=80, writer='Pillow')
