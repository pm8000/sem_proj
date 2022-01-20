#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 15:20:13 2021

@author: pascal
"""
import matplotlib.animation as anim
import matplotlib.pyplot as plt

def create_animation(i, x, u, y_min,y_max,t):
    plt.cla()
    plt.xlim(0,1)
    plt.ylim(y_min,y_max)
    plt.plot(x,u[i,:],label=t[i])
    plt.legend()
    
def create_gif(x,nt,u,fname,y_min,y_max,t):
    fig=plt.figure()
    animation=anim.FuncAnimation(fig,create_animation, fargs=(x,u,y_min,y_max,t),frames=range(0,nt,200),blit=False, interval=1)
    animation.save(fname, dpi=80, writer='Pillow')
