#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 19:35:08 2022

@author: pascal
"""

from mpi4py import MPI
import numpy as np

class parallel:
    def __init__(self):
        self.comm=MPI.COMM_WORLD
        self.rank=self.comm.Get_rank()
        self.num_procs=self.comm.Get_size()

def setup_parallelisation():
    comm=MPI.COMM_WORLD
    rank=comm.Get_rank()
    num_procs=comm.Get_size()
    return (rank, num_procs)

def get_N(N_tot, parallel):
    N=N_tot*1.0/parallel.num_procs
    if N<2:
        print("Too many threads, number of threads must be at most half the number of cells")
        assert(False)
    remainder=N_tot%parallel.num_procs
    if parallel.rank<remainder:
        return int(N+1)
    else:
        return int(N)
    
def exchange_border(fields, par):
    exchange_lower=np.zeros([3,2])
    exchange_upper=np.zeros(3)
    
    exchange_lower[0,:]=fields[0,0:2]
    exchange_lower[1,:]=fields[1,0:2]/fields[0,0:2]
    exchange_lower[2,:]=fields[3,0:2]
    
    exchange_upper[0]=fields[0,-1]
    exchange_upper[1]=fields[1,-1]/fields[0,-1]
    exchange_upper[2]=fields[3,-1]
    
    if par.rank%2==0:
        if par.rank!=0:
            par.comm.send(exchange_lower, dest=par.rank-1, tag=par.rank*10+1)
        if par.rank!=par.num_procs-1:
            border_upper=par.comm.recv(source=par.rank+1, tag=(par.rank+1)*10+1)
            par.comm.send(exchange_upper, dest=par.rank+1, tag=par.rank*10+2)
        if par.rank!=0:
            border_lower=par.comm.recv(source=par.rank-1, tag=(par.rank-1)*10+2)
    else:
        if par.rank != par.num_procs-1:
            border_upper=par.comm.recv(source=par.rank+1, tag=(par.rank+1)*10+1)
        par.comm.send(exchange_lower, dest=par.rank-1, tag=par.rank*10+1)
        border_lower=par.comm.recv(source=par.rank-1, tag=(par.rank-1)*10+2)
        if par.rank!=par.num_procs-1:
            par.comm.send(exchange_upper, dest=par.rank+1, tag=par.rank*10+2)
            
    if par.rank==0:
        border_lower=None
    if par.rank==par.num_procs-1:
        border_upper=None
    
    return(border_upper, border_lower)

def exchange_fluxes(fluxes, par):
    if par.rank%2==0:
        if par.rank<par.num_procs-1:
            par.comm.send(fluxes[:,-1], dest=par.rank+1, tag=par.rank*10+3)
        if par.rank!=0:
            flux=par.comm.recv(source=par.rank-1, tag=(par.rank-1)*10+3)
        else:
            flux=fluxes[:,0]
    else:
        flux=par.comm.recv(source=par.rank-1, tag=(par.rank-1)*10+3)
        if par.rank<par.num_procs-1:
            par.comm.send(fluxes[:,-1], dest=par.rank+1, tag=par.rank*10+3)
        
    return flux

def merge_results(result, par):
    result=par.comm.gather(result, root=0)
    if par.rank==0:
        list=[]
        for i in range(par.num_procs):
            list.append(result[i])
        result=np.concatenate(list, axis=1)
    return result