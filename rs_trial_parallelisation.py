#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 19:35:08 2022

@author: pascal
"""

from mpi4py import MPI
import numpy as np

class parallel:
    #object with info for parallelisation
    #simplifies argument transfer to functions
    def __init__(self):
        self.comm=MPI.COMM_WORLD
        self.rank=self.comm.Get_rank() #current thread
        self.num_procs=self.comm.Get_size() #amount of treads

def setup_parallelisation():
    #output: int current thread and number of treads
    #not in use, replaced by class parallel
    comm=MPI.COMM_WORLD
    rank=comm.Get_rank()
    num_procs=comm.Get_size()
    return (rank, num_procs)

def get_N(N_tot, parallel):
    #split up domain to multiple threads and allocate amount of cells to a specific thread
    #inputs: int N_tot(amount of grid cells in the domain)
    #        parallel class parallel
    #output: int N (amount of grid cells for current thread)
    N=N_tot*1.0/parallel.num_procs
    if N<4:
        print("Too many threads, number of threads must be at most four times the number of cells")
        assert(False)
    remainder=N_tot%parallel.num_procs
    if parallel.rank<remainder:
        return int(N+1)
    else:
        return int(N)
    
def exchange_border(fields, par):
    #exchange field variables in border cells with neighbouring thread
    #inputs: array field(6xN array containg field variables)
    #        parallel class par
    #output: array border_upper(3, containing conservative field variables in cell to the right)
    #        array border_lower(3x2, containing conservative field variables in cells to the left)
    exchange_lower=np.zeros([3,2])
    exchange_upper=np.zeros(3)
    
    exchange_lower[0,:]=fields[0,0:2]
    exchange_lower[1,:]=fields[1,0:2]/fields[0,0:2]
    exchange_lower[2,:]=fields[3,0:2]
    
    exchange_upper[0]=fields[0,-1]
    exchange_upper[1]=fields[1,-1]/fields[0,-1]
    exchange_upper[2]=fields[3,-1]
    
    if par.rank%2==0:
        #to give order to communication even threads first send, while odd threads first wait to recieve
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
    
    #the threads most to the left and to the right have no neighbours, assign them object of type None
    if par.rank==0:
        border_lower=None
    if par.rank==par.num_procs-1:
        border_upper=None
    
    return(border_upper, border_lower)

def exchange_fluxes(fluxes, par):
    #communicate computed flux at the thread boundary to neighbouring thread
    #inputs: array fluxes (3xN, contains fluxes in thread domain)
    #        parallel object par
    #output: array flux (3, contains fluxes to the lower boundary cell of thread)
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
    #merge field variables from different threads to one array, containing all information
    #inputs: array result(ntx6xN, results within the domain of one thread)
    #        parallel class par
    #output: array result (ntx6xN_tot, results of the whole domain (of all threads))
    result=par.comm.gather(result, root=0) #rank 0 collects the array result of all the other ranks
    if par.rank==0:
        list=[]
        for i in range(par.num_procs):
            list.append(result[i])
        result=np.concatenate(list, axis=1) #merge all the different result arrays to one big array
    return result #rank 0 has gathered all the results, the other ranks result array is replaced by None type object