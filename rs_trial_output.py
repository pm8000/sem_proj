#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:40:20 2021

@author: pascal
"""

import csv

def write_step(fields, i, dt, appendix):
    #CAVEAT: only to be performed when running on single core
    #current step is written/appended to .csv output file
    #files are saved in folder result (make sure folder exists)
    #inputs: array fields(6xN, contains all field variables at current time step)
    #        int i (time step)
    #        float dt (time step size)
    #        string appendix(part of filename to label simulation)
    #output: void (.csv files are saved/updated in folder result)
    if i==0:
        fname='result/density_'+appendix+'.csv'
        with open(fname, mode='w', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[0,:])
        fname='result/momentum_'+appendix+'.csv'
        with open(fname, mode='w', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[1,:])
        fname='result/energy_'+appendix+'.csv'
        with open(fname, mode='w', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[2,:])
        fname='result/pressure_'+appendix+'.csv'
        with open(fname, mode='w', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[3,:])
        fname='result/temperature_'+appendix+'.csv'
        with open(fname, mode='w', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[4,:])
        fname='result/velocity_'+appendix+'.csv'
        with open(fname, mode='w', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[5,:])
        fname='result/time_'+appendix+'.csv'
        with open(fname, mode='w', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(['time'] + [(i+1)*dt])
    else:
        fname='result/density_'+appendix+'.csv'
        with open(fname, mode='a+', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[0,:])
        fname='result/momentum_'+appendix+'.csv'
        with open(fname, mode='a+', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[1,:])
        fname='result/energy_'+appendix+'.csv'
        with open(fname, mode='a+', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[2,:])
        fname='result/pressure_'+appendix+'.csv'
        with open(fname, mode='a+', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[3,:])
        fname='result/temperature_'+appendix+'.csv'
        with open(fname, mode='a+', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[4,:])
        fname='result/velocity_'+appendix+'.csv'
        with open(fname, mode='a+', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(fields[5,:])
        fname='result/time_'+appendix+'.csv'
        with open(fname, mode='a+', newline='') as csvfile:
            visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            visualwriter.writerow(['time'] + [(i+1)*dt])
            
def write_all(data, timestamp, appendix):
    #CAVEAT: if more than approx. 1 million time steps are written the .csv cannot be opened with MS Excel or simillar
    # write all data to a .csv output file
    #inputs: array data (ntxNx6, data to be written to output file)
    #        array timestamp (timestamp.shape=data.shape[0], contains time to allocate data values to a certain time in plot_data.py)
    #        string appendix(part of filename to label simulation)
    #output: void (.csv files are saved in folder result)
    for i in range(data.shape[0]):
        if i==0:
            fname='result/density_'+appendix+'.csv'
            with open(fname, mode='w', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[0,:,0])
            fname='result/momentum_'+appendix+'.csv'
            with open(fname, mode='w', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[0,:,1])
            fname='result/energy_'+appendix+'.csv'
            with open(fname, mode='w', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[0,:,2])
            fname='result/pressure_'+appendix+'.csv'
            with open(fname, mode='w', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[0,:,3])
            fname='result/temperature_'+appendix+'.csv'
            with open(fname, mode='w', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[0,:,4])
            fname='result/velocity_'+appendix+'.csv'
            with open(fname, mode='w', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[0,:,5])
            fname='result/time_'+appendix+'.csv'
            with open(fname, mode='w', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(['time'] + [timestamp[i]])
                
        else:
            fname='result/density_'+appendix+'.csv'
            with open(fname, mode='a+', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[i,:,0])
            fname='result/momentum_'+appendix+'.csv'
            with open(fname, mode='a+', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[i,:,1])
            fname='result/energy_'+appendix+'.csv'
            with open(fname, mode='a+', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[i,:,2])
            fname='result/pressure_'+appendix+'.csv'
            with open(fname, mode='a+', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[i,:,3])
            fname='result/temperature_'+appendix+'.csv'
            with open(fname, mode='a+', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[i,:,4])
            fname='result/velocity_'+appendix+'.csv'
            with open(fname, mode='a+', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(data[i,:,5])
            fname='result/time_'+appendix+'.csv'
            with open(fname, mode='a+', newline='') as csvfile:
                visualwriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                visualwriter.writerow(['time'] + [timestamp[i]])
                
def write_runtime(runtime,appendix):
    #write the run time (time in loop) to a .txt file and save in result folder (folder must exist)
    #inputs: float runtime
    #        string appendix(part of filename to label simulation)
    #output: void (.csv files are saved/updated in folder result)
    fname='result/runtime_'+appendix+'.txt'            
    with open(fname, 'w') as f:
        f.write(str(runtime)+' sec')
        