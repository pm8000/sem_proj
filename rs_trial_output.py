#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:40:20 2021

@author: pascal
"""

import csv

def write_step(fields, i, dt, appendix):
    #only to be performed when runing on single core
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
            
def write_all(data, dt, appendix):
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
                visualwriter.writerow(['time'] + [dt])
                
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
                visualwriter.writerow(['time'] + [(i+1)*dt])
                
def write_runtime(runtime,appendix):
    fname='result/runtime_'+appendix+'.txt'            
    with open(fname, 'w') as f:
        f.write(str(runtime)+' sec')
        