#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 09:54:21 2019

@author: owang
"""

# %%
from __future__ import division

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from copy import deepcopy 
import os
import glob
import pylab as py
import re
import time


# dictionary that will hold the different prof_stats 
prof_stats = {}

# legend list
lgnd = []
run_codes = [0,1]

# loop through the different run codes, load the prof_stats for each
for run_code in run_codes:
    if run_code == 0:
        # ECCO v4 R3 original
        output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/h8/h8i_i48/python_output/prof_cost_comparisons/'
        output_name = 'ECCOv4r3_original_mon'
    
        
    elif run_code == 1:
        output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/r042/r042f/run.v4_rls2.042ff10.iter115/python_output/prof_cost_comparisons/'
        output_name = 'ECCOv4r3_extension_it115_mon'
        
        
    else:
        print('not a valid run code ', run_code)
        sys.exit()
        
   
    os.chdir(output_dir)
    prof_stats[run_code] = np.load(output_name + '.npy').item()
    
    lgnd.append(prof_stats[run_code]['run_name'])
        
    #%%
min_year = 1992
max_year = 2018


# a counter of the total number of figures that we've made     
total_figures = 0

# a dictionary matching profile type with figure numbers for T and S
T_fig_nums = {}
S_fig_nums = {}


# loop through the different prof_stats objects
costs_keys = []

for run_code in run_codes:
    # pull the keys out for the different profile types for T and S
    costs_keys = (set(prof_stats[run_code]['costs_T'].keys()) | set(costs_keys))
    costs_keys = (set(prof_stats[run_code]['costs_S'].keys()) | set(costs_keys))
    
print costs_keys
#%%    
# start the plotting by closing all existing plots
plt.close('all')

for run_code in run_codes:
    
    years = np.array(prof_stats[run_code]['years'])
    num_months=12
    months=range(1,num_months+1)
    num_years=len(years)
    yearsinmon = np.zeros((num_years,num_months))
    tot_mon=num_years*num_months

    for yi in range(num_years):
        for mi in range(num_months):
            yearsinmon[yi,mi]=years[yi]+mi/12
        
    # loop through each profile type that has a T cost
    #%%
    for prof_type in costs_keys:
        print ('PROFILE TYPE: ', prof_type)
        
        if prof_type in prof_stats[run_code]['costs_T'].keys():
            cost_T  = prof_stats[run_code]['costs_T'][prof_type]
            count_T = prof_stats[run_code]['count_T'][prof_type]
            #print ('count_T: ', count_T)
            #if count_T != 0:
            cpd_T   = cost_T / count_T

        else:
            cost_T = yearsinmon*0
            count_T = yearsinmon*0
            cpd_T = yearsinmon*0
        
        # see if we have a figure associated with this profile type
        # already (we probably will if we have more than one run code)
        # if so, get the figure number
        # if not, make a new figure number
        if prof_type in T_fig_nums:
            fig_num = T_fig_nums[prof_type]
        else:
            total_figures += 1
            fig_num = total_figures
            # add this figure number to the dictionary
            T_fig_nums[prof_type] = fig_num
        
        #go to figure number fig_num
        plt.figure(fig_num)
        # make 3 subplots
        ax=plt.subplot(3, 1, 1)
        # first plot cost
        plt.plot(np.reshape(yearsinmon,tot_mon), np.reshape(cost_T,tot_mon))
        plt.grid(True)
        plt.ylabel('cost')
        plt.title('Temperature : ' + prof_type)
        plt.legend(lgnd, loc='best')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.set_xlim([min_year, max_year])
        
        # then plot count
        ax=plt.subplot(3,1,2)
        plt.plot(np.reshape(yearsinmon,tot_mon), np.reshape(count_T,tot_mon))
        plt.grid(True)
        plt.ylabel('count')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.set_xlim([min_year, max_year])
        
        # then plot normalized cost
        ax=plt.subplot(3,1,3)
        plt.plot(np.reshape(yearsinmon,tot_mon), np.reshape(cpd_T,tot_mon))
        plt.ylabel('T cost per datum')
        ax.set_xlim([min_year, max_year])
        #[a,b] = ax.get_ylim()
        #ax.set_ylim([0,b])
        #print b
        
        plt.grid(True)
        plt.show()
        py.savefig(output_dir + '/fig/' + 'T_'+prof_type+'.png') 

    
        # repeat as above for salinity
        if prof_type in prof_stats[run_code]['costs_S'].keys():
            cost_S  = prof_stats[run_code]['costs_S'][prof_type]
            count_S = prof_stats[run_code]['count_S'][prof_type]
            #print ('count_S: ', count_S)
            #if count_S != 0:
            cpd_S   = cost_S / count_S

        else:
            cost_S = yearsinmon*0
            count_S = yearsinmon*0
            cpd_S = yearsinmon*0
        
        # see if we have a figure associated with this profile type
        # already (we probably will if we have more than one run code)
        # if so, get the figure number
        # if not, make a new figure number
        if prof_type in S_fig_nums:
            fig_num = S_fig_nums[prof_type]
        else:
            total_figures += 1
            fig_num = total_figures
            # add this figure number to the dictionary
            S_fig_nums[prof_type] = fig_num
        
        #go to figure number fig_num
        plt.figure(fig_num)
        # make 3 subplots
        ax=plt.subplot(3, 1, 1)
        # first plot cost
        plt.plot(np.reshape(yearsinmon,tot_mon), np.reshape(cost_S,tot_mon))
        plt.grid(True)
        plt.ylabel('cost')
        plt.title('Salinity : ' + prof_type)
        plt.legend(lgnd, loc='best')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.set_xlim([min_year, max_year])
        
        # then plot count
        ax=plt.subplot(3,1,2)
        plt.plot(np.reshape(yearsinmon,tot_mon), np.reshape(count_S,tot_mon))
        plt.grid(True)
        plt.ylabel('count')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.set_xlim([min_year, max_year])
        
        # then plot normalized cost
        ax=plt.subplot(3,1,3)
        plt.plot(np.reshape(yearsinmon,tot_mon), np.reshape(cpd_S,tot_mon))
        plt.ylabel('S cost per datum')
        ax.set_xlim([min_year, max_year])
        #[a,b] = ax.get_ylim()
        #ax.set_ylim([0,b])
        #print b
        
        plt.grid(True)
        plt.show()
        py.savefig(output_dir + '/fig/' + 'S_'+prof_type+'.png') 
        
        

