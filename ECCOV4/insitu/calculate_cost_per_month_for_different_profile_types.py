#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:50:35 2019

@author: ifenty
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
#%%


#run_code = int(sys.argv[1])
run_code = 0

print(' run_code = ', run_code)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%

# default parameters
remove_area_gamma = 0
remove_spatial_scaling_factor = 0
   
y1 = 1992
y2 = 2018

if run_code == 0:
    # ECCO v4 R3 original
    prof_dir= '/mnt/intraid/ecco-intrnl/data01/ecco1/data4/Version4/Release3/profiles/'
    output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/h8/h8i_i48/python_output/prof_cost_comparisons/'
    output_name = 'ECCOv4r3_original_mon'

    apb = ['seals']
    argo =['argo']
    ctd = ['ctd','climode','ices']
    itp = ['itp']
    mooring = ['*mooring']
    xbt = ['xbt']
       
    # make a dictionary, one entry for each profile type
    profs = {}
    profs['apb'] = apb
    profs['argo'] = argo
    profs['ctd'] = ctd
    profs['itp'] = itp
    profs['mooring'] = mooring
    profs['xbt'] = xbt
    
    # set parameters for this particular run code
    remove_area_gamma = 1
    remove_spatial_scaling_factor = 0
    
    y1 = 1992
    y2 = 2018
    
elif run_code == 1:
    print ('not yet')
    prof_dir= '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/JPL/r4/release2/r042/r042f/run.v4_rls2.042ff10.iter115/profiles/'
    output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/r042/r042f/run.v4_rls2.042ff10.iter115/python_output/prof_cost_comparisons/'
    output_name = 'ECCOv4r3_extension_it115_mon'
    
    apb = ['APB']
    argo =['ARGO']
    ctd = ['CTD','GOSHIP','CCHDO']
    itp = ['ITP']
    mooring = ['MRB']
    xbt = ['XBT']
    glider = ['GLD']
       
    # make a dictionary, one entry for each profile type
    profs = {}
    profs['apb'] = apb
    profs['argo'] = argo
    profs['ctd'] = ctd
    profs['glider'] = glider
    profs['itp'] = itp
    profs['mooring'] = mooring
    profs['xbt'] = xbt

    # set parameters for this particular run code
    remove_area_gamma = 1
    remove_spatial_scaling_factor = 1
    
    y1 = 1992
    y2 = 2018
else:
    print('not a valid run code ', run_code)
    sys.exit()
# %%

all_profs_costs_T = {}
all_profs_count_T = {}
all_profs_costs_S = {}
all_profs_count_S = {}
num_months=12
months=range(1,num_months+1)

#%%
for prof_type in profs.keys():
    print ('PROFILE TYPE: ', prof_type)
    
    os.chdir(prof_dir)

    prof_type_file_globs = profs[prof_type]

    years = range(y1,y2+1)
    num_years = len(years)
    
    prof_costs_T = np.zeros((num_years,12,1))
    prof_costs_S = np.zeros((num_years,12,1))
    prof_count_T = np.zeros((num_years,12,1))
    prof_count_S = np.zeros((num_years,12,1))
        
    # loop through the different filenames that 
    # contain files of type 'prof_type'
    for ptfgi in prof_type_file_globs:
        print ptfgi

        files = glob.glob(ptfgi + "*nc")
        
        
        for f in files:
            tmp = xr.open_dataset(f)
            # by default we do not assume that the profile file
            # has T or S data
            has_T_data = 0;
            has_S_data = 0;
            
            # loop through the fields and look to see if the prof_T 
            # or prof_S key exists
            for k in tmp.data_vars.keys():
                if 'prof_T' in k:
                    has_T_data = 1
                if 'prof_S' in k:
                    has_S_data = 1
                    
            # find the years of each profile in this file
            prof_years = np.floor(tmp.prof_YYYYMMDD.data / 1e4)
            prof_months = np.floor((tmp.prof_YYYYMMDD.data-prof_years*1e4) / 1e2)
      
            # loop through all of the years under consideration
            for yi in range(num_years):
                
                # pull the current year
                cur_year = years[yi]
     
                for mi in range(num_months):
                    cur_month=months[mi]
                    print('cur_month = ' + str(cur_month))
                    # find the indices of the profiles that fall within cur_yearmpiprof.test01.out
                    ix= np.isin(prof_years*1e2+prof_months*1, cur_year*1e2+cur_month)        
            
                    # continue if at least one profile is within this cur_year
                    if sum(ix) > 0:
                        print('found profiles of year,month ' + str(cur_year) + ' ' + str(cur_month) + ' in file ' + f)
            
                        if has_T_data:
                            # find the number of profiles within this year 
                            # that have temperature weight > 0
                            T_data_count = len(np.where(tmp.prof_Tweight[ix,:] > 0)[0])
                            
                            # model data difference
                            T_diff = tmp.prof_Testim[ix,:] - tmp.prof_T[ix,:]
                            
                            T_weight = tmp.prof_Tweight[ix,:]
                            
                            # if necessary, remove the gamma factor
                            if remove_area_gamma:
                                T_weight = T_weight / tmp.prof_area_gamma[ix]
                            # if necessary, remove the spatial scaling factor
                            if remove_spatial_scaling_factor:
                                T_weight = T_weight / tmp.prof_spatial_scaling_factor[ix]
                                
                            T_cost = (T_diff **2) * T_weight
                                
                            # accumulate costs and counts to prof_cost and count       
                          
                            prof_costs_T[yi,mi] += np.nansum(T_cost)
                            prof_count_T[yi,mi] += T_data_count
                            
                        if has_S_data:
                                
                            S_data_count = len(np.where(tmp.prof_Sweight[ix,:] > 0)[0])
            
                            S_diff = tmp.prof_Sestim[ix,:] - tmp.prof_S[ix,:]
                            S_weight = tmp.prof_Sweight[ix,:]
                            
                            # if necessary, remove the gamma factor
                            if remove_area_gamma:
                                S_weight = S_weight / tmp.prof_area_gamma[ix]
                            # if necessary, remove the spatial scaling factor
                            if remove_spatial_scaling_factor:
                                S_weight = S_weight / tmp.prof_spatial_scaling_factor[ix]
                                
                            S_cost = (S_diff **2) * S_weight
                                
                            prof_costs_S[yi,mi] += np.nansum(S_cost)
                            prof_count_S[yi,mi] += S_data_count

    #%%    
    all_profs_costs_S[prof_type] = prof_costs_S
    all_profs_costs_T[prof_type] = prof_costs_T
    all_profs_count_S[prof_type] = prof_count_S
    all_profs_count_T[prof_type] = prof_count_T

#%%
prof_stats = {}   
prof_stats['costs_S'] = all_profs_costs_S
prof_stats['costs_T'] = all_profs_costs_T

prof_stats['count_S'] = all_profs_count_S
prof_stats['count_T'] = all_profs_count_T

prof_stats['years'] = years
prof_stats['prof_dir'] = prof_dir
prof_stats['run_name'] = output_name

try:
    os.mkdir(output_dir)
except OSError:
    print('directory exists ', output_dir)
    
print('saving ', output_name)
print('to ', output_dir)
    
np.save(output_dir + '/' + output_name, prof_stats)
