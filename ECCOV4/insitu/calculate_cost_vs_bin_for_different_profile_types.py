#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 08:58:52 2019

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
import string as str
#%%


#run_code = int(sys.argv[1])
run_code = 2

print(' run_code = ', run_code)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%

# default parameters
remove_area_gamma = 0
remove_spatial_scaling_factor = 0
   
y1 = 1992
y2 = 2018

if run_code == 0:
    # ECCO v4 R3 original
    profdata_dir='/mnt/intraid/ecco-intrnl/data07/owang/Ecco_data/v4/JPL/r4/release2/input/data/profiles/ECCO_v4r2/llc90_20160308/'
    prof_dir= '/mnt/intraid/ecco-intrnl/data01/ecco1/data4/Version4/Release3/profiles/'
    output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/h8/h8i_i48/python_output/prof_cost_comparisons/'
    output_name = 'ECCOv4r3_original_bin'

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

elif run_code == 1 or run_code == 2:
    print ('not yet')
    profdata_dir= '//mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/JPL/r4/release3/input/data/profiles/ECCO_v4r2/20190131/llc90/post-spatial-scaling/merged/output/'

    if run_code == 1:
        # iter 115
        prof_dir= '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/JPL/r4/release2/r042/r042f/run.v4_rls2.042ff10.iter115/profiles/'
        output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/r042/r042f/run.v4_rls2.042ff10.iter115/python_output/prof_cost_comparisons/'
        output_name = 'ECCOv4r3_extension_it115'
    elif run_code == 2:
    # iter 118
        prof_dir= '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/JPL/r4/release2/r042/r042f/run.v4_rls2.042ff10.iter118/profiles/'
        output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/r042/r042f/run.v4_rls2.042ff10.iter118/python_output/prof_cost_comparisons/'
        output_name = 'ECCOv4r3_extension_it118'
        
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
nbins=10242

#%%
for prof_type in profs.keys():
    print ('PROFILE TYPE: ', prof_type)
   
    os.chdir(prof_dir)
    prof_type_file_globs = profs[prof_type]

    #prof_type_file_globs = profs[prof_type]

    print prof_type_file_globs
    print "==========="
    prof_costs_bin_T = np.zeros((nbins,1))
    prof_costs_bin_S = np.zeros((nbins,1))
    prof_count_bin_T = np.zeros((nbins,1))
    prof_count_bin_S = np.zeros((nbins,1))      
    # loop through the different filenames that 
    # contain files of type 'prof_type'
    for ptfgi in prof_type_file_globs:
     #   print "test ====="
 
        files = glob.glob(ptfgi + "*nc")
        print files

        #tmpstr=str.replace(files,"_model.nc","")
        #filesdata=glob.glob(profdata_dir+tmpstr + "*nc")
        #print filesdata

       # print files    
        #%%
        years = range(y1,y2+1)
        num_years = len(years)
        

        #%%
        
        for ifidx in range(len(files)):
        
            #fdata=filesdata[ifidx] 
            f=files[ifidx]
        
            tmpstr=str.replace(f,"_model.nc","")
            fdata=profdata_dir+tmpstr + ".nc"
            print '======================='
            
            print f
            print tmpstr
            print profdata_dir+tmpstr + ".nc"
            print fdata
            
            #%%
            
            tmpdata = xr.open_dataset(fdata)
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
            #prof_months = np.floor((tmp.prof_YYYYMMDD.data-prof_years*1e4) / 1e2)
            prof_bins = tmpdata.prof_bin_id_a.data

            

#%%

            ix=np.ones(len(prof_years),dtype=bool)
            ix.shape
            #%%
            
            # continue if at least one profile is within this cur_year
            if sum(ix) > 0:
                #print('found profiles of year,month ' + str(cur_year) + ' ' + str(cur_month) + ' in file ' + f)
    
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
                    for ibin in range(nbins):
                              
                        # pull the current year
                        cur_bin = ibin+1
                        if (cur_bin % 5000)==1:
                            print cur_bin
                        #for mi in range(num_months):
                        #cur_month=months[mi]
                        #print('cur_bin = ' + str(cur_bin))
                        # find the indices of the profiles that fall within cur_yearmpiprof.test01.out
        #                   ix= np.isin(prof_years*1e2+prof_months*1, cur_year*1e2+cur_month)      
                        #idix= np.isin(prof_bins, cur_bin)    
                        #ixtmp= np.logical_and(prof_years>=y1,prof_years<y2)
                        #ix=np.logical_and(prof_bins==(ibin+1),ixtmp)
                        ix1=np.equal(prof_bins,(cur_bin))
                        
                        #ix=np.logical_and(prof_bins==15000,ixtmp)
        
                        ix1.shape
                        
                        prof_costs_bin_T[ibin] += np.nansum(T_cost[ix1])
                        prof_count_bin_T[ibin] += len(np.where(T_cost[ix1,:] > 0)[0])
                    
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
                        
                    S_cost = (S_diff **2) * S_weight                    # accumulate costs and counts to prof_cost and count       
                    for ibin in range(nbins):
                              
                        # pull the current year
                        cur_bin = ibin+1
                        if (cur_bin % 5000)==1:
                            print cur_bin
                        #for mi in range(num_months):
                        #cur_month=months[mi]
                        #print('cur_bin = ' + str(cur_bin))
                        # find the indices of the profiles that fall within cur_yearmpiprof.test01.out
        #                   ix= np.isin(prof_years*1e2+prof_months*1, cur_year*1e2+cur_month)      
                        #idix= np.isin(prof_bins, cur_bin)    
                        #ixtmp= np.logical_and(prof_years>=y1,prof_years<y2)
                        #ix=np.logical_and(prof_bins==(ibin+1),ixtmp)
                        ix1=np.equal(prof_bins,cur_bin)
        
                        #ix=np.logical_and(prof_bins==15000,ixtmp)
        
                        ix1.shape
                        
                        prof_costs_bin_S[ibin] += np.nansum(S_cost[ix1])
                        prof_count_bin_S[ibin] += len(np.where(S_cost[ix1,:] > 0)[0])
                    
    #%%    
    all_profs_costs_S[prof_type] = prof_costs_bin_S
    all_profs_costs_T[prof_type] = prof_costs_bin_T
    all_profs_count_S[prof_type] = prof_count_bin_S
    all_profs_count_T[prof_type] = prof_count_bin_T

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
