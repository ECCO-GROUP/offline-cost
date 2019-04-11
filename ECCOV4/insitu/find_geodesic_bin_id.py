#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 13:50:21 2019

@author: ifenty
"""


# %%
from __future__ import division

import numpy as np
import sys
import os
import glob
import pylab as py
import pyresample as pr
import csv
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from mpl_toolkits.mplot3d import Axes3D
#%%

def find_geo_bin_id(geo_bin_fname,
                    geo_bin_dir, 
                    sample_lon_lat):
    
    #%%
    csv_file = csv.reader(open(bin_dir + '/' + fname, 'r'))
    
    lons = []
    lats = []
    
    grid_lon_lat = np.empty((10242,2))
    k=0
    for row in csv_file:
        # lon first
        grid_lon_lat[k,0] = np.float(row[0])
        # lat second
        grid_lon_lat[k,1] = np.float(row[1])
        k+=1

    geo_bin_ids = np.arange(1,10243)

    clat = np.cos(grid_lon_lat[:,1]* np.pi / 180.)
    clon = np.cos(grid_lon_lat[:,0]* np.pi / 180.)
    slat = np.sin(grid_lon_lat[:,1]* np.pi / 180.)
    slon = np.sin(grid_lon_lat[:,0]* np.pi / 180.)
    
    xg = clat * clon
    yg = clat * slon
    zg = slat
    geo_xx=np.column_stack((xg,yg,zg))
    
    clat = np.cos(sample_lon_lat[:,1]* np.pi / 180.)
    clon = np.cos(sample_lon_lat[:,0]* np.pi / 180.)
    slat = np.sin(sample_lon_lat[:,1]* np.pi / 180.)
    slon = np.sin(sample_lon_lat[:,0]* np.pi / 180.)
    
    xs = clat * clon
    ys = clat * slon
    zs = slat
    sample_xx=np.column_stack((xs,ys,zs))
    
    fig = plt.figure()
    ax = fig.add_subplot(111,projection= '3d')
    ax.scatter(xg,yg,zs=zg,s=10,c='k')
    ax.scatter(xs,ys,zs=zs,s=10,c='r')
    #%%
    geo_bin_ids_at_sample_locations = \
        interp.griddata(geo_xx, 
                        geo_bin_ids, 
                        sample_xx, 
                        method='nearest')
    print geo_bin_ids_at_sample_locations  
    #%%
    return geo_bin_ids_at_sample_locations

   # In[]
if __name__ == '__main__':
    
    #%%
    geo_bin_dir ='/home/ifenty/git_repo_mine/offline-cost/ECCOV4/insitu'
    geo_bin_fname = '10242_bin_locations.csv'
    
    sample_lon_lat = np.array([[-40, 40],[-20, 20]])
    find_geo_bin_id(geo_bin_fname, geo_bin_dir, sample_lon_lat)
    
