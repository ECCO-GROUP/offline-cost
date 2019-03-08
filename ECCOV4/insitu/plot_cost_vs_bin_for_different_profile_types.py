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
import struct

sys.path.append('/home/owang/CODE/projects/Map_projection/')
import map_ploting_demo
import pyresample as pr

from  pyproj import Proj, transform
import matplotlib.path as mpath

import cartopy.crs as ccrs
import cartopy.feature as cfeature
 
#%%  HERE THE XX AND YY ARE IN THE PROJECTION CODE SPECIFIED BY the EPSG data_projection_code
def plot_pcolormesh_polar_stereographic(xx,yy, data, data_projection_code, \
                                   north_or_south, lat_lim, cmin, cmax, subplots=[1,1,1]):
    plt.subplots=subplots
    if north_or_south == 0: # south
        # 0 is south
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90, lat_lim], ccrs.PlateCarree())

    else: # north 
        # 1 is north
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, lat_lim, 90], ccrs.PlateCarree())


    #ax.add_feature(cfeature.LAND)
    #ax.add_feature(cfeature.OCEAN)

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)

    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=1, color='white', alpha=0.5, linestyle='--')

    ease_crs=ccrs.epsg(data_projection_code)

    p=plt.pcolormesh(xx, yy, data, transform=ease_crs, vmin=cmin, vmax=cmax,
                   cmap='jet')

    ax.add_feature(cfeature.LAND)
    plt.colorbar()
    ax.coastlines('110m', linewidth=0.8)
#    plt.show()
#%%
# global fields pass xx_ul, yy_ul not the full xx and yy with corners
def plot_pcolormesh_global(xx,yy, data, data_projection_code,
                      global_code, cmin, cmax, subplots=[1,1,1]):
    plt.subplots=subplots
    #  https://scitools.org.uk/cartopy/docs/latest/crs/projections.html
    if global_code == 0:
        ax = plt.axes(projection=ccrs.Robinson())
    elif global_code == 1:
        ax = plt.axes(projection=ccrs.EqualEarth())
    elif global_code == 2:
        ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
    elif global_code == 3:
        ax = plt.axes(projection=ccrs.PlateCarree()) # lat/lon


    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=1, color='white', alpha=0.5, linestyle='--')

    projection =ccrs.epsg(data_projection_code)

    plt.pcolormesh(xx, yy, data, transform=projection, vmin=cmin, vmax=cmax,
                   cmap='jet')
    plt.colorbar()
    ax.coastlines('110m', linewidth=0.8)
    ax.add_feature(cfeature.LAND)
#    plt.show()



#%%

#%%
def plot_field(xx,yy,data, cmin, cmax, n_lat_lim = -45):

    #fig = plt.figure(figsize=[7,7])
    ax = plt.axes(projection=ccrs.SouthPolarStereo())

    # Limit the map to -60 degrees latitude and below.
    ax.set_extent([-180, 180, -90, n_lat_lim], ccrs.PlateCarree())
    #ax.add_feature(cfeature.LAND)
    #ax.add_feature(cfeature.OCEAN)

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)

    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=1, color='white', alpha=0.5, linestyle='--')


    ease_crs=ccrs.epsg(6932)

    p=plt.pcolormesh(xx, yy, data, transform=ease_crs, vmin=cmin, vmax=cmax,
                   cmap='jet')

    ax.add_feature(cfeature.LAND)

    ax.coastlines('110m', linewidth=0.8)



#    ax.gridlines()
    plt.colorbar(p)#, format='%.0e')
    plt.show()

#%%
def plot_field_contourf(xx,yy,data, levels, cmin, cmax, n_lat_lim = -45, cmap='jet'):

    #fig = plt.figure(figsize=[7,7])
    ax = plt.axes(projection=ccrs.SouthPolarStereo())

    # Limit the map to -60 degrees latitude and below.
    ax.set_extent([-180, 180, -90, n_lat_lim], ccrs.PlateCarree())
    #ax.add_feature(cfeature.LAND)
    #ax.add_feature(cfeature.OCEAN)

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)

    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=1, color='black', alpha=0.5, linestyle='--')


    ease_crs=ccrs.epsg(6932)

    plt.contourf(xx, yy, data, levels, transform=ease_crs, vmin=cmin, vmax=cmax,
                   cmap=cmap)

    ax.add_feature(cfeature.LAND)

    ax.coastlines('110m', linewidth=0.8)



#    ax.gridlines()
    plt.colorbar()
    plt.show()

#%%
def plot_empty_so():

    ax = plt.axes(projection=ccrs.SouthPolarStereo())

    # Limit the map to -60 degrees latitude and below.
    ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)

    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=.5, color='gray', alpha=0.5, linestyle='--')

    ax.coastlines('110m', linewidth=0.8)

    plt.show()

#%%       
def plot_points_so(lons, lats, markersize):

    ax = plt.axes(projection=ccrs.SouthPolarStereo())

    # Limit the map to -60 degrees latitude and below.
    ax.set_extent([-180, 180, -90, -45], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)

    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=.5, color='gray', alpha=0.5, linestyle='--')

    ax.coastlines('110m', linewidth=0.8)

    plt.scatter(lons, lats, marker='.', color='red', s=markersize,
             alpha=0.7, transform=ccrs.Geodetic())

#    plt.show()

#%%
def plot_points_scatter_so(lons, lats, markersize, colors):

    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Limit the map to -60 degrees latitude and below.
    ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)

    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=.5, color='gray', alpha=0.5, linestyle='--')

    ax.coastlines('110m', linewidth=0.8)

    plt.scatter(lons, lats, c=colors, s=markersize, marker='.',
             alpha=0.7, transform=ccrs.Geodetic())

#    plt.show()

#%%
def define_grid(epsg_code):    
    #epsg_code=6932
    # Download masks from 
    # ftp://sidads.colorado.edu/pub/DATASETS/nsidc0609_loci_ease2/
    
    # see  https://nsidc.org/ease/clone-ease-grid-projection-gt
    
    ease_grid_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/JPL/input/mask/nsidc0609_loci_ease2/'
    # Load NSIDC EASE-2 Grid for the southern ocean, 25km grid resolution.
    #    https://nsidc.org/data/ease/ease_grid2.html
    # https: // daacdata.apps.nsidc.org / pub / DATASETS / nsidc0609_loci_ease2 / south /
    s25_Proj = Proj(init='epsg:'+str(epsg_code)) #
    latlonProj = Proj(init='epsg:4326') # lat lon
    proj_id = 'EPSG:'+str(epsg_code)
#%%
    if epsg_code == 6932:
    #https: // epsg.io / 6932
        area_id = 'WGS 84 / NSIDC EASE-Grid 2.0 South'
        area_name = 'WGS 84 / NSIDC EASE-Grid 2.0 South'
        proj4_str = '+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
        nx=720
        ny=720
        s25_ease_grid_fname = 'EASE2_S25km.LOCImask_land50_coast0km.'+str(nx)+'x'+str(ny)+'.bin'
    
        # the x and y limits in the NW corner [m]
        s25_ulxmap = -8987500
        s25_ulymap =  8987500
        # the spacing in x and y [m]
        s25_ease_dx = 25000
        s25_ease_dy = 25000
#%%
    elif epsg_code == 6931:
    #https: // epsg.io / 6931
        area_id = 'WGS 84 / NSIDC EASE-Grid 2.0 North'
        area_name = 'WGS 84 / NSIDC EASE-Grid 2.0 North'
        proj4_str = '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
        nx=720
        ny=720
        s25_ease_grid_fname = 'EASE2_N25km.LOCImask_land50_coast0km.'+str(nx)+'x'+str(ny)+'.bin'
        
        # the x and y limits in the NW corner [m]
        s25_ulxmap = -8987500
        s25_ulymap =  8987500
        # the spacing in x and y [m]
        s25_ease_dx = 25000
        s25_ease_dy = 25000
    #%% GLOBAL 25KM MAP

    elif epsg_code == 6933:
    #https: // epsg.io / 6933
        area_id = 'WGS 84 / NSIDC EASE-Grid 2.0 Global'
        area_name = 'WGS 84 / NSIDC EASE-Grid 2.0 Global'
        proj4_str = '+proj=cea +lat 0=0  lon 0=0 +lat 1=30 +x 0=0 +y 0=0 +ellps=WGS84 +datum=WGS84 +units=m'
        nx=1388
        ny=584
        s25_ease_grid_fname = 'EASE2_M25km.LOCImask_land50_coast0km.'+str(nx)+'x'+str(ny)+'.bin'
    
        # the x and y limits in the NW corner [m]
        s25_ulxmap = -1.7367530439999998E7
        s25_ulymap =  7307375.92
        # the spacing in x and y [m]
        s25_ease_dx = 25025.26
        s25_ease_dy = 25025.26
    else:
        print("incorrect epsg_code")
        sys.exit()
    #%%
    nxp1=nx+1
    nyp1=ny+1
    ease_grid_full_path = ease_grid_dir + '/' + s25_ease_grid_fname
    
    # load the EASE-2 Southern Ocean Mask
    with open(ease_grid_full_path, 'rb') as f:
        s25_LOCImask= np.fromfile(f, dtype=np.uint8)
    
    s25_LOCImask = np.reshape(s25_LOCImask, (nx, ny))
    
    # the x and y of the upper left corners of the grid cells
    s25_x_ul = np.linspace(s25_ulxmap, -s25_ulxmap, nx)
    s25_y_ul = np.linspace(s25_ulymap, -s25_ulymap, ny)
    s25_xx_ul, s25_yy_ul = np.meshgrid(s25_x_ul,s25_y_ul)
    
    # the x and y of the center of the grid cell
    s25_x = s25_x_ul + s25_ease_dx/2.0
    s25_y = s25_y_ul + s25_ease_dy/2.0
    s25_xx,s25_yy = np.meshgrid(s25_x,s25_y)
    
    # the lons and lats of the grid cell centers 
    s25_lon, s25_lat = transform(s25_Proj, latlonProj, s25_xx.ravel(), s25_yy.ravel())
    s25_lon = np.reshape(s25_lon, s25_xx.shape)
    s25_lat = np.reshape(s25_lat, s25_yy.shape)
    
    # the x and y of the lower right corners of the grid cells
    s25_x_lr = s25_x_ul + s25_ease_dx
    s25_y_lr = s25_y_ul + s25_ease_dy
    s25_xx_lr, s25_yy_lr = np.meshgrid(s25_x_lr, s25_y_lr)
    
    s25_x_corners = np.linspace(s25_ulxmap, -s25_ulxmap+s25_ease_dx, nxp1)
    s25_y_corners = np.linspace(s25_ulymap, -s25_ulymap+s25_ease_dx, nyp1)
    s25_xx_corners, s25_yy_corners = np.meshgrid(s25_x_corners,s25_y_corners)
    
    return s25_lon, s25_lat, s25_xx_corners, s25_yy_corners, s25_xx_ul, s25_yy_ul
        

# dictionary that will hold the different prof_stats 
prof_stats = {}
nbins=10242

# legend list
lgnd = []
run_codes = [2]
diravgbin = '/mnt/intraid/ecco-intrnl/data07/owang/Ecco_data/v4/JPL/r4/release2/input/data/profiles/ECCO_v4r2/geodesic_bins/'
fbinlat = diravgbin + 'sphere_point_n_10242_lats.bin'
fbinlon = diravgbin + 'sphere_point_n_10242_lons.bin'
with open(fbinlat, mode='rb') as file: # b is important -> binary
    binlattmp = file.read()
    file.seek(0)
    binlat=np.fromfile(file,dtype='>f',count=nbins)
    #binlat=struct.unpack('d',binlattmp[0:8])
with open(fbinlon, mode='rb') as file: # b is important -> binary
    binlontmp = file.read()
    file.seek(0)
    binlon=np.fromfile(file,dtype='>f',count=nbins)
orig_lats_1d=binlat
orig_lons_1d=binlon

radius_of_influence=150000
nprocs_user=1

#%%
#%% grid
epsg_code_north=6931
s25_lon_north, s25_lat_north, s25_xx_corners_north, s25_yy_corners_north, s25_xx_ul_north, s25_yy_ul_north = define_grid(epsg_code_north)

epsg_code_south=6932
s25_lon_south, s25_lat_south, s25_xx_corners_south, s25_yy_corners_south, s25_xx_ul_south, s25_yy_ul_south = define_grid(epsg_code_south)

epsg_code_glb=6933
s25_lon_glb, s25_lat_glb, s25_xx_corners_glb, s25_yy_corners_glb,  s25_xx_ul_glb, s25_yy_ul_glb= define_grid(epsg_code_glb)



#%%
#%%

orig_grid = pr.geometry.SwathDefinition(lons=orig_lons_1d,
                                        lats=orig_lats_1d)

# =============================================================================
# # loop through the different run codes, load the prof_stats for each
for run_code in run_codes:
    if run_code == 0:
        # ECCO v4 R3 original
        output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/h8/h8i_i48/python_output/prof_cost_comparisons/'
        output_name = 'ECCOv4r3_original_bin'
    
    elif run_code == 1 or run_code == 2:
        if run_code == 1:
            # iter 115
            output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/r042/r042f/run.v4_rls2.042ff10.iter115/python_output/prof_cost_comparisons/'
            output_name = 'ECCOv4r3_extension_it115'
        elif run_code == 2:
        # iter 118
            output_dir = '/mnt/intraid/ecco-intrnl/data10/owang/Ecco_data/v4/Processed/r4/release3/r042/r042f/run.v4_rls2.042ff10.iter118/python_output/prof_cost_comparisons/'
            output_name = 'ECCOv4r3_extension_it118'
    else:
        print('not a valid run code ', run_code)
        sys.exit()
        
   
    os.chdir(output_dir)
    prof_stats[run_code] = np.load(output_name + '.npy').item()
    
    lgnd.append(prof_stats[run_code]['run_name'])
        
    #%%
min_year = 1992
max_year = 2018

# 
# # a counter of the total number of figures that we've made     
total_figures = 0
# 
# # a dictionary matching profile type with figure numbers for T and S
T_fig_nums = {}
S_fig_nums = {}
# 
# 
# # loop through the different prof_stats objects
costs_keys = []
# 
for run_code in run_codes:
    # pull the keys out for the different profile types for T and S
    costs_keys = (set(prof_stats[run_code]['costs_T'].keys()) | set(costs_keys))
    costs_keys = (set(prof_stats[run_code]['costs_S'].keys()) | set(costs_keys))
    
#print costs_keys
#%%    
# start the plotting by closing all existing plots
plt.close('all')
# 
#for run_code in run_codes:
#         
#     # loop through each profile type that has a T cost
#     
for prof_type in costs_keys:
        
    for run_code in run_codes:
        print ('run_code PROFILE TYPE: ', run_code, prof_type)
#         
        if prof_type in prof_stats[run_code]['costs_T'].keys():
            cost_T  = prof_stats[run_code]['costs_T'][prof_type]
            count_T = prof_stats[run_code]['count_T'][prof_type]
            #print ('count_T: ', count_T)
            #if count_T != 0:
            cpd_T   = cost_T / count_T

        else:
            cost_T = np.zeros(nbins)
            count_T = np.zeros(nbins)
            cpd_T = np.zeros(nbins)
        orig_field=cpd_T
# 
#         # see if we have a figure associated with this profile type
#         # already (we probably will if we have more than one run code)
#         # if so, get the figure number
#         # if not, make a new figure number
        if prof_type in T_fig_nums:
            fig_num = T_fig_nums[prof_type]
        else:
            total_figures += 1
            fig_num = total_figures
            # add this figure number to the dictionary
            T_fig_nums[prof_type] = fig_num
                #print "============="     
        fig = plt.figure(num=fig_num, figsize=[8,8]);fig.clf()
        cmin=0
        cmax=5
        if prof_type=='itp':
#        tyty1=data_latlon_projection
            cmin=0
            cmax=20
            north_or_south = 1;
            lat_limit = 65    
            new_grid  = pr.geometry.GridDefinition(lons=s25_lon_north,
                                           lats=s25_lat_north)
            data_latlon_projection = \
                pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
                                            radius_of_influence=radius_of_influence, \
                                            fill_value=None, nprocs=nprocs_user)            
            
            tyty1=np.copy(data_latlon_projection)
            tyty1=np.reshape(tyty1,s25_lon_north.shape)
            plot_pcolormesh_polar_stereographic(s25_xx_corners_north,s25_yy_corners_north, tyty1,
                                        epsg_code_north, north_or_south, lat_limit, cmin, cmax)
            plt.show()
        else:
 #           plt.subplot(2,1,1)
# =============================================================================
#             new_grid  = pr.geometry.GridDefinition(lons=s25_lon_glb,
#                                            lats=s25_lat_glb)
#             data_latlon_projection = \
#                 pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
#                                             radius_of_influence=radius_of_influence, \
#                                             fill_value=None, nprocs=nprocs_user)            
#             
#             tyty1=np.copy(data_latlon_projection)
#             tyty1=np.reshape(tyty1,s25_lon_glb.shape)
#             global_code=1
#             plot_pcolormesh_global(s25_xx_ul_glb,s25_yy_ul_glb, tyty1, epsg_code_glb,
#                       global_code, cmin, cmax)
# =============================================================================
# =============================================================================
#%%
#            plt.subplot(2,2,3)
            new_grid  = pr.geometry.GridDefinition(lons=s25_lon_north,
                                           lats=s25_lat_north)
            data_latlon_projection = \
                pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
                                            radius_of_influence=radius_of_influence, \
                                            fill_value=None, nprocs=nprocs_user)            
            
            tyty1=np.copy(data_latlon_projection)
            tyty1=np.reshape(tyty1,s25_lon_north.shape)            
            north_or_south = 1;
            lat_limit = 45    
            plot_pcolormesh_polar_stereographic(s25_xx_corners_north,s25_yy_corners_north, tyty1,
                                        epsg_code_north, north_or_south, lat_limit, cmin, cmax)

#%%
# =============================================================================
# #            plt.subplot(2,2,4)
#             new_grid  = pr.geometry.GridDefinition(lons=s25_lon_south,
#                                            lats=s25_lat_south)
#             data_latlon_projection = \
#                 pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
#                                             radius_of_influence=radius_of_influence, \
#                                             fill_value=None, nprocs=nprocs_user)            
#             
#             tyty1=np.copy(data_latlon_projection)
#             tyty1=np.reshape(tyty1,s25_lon_south.shape)                        
#             north_or_south = 0;
#             lat_limit = -45    
#             plot_pcolormesh_polar_stereographic(s25_xx_corners_south,s25_yy_corners_south, tyty1,
#                                         epsg_code_south, north_or_south, lat_limit, cmin, cmax)  
# #%%
# =============================================================================
            plt.show()
        #plt.title('Temperature : ' + prof_type + str(run_code))
        subnm='_v4r3'
        if run_code == 1:
            subnm='_v4ext'
        plt.title('Temperature : ' + prof_type + subnm)
#        py.savefig(output_dir + '/fig/vsbin/' + 'T_'+prof_type+subnm+'.png')
#        py.savefig(output_dir + '/fig/vsbin/' + 'T_'+prof_type+subnm+'_Antarctic.png')
        py.savefig(output_dir + '/fig/vsbin/' + 'T_'+prof_type+subnm+'_Arctic.png')

        print ('run_code PROFILE TYPE: ', run_code, prof_type)
#         
        if prof_type in prof_stats[run_code]['costs_S'].keys():
            cost_S  = prof_stats[run_code]['costs_S'][prof_type]
            count_S = prof_stats[run_code]['count_S'][prof_type]
            #print ('count_T: ', count_T)
            #if count_T != 0:
            cpd_S   = cost_S / count_S

        else:
            cost_S = np.zeros(nbins)
            count_S = np.zeros(nbins)
            cpd_S = np.zeros(nbins)
#         
#         #%%
#         #cost_T=cost_T.data
#         #count_T=count_T.data
#         #cpd_T=cpd_T.data
#         
        orig_field=cpd_S

    
# 
#%%
#        print  run_code
#         #orig_field = orig_field.values
# 
#         # see if we have a figure associated with this profile type
#         # already (we probably will if we have more than one run code)
#         # if so, get the figure number
#         # if not, make a new figure number
        if prof_type in S_fig_nums:
            fig_num = S_fig_nums[prof_type]
        else:
            total_figures += 1
            fig_num = total_figures
            # add this figure number to the dictionary
            S_fig_nums[prof_type] = fig_num
                #print "============="     
        fig = plt.figure(num=fig_num, figsize=[8,8]);fig.clf()
        cmin=0
        cmax=5
        if prof_type=='itp':
#        tyty1=data_latlon_projection
            cmin=0
            cmax=20
            north_or_south = 1;
            lat_limit = 65    
            new_grid  = pr.geometry.GridDefinition(lons=s25_lon_north,
                                           lats=s25_lat_north)
            data_latlon_projection = \
                pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
                                            radius_of_influence=radius_of_influence, \
                                            fill_value=None, nprocs=nprocs_user)            
            
            tyty1=np.copy(data_latlon_projection)
            tyty1=np.reshape(tyty1,s25_lon_north.shape)
            plot_pcolormesh_polar_stereographic(s25_xx_corners_north,s25_yy_corners_north, tyty1,
                                        epsg_code_north, north_or_south, lat_limit, cmin, cmax)
            plt.show()
        else:
 #           plt.subplot(2,1,1)
# =============================================================================
#             new_grid  = pr.geometry.GridDefinition(lons=s25_lon_glb,
#                                            lats=s25_lat_glb)
#             data_latlon_projection = \
#                 pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
#                                             radius_of_influence=radius_of_influence, \
#                                             fill_value=None, nprocs=nprocs_user)            
#             
#             tyty1=np.copy(data_latlon_projection)
#             tyty1=np.reshape(tyty1,s25_lon_glb.shape)
#             global_code=1
#             plot_pcolormesh_global(s25_xx_ul_glb,s25_yy_ul_glb, tyty1, epsg_code_glb,
#                       global_code, cmin, cmax)
# =============================================================================
#%%
#            plt.subplot(2,2,3)
            new_grid  = pr.geometry.GridDefinition(lons=s25_lon_north,
                                           lats=s25_lat_north)
            data_latlon_projection = \
                pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
                                            radius_of_influence=radius_of_influence, \
                                            fill_value=None, nprocs=nprocs_user)            
            
            tyty1=np.copy(data_latlon_projection)
            tyty1=np.reshape(tyty1,s25_lon_north.shape)            
            north_or_south = 1;
            lat_limit = 45    
            plot_pcolormesh_polar_stereographic(s25_xx_corners_north,s25_yy_corners_north, tyty1,
                                        epsg_code_north, north_or_south, lat_limit, cmin, cmax)

#%%
#            plt.subplot(2,2,4)
# =============================================================================
#             new_grid  = pr.geometry.GridDefinition(lons=s25_lon_south,
#                                            lats=s25_lat_south)
#             data_latlon_projection = \
#                 pr.kd_tree.resample_nearest(orig_grid, orig_field, new_grid, \
#                                             radius_of_influence=radius_of_influence, \
#                                             fill_value=None, nprocs=nprocs_user)            
#             
#             tyty1=np.copy(data_latlon_projection)
#             tyty1=np.reshape(tyty1,s25_lon_south.shape)                        
#             north_or_south = 0;
#             lat_limit = -45    
#             plot_pcolormesh_polar_stereographic(s25_xx_corners_south,s25_yy_corners_south, tyty1,
#                                         epsg_code_south, north_or_south, lat_limit, cmin, cmax)
# 
# =============================================================================
            plt.show()
        #plt.title('Temperature : ' + prof_type + str(run_code))
        subnm='_v4r3'
        if run_code == 1:
            subnm='_v4ext_it03'
        elif run_code == 2:
            subnm='_v4ext_it06'            
        plt.title('Salinity : ' + prof_type + subnm)
#        py.savefig(output_dir + '/fig/vsbin/' + 'S_'+prof_type+subnm+'png')        
        #py.savefig(output_dir + '/fig/vsbin/' + 'S_'+prof_type+subnm+'_Antarctic.png')
        py.savefig(output_dir + '/fig/vsbin/' + 'S_'+prof_type+subnm+'_Arctic.png')
            
# =============================================================================
#
