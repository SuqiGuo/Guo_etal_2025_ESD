#! /usr/bin/env python

# This script prepares ESM output with LCLM chessboard pattern 
# to separate into local, nonlocal, and total effects 

# It is based on the scripts of Johannes Winckler for his 2017 study
# doi: 10.1175/JCLI-D-16-0067.1


# This function takes a difference map and a chessboard pattern as input.
# Calculate everything on extended grid, such that values around zero meridian can be interpolated correctly.
# Output are interpolated local, non-local, and total effects.


import numpy as np
#from scipy.interpolate import RectBivariateSpline
import scipy.interpolate as interpolate
import datetime as dt
import os
from copy import deepcopy as cp
import netCDF4 as nc
import sys
from dask import delayed

# GET INFO FROM INPUT
if os.path.isfile(sys.argv[1]):
    FILE_DIFFERENCEMAP = sys.argv[1]
    VAR = sys.argv[2]
    
elif os.path.isfile(sys.argv[2]):
    
    FILE_DIFFERENCEMAP = sys.argv[2]
    VAR = sys.argv[1]
    
else:    
    print(sys.argv[:], " there is no file in the parsed arguments...")
    sys.exit( 2)

# GET INFO ABOUT THE MODEL
if "cesm" in FILE_DIFFERENCEMAP:
    MODEL = "cesm"
    FILL_VALUE = 1e36
elif "mpiesm" in FILE_DIFFERENCEMAP:
    MODEL = "mpiesm"
    FILL_VALUE = 1e20
elif "ecearth" in FILE_DIFFERENCEMAP:
    MODEL = "ecearth"    
else:
    print("The model could not be identified... :|")
    sys.exit( 2)
    

#INPUT_WD = "/home/felixh/data/mpiesm/chessboard/"
INPUT_WD = "/home/b/b380948/lamaclima/chessboard_tests/"
DATA_WD = "/home/felixh/data/mpiesm/differencemaps/ctl-crop/"

# IMPORT GLOBAL CHESSBOARD DISTRIBUTION ON LAND
#ifile = INPUT_WD + "lamaclima_experiments_chessboard_pattern_no-missing.nc"
ifile = INPUT_WD + "lamaclima_experiments_chessboard_pattern_" + MODEL + ".nc"
f = nc.Dataset(ifile, 'r')
CHESSBOARD = f.variables['chessboard_pattern'][:]
LAT = f.variables['lat'][:]
print(LAT)
LON = f.variables['lon'][:]
f.close()

# IMPORT GLOBAL CHESSBOARD DISTRIBUTION GLOBALLY
ifile = INPUT_WD + "lamaclima_experiments_global_chessboard_pattern_" + MODEL + ".nc"
f = nc.Dataset(ifile, 'r')
GLOBAL_CHESSBOARD = f.variables['chessboard_pattern'][:]
f.close()

# IMPORT SEA LAND MASK
if MODEL == "mpiesm":
    ifile = INPUT_WD + "jsbach_T63GR15_11tiles_5layers_2005_dynveg_slm_glac.nc"
    f = nc.Dataset(ifile, 'r')
    SEA_LAND_MASK = f.variables['slm'][:]
    GLACIER = f.variables['glac'][:]    

elif MODEL == "cesm":
    ifile = INPUT_WD + "landmask_glacier_cesm.nc"
    f = nc.Dataset(ifile, 'r')
    #SEA_LAND_MASK = f.variables['LANDFRAC_PFT'][:]
    SEA_LAND_MASK = f.variables['landmask'][:] * 1.0
    GLACIER = f.variables['GLACIER_REGION'][:] * 1.0 # simple transformation to float...
    
f.close()    

# IMPORT DIFFERENCEMAP
#ifile = DATA_WD + "nep_Emon_MPI-ESM1-2-LR_215501-217412_ctl-crop_timavg.nc"
ifile = FILE_DIFFERENCEMAP
f = nc.Dataset(ifile, 'r+')
# CHECK IF INVERTED LAT VALUES IN ORIGINAL FILE AND KEEP THE ORIGINAL CONVENTION
LAT_DIFF = f.variables["lat"][:]
print(LAT_DIFF)
LON_DIFF = f.variables["lon"][:]
if np.around( LAT_DIFF[0], 5) == np.around( LAT[0], 5):
    i = 0
    j = 1    

elif np.around( LAT_DIFF[0], 5) == np.around( -1 * LAT[0], 5):
    i = -1
    j = -1
    
else:
    print("Latitude information of used masks doesn't match... :|")
    sys.exit( 2)
print(i)
print(j)
#UNIT = "m/s"
UNIT = "/"
#UNIT = f.variables[VAR].units
NLAT = len(LAT_DIFF)
NLON = len(LON_DIFF)
#DIFFERENCEMAP_ALLTIMESTEPS = f.variables[VAR][:,-1::-1,:] # time, lat, lon
DIFFERENCEMAP_ALLTIMESTEPS = f.variables[VAR][:,i::j,:] # time, lat, lon
NSTEPS = len( DIFFERENCEMAP_ALLTIMESTEPS[:,0,0])
#f.close()


# ----- CREATE CONSTANT INPUT MASKS

# 4 % EXTENDED NUMBER OF GRIDPOINTS IN LON DIRECTION NEEDED FOR INTERPOLATION AROUND 0 DEGREE W
# WILL BE RESTRICTED TO ORIGINAL GRID LATER.

# FOR MPI-ESM: INSTEAD OF 192 LONGITUDE ROWS WE ADD 8 LONGITUDE ROWS AND USE 200 LONGITUDE ROWS
#EXTEND_LON_NUMBER = int( np.around( len(LON) * 1.04, 0))
EXTEND_LON_NUMBER = len(LON) + 8

## CREATE EXTENDED CHESSBOARD MAP ON ALL LAND AND GLOBAL GRID POINTS -- non masked arrays
#CHESSBOARD_EXTENDED = np.tile( CHESSBOARD,2)[:,0:EXTEND_LON_NUMBER]
#GLOBAL_CHESSBOARD_EXTENDED = np.tile( GLOBAL_CHESSBOARD,2)[:,0:EXTEND_LON_NUMBER]

# CREATE EXTENDED CHESSBOARD MAP ON ALL LAND AND GLOBAL GRID POINTS -- masked arrays
CHESSBOARD_EXTENDED = np.ma.masked_equal( np.tile( CHESSBOARD,2)[:,0:EXTEND_LON_NUMBER], value=FILL_VALUE)
GLOBAL_CHESSBOARD_EXTENDED = np.ma.masked_equal( np.tile( GLOBAL_CHESSBOARD,2)[:,0:EXTEND_LON_NUMBER], value=FILL_VALUE)

# REVERSE GLOBAL CHESSBOARD GRID POINTS (1-->0; 0-->1)
GLOBAL_CHESSBOARD_EXTENDED_REVERSE = np.abs( GLOBAL_CHESSBOARD_EXTENDED - 1)

## CREATE EXTENDED SEA LAND MASK AND GLACIER MASK -- non masked arrays
#LAND_EXTENDED = np.tile( SEA_LAND_MASK,2)[:,0:EXTEND_LON_NUMBER]
#GLACIER_EXTENDED = np.tile( GLACIER,2)[:,0:EXTEND_LON_NUMBER]

# CREATE EXTENDED SEA LAND MASK AND GLACIER MASK -- masked arrays
LAND_EXTENDED = np.ma.masked_equal( np.tile( SEA_LAND_MASK,2)[:,0:EXTEND_LON_NUMBER], value=FILL_VALUE)
GLACIER_EXTENDED = np.ma.masked_equal( np.tile( GLACIER, 2)[:,0:EXTEND_LON_NUMBER], value=FILL_VALUE)    

def interpolate_nonlocal( chessboard_extended, differencemap_extended, land_extended, coastboxes):
    
    #---------coordinate of the grid boxes that are interpolated
    points = np.where( np.all( np.array(( chessboard_extended == 0, land_extended == 1)), axis = 0))
    values = differencemap_extended[np.all( np.array(( chessboard_extended == 0, land_extended == 1)), axis = 0)]
    #values = differencemap_extended[points]
    
    #grid_x, grid_y = np.mgrid[0:96:1, 0:land_extended.shape[1]:1]
    grid_x, grid_y = np.mgrid[0:land_extended.shape[0]:1, 0:land_extended.shape[1]:1]
    
    interped = interpolate.griddata( np.array(( points[1] * 1.0, points[0] * 1.0)).T, values,\
        (grid_y,grid_x), method='linear')
    interped_nearest = interpolate.griddata( np.array(( points[1] * 1.0, points[0] * 1.0)).T, values,\
        (grid_y,grid_x), method='nearest')
    
    non_local = cp( interped)
    non_local[coastboxes == 1] = interped_nearest[coastboxes == 1]
    non_local[land_extended == 0] = differencemap_extended[land_extended == 0]

    return non_local


def interpolate_local( chessboard_extended, local_extended, land_extended, glacier_extended, coastboxes):
    
    #---------coordinate of the grid boxes that are interpolated
    points = np.where( chessboard_extended == 1)
    values = local_extended[ chessboard_extended == 1]
    
    #grid_x, grid_y = np.mgrid[0:len(LAT):1, 0:land_extended.shape[1]:1]
    grid_x, grid_y = np.mgrid[0:land_extended.shape[0]:1, 0:land_extended.shape[1]:1]
    
    interped = interpolate.griddata( np.array(( points[1] * 1.0, points[0] * 1.0)).T, values,\
        (grid_y,grid_x), method='linear')
    interped_nearest = interpolate.griddata( np.array(( points[1] * 1.0, points[0] * 1.0)).T, values,\
        (grid_y,grid_x), method='nearest')
        
    local_interp = cp( interped)
    local_interp[coastboxes == 1] = interped_nearest[coastboxes == 1]
    
    # MASK OUT WATER BODIES AND GLACIER TOGETHER
    local_interp[land_extended == 0] = 0
    local_interp[glacier_extended == 1] = 0

    return local_interp

def subtract(extended_diff_map, remote_extended):
    return extended_diff_map - remote_extended

    
def determine_coast_boxes(chessboard_extended, global_chessboard_extended, global_chessboard_extended_reverse):
    # coastpixels are boxes that are ocean or that are near the ocean,
    # such that an ocean grid box would be used for interpolating the lcc grid boxes
    coastboxes = np.zeros( chessboard_extended.shape)
    #coastboxes_2 = np.zeros( chessboard_extended.shape)
    
    for latindex in range( chessboard_extended.shape[0]):
        for lonindex in range(chessboard_extended.shape[1]):
            # in a 9 boxes window around the grid box, look if there are more global chessboard pixels
            # than land chessboard pixels
            # Johannes took global_chessboard_extended_reverse for his calculation....
            if (np.sum( global_chessboard_extended[latindex-2:latindex+3,lonindex-2:lonindex+3]) -
                np.sum( chessboard_extended[latindex-2:latindex+3,lonindex-2:lonindex+3]) > 0):
                coastboxes[latindex,lonindex]=1
            #if (np.sum( global_chessboard_extended_reverse[latindex-2:latindex+3,lonindex-2:lonindex+3]) -
                #np.sum( chessboard_extended[latindex-2:latindex+3,lonindex-2:lonindex+3]) > 0):
                #coastboxes_2[latindex,lonindex]=1       
    
    return coastboxes#, coastboxes_2


def write_netcdf( non_local_timeseries,local_timeseries, f, UNIT, LAT_DIFF, NSTEPS, NLAT, NLON):
    # USE COPY OF FILE TO CREATE NEW NETCDF FILE
    print("Writing netCDF4 file...")
   # ifile = FILE_DIFFERENCEMAP   #.replace(".nc","") + "_signal-separated.nc"
   # f = nc.Dataset( ifile, 'r+') #, format='NETCDF4')

    # Define global attributes
    f.creation_date=str(dt.datetime.today())
    f.contact='felix.havermann@lmu.de; Ludwig-Maximilians-University Munich'
    f.comment='Produced with Script ' + os.path.basename(__file__)
    #f.title='Land cover fractions of MPIESM-1.2-LR LAMACLIMA ' + SIMULATION_TYPE + ' simulation'
    #f.subtitle=''
    #f.anything=''
    #f_restart.comment='Only the variable cover_fract_pot was changed to a 100 % forest world and all other cover types are set to zero.\n'\
    #f.comment='Only the variable cover_fract_pot was changed to a 100 % forest world and all other cover types are set to fract_small (1e-10).\n'\
        #'This file is used for a 100 % forest simulation within the LAMACLIMA project\n'\
        #'Produced with Script ' + os.path.basename(__file__)
    #f.createVariable(LCT_CONVERSION[key][1], np.float32,('lat','lon'), fill_value=FILL_VALUE)

    # NON-LOCAL SIGNAL 
    f.createVariable( VAR + "_nonlocal", np.float64,('time','lat','lon'), fill_value=FILL_VALUE)
    #f.variables[VAR + "_nonlocal"][0,:,:] = non_local_timeseries[:,i::j,:]
    #f.variables[VAR + "_nonlocal"][1:,:,:] = np.zeros(( NSTEPS-1, NLAT, NLON))
    f.variables[VAR + "_nonlocal"][:,:,:] = non_local_timeseries[:,i::j,:]
    f.variables[VAR + "_nonlocal"].longname = 'Non-local (=remote) effect of LCLM change'
    f.variables[VAR + "_nonlocal"].units = UNIT
    f.variables[VAR + "_nonlocal"].grid_type = 'gaussian'

    # TOTAL SIGNAL 
    total = local_timeseries[:,i::j,:] + non_local_timeseries[:,i::j,:]
    f.createVariable( VAR + "_total", np.float64,('time','lat','lon'), fill_value=FILL_VALUE)
    #f.variables[VAR + "_total"][0,:,:] = local_interpolated[:,i::j,:] + non_local_interpolated[:,i::j,:]
    f.variables[VAR + "_total"][:,:,:] = total
    #f.variables[VAR + "_total"][1:,:,:] = np.zeros(( NSTEPS-1, NLAT, NLON))
    f.variables[VAR + "_total"].longname = 'Total (=local + non-local) effect of LCLM change'
    f.variables[VAR + "_total"].units = UNIT
    f.variables[VAR + "_total"].grid_type = 'gaussian'

    # LOCAL SIGNAL 
    MASK = np.tile( CHESSBOARD[i::j,:], (NSTEPS,1,1))
    f.createVariable( VAR + "_local", np.float64,('time','lat','lon'), fill_value=FILL_VALUE)
    #f.variables[VAR + "_local"][0,:,:] = local_interped[:,i::j,:]
    #f.variables[VAR + "_local"][0,:,:] = local_timeseries[i::j,:]
    #f.variables[VAR + "_local"][1:,:,:] = np.zeros(( NSTEPS-1, NLAT, NLON))
    f.variables[VAR + "_local"][:,:,:] = np.ma.masked_array( local_timeseries[:,i::j,:], MASK.mask)
    f.variables[VAR + "_local"].longname = 'Local effect of LCLM change'
    f.variables[VAR + "_local"].units = UNIT
    f.variables[VAR + "_local"].grid_type = 'gaussian'

    f.close()
    print("writing finished.")


def get_original_grid_nonlocal(grid_extended, LON, step, non_local_timeseries):
    tmp1 = grid_extended[:, 4:len(LON)]    
    tmp2 = grid_extended[:, len(LON):len(LON)+4]
    non_local_timeseries[step,:,:] = np.hstack((tmp2, tmp1))
    
    return non_local_timeseries

def get_original_grid_local(grid_extended, LON, step, local_timeseries):
    tmp1 = grid_extended[:, 4:len(LON)]    
    tmp2 = grid_extended[:, len(LON):len(LON)+4]
    local_timeseries[step,:,:] = np.hstack((tmp2, tmp1))
    
    return local_timeseries

# DETERMINE COASTBOXES TO APPLY NEAREST NEIGHBOUR CALCULATION FOR THOSE BOXES      
coastboxes = determine_coast_boxes( CHESSBOARD_EXTENDED, GLOBAL_CHESSBOARD_EXTENDED, GLOBAL_CHESSBOARD_EXTENDED_REVERSE)

non_local_timeseries = np.ma.masked_all((NSTEPS, NLAT, NLON))
local_timeseries = np.ma.masked_all((NSTEPS, NLAT, NLON))
for step in range( NSTEPS):
    #print("Step", step, "of")
        
    differencemap = DIFFERENCEMAP_ALLTIMESTEPS[step,:,:]
    
    ## CREATE EXTENDED DIFFERENCE MAPS FROM CTL, CROP, IRR, FRST, AND HARV EXPERIMENTS -- non masked arrays
    #DIFFERENCEMAP_EXTENDED = np.tile( differencemap,2)[:,0:EXTEND_LON_NUMBER]

    # CREATE EXTENDED DIFFERENCE MAPS FROM CTL, CROP, IRR, FRST, AND HARV EXPERIMENTS -- masked arrays
    DIFFERENCEMAP_EXTENDED = np.ma.masked_equal( np.tile( differencemap, 2)[:,0:EXTEND_LON_NUMBER], value=FILL_VALUE)

    # CALCULATE NON-LOCAL EFFECTS
    #nonlocal_extended = interpolate_nonlocal( CHESSBOARD_EXTENDED, DIFFERENCEMAP_EXTENDED, LAND_EXTENDED, coastboxes)
    nonlocal_extended = delayed(interpolate_nonlocal)( CHESSBOARD_EXTENDED, DIFFERENCEMAP_EXTENDED, LAND_EXTENDED, coastboxes)

    # CALCULATE LOCAL EFFECTS
    #local_extended = DIFFERENCEMAP_EXTENDED - nonlocal_extended
    local_extended = delayed(subtract)( DIFFERENCEMAP_EXTENDED, nonlocal_extended)

    # INTERPOLATE LOCAL EFFECTS TO NO-LCC GRID BOXES    
    #local_interped_extended = interpolate_local( CHESSBOARD_EXTENDED, local_extended, LAND_EXTENDED, GLACIER_EXTENDED, coastboxes)
    local_interped_extended = delayed(interpolate_local)( CHESSBOARD_EXTENDED, local_extended, LAND_EXTENDED, GLACIER_EXTENDED, coastboxes)

    #non_local[:, 0:4] = nonlocal_extended[:, 192:196]

    # RESTRICT TO ORIGINAL GRID 
    #non_local=nonlocal_extended[:, 0:len(LON)]    
    #non_local[:, 0:4] = nonlocal_extended[:, len(LON):len(LON)+4]    
    ##non_local_masked = np.ma.masked_equal(non_local, value=FILL_VALUE)    
    #non_local_timeseries[step,:,:] = non_local
    non_local_timeseries = delayed(get_original_grid_nonlocal)(nonlocal_extended, LON, step, non_local_timeseries)
   

    #local = local_extended[:,0:len(LON)]
    #local[:, 0:4] = local_extended[:, len(LON):len(LON)+4]

    #local_interped = local_interped_extended[:,0:len(LON)]
    #local_interped[:, 0:4] = local_interped_extended[:, len(LON):len(LON)+4]   
    ##local_interped_masked = np.ma.masked_equal(local_interped, value=0)
    #local_timeseries[step,:,:] = local_interped
    local_timeseries = delayed(get_original_grid_local)(local_interped_extended, LON, step, local_timeseries)
    
    #return local_interped, local, non_local


total = delayed(write_netcdf)( non_local_timeseries, local_timeseries, f, UNIT, LAT_DIFF, NSTEPS, NLAT, NLON)
total.compute()

    #if False:
    
        #fig = plt.figure()  
        #plt.imshow(non_local, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(local, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()

        #fig = plt.figure()  
        #plt.imshow(local_interped, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(local_interped_masked, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(SEA_LAND_MASK, cmap="Wistia",interpolation= 'None')
        ##cb = plt.colorbar(im, orientation='horizontal')
        #plt.colorbar()
        #plt.show()    
        
        #fig = plt.figure()  
        #plt.imshow(LAND_EXTENDED, cmap="Wistia",interpolation= 'None')
        ##cb = plt.colorbar(im, orientation='horizontal')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(CHESSBOARD_EXTENDED, cmap="Wistia",interpolation= 'None')
        ##cb = plt.colorbar(im, orientation='horizontal')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(GLOBAL_CHESSBOARD_EXTENDED, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(GLOBAL_CHESSBOARD_EXTENDED_REVERSE, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(coastboxes, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()
        
        #fig = plt.figure()  
        #plt.imshow(coastboxes_2, cmap="Wistia", interpolation= 'None')
        #plt.colorbar()
        #plt.show()



    #if False:
        #fig = plt.figure()
        #plt.subplot(2,2,1)
        #plt.imshow(local_interped, interpolation= 'None')#,vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.subplot(2,2,2)
        #plt.imshow(local_extended, interpolation= 'None')#,vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.subplot(2,2,3)
        #plt.imshow(local_interped_extended, interpolation= 'None')#,vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.subplot(2,2,4)
        #plt.imshow(CHESSBOARD_EXTENDED, interpolation= 'None')#,vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.show()    

    #if False:
        #fig = plt.figure()
        #plt.subplot(2,2,1)
        #plt.imshow(local_interped, interpolation= 'None',vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.subplot(2,2,2)
        #plt.imshow(local, interpolation= 'None',vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.subplot(2,2,3)
        #plt.imshow(non_local, interpolation= 'None',vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.subplot(2,2,4)
        #plt.imshow(interped_interp1, interpolation= 'None',vmin=-2,vmax=2)
        #plt.colorbar()
        #plt.show()
        
    
    #m = Basemap(projection='cyl', lon_0 = 0, resolution='c')
    ##meshlon, meshlat = m(LON, LAT)    
    
    #fig = plt.figure()  
    ##m.drawcoastlines(linewidth=0.5)
    ##m.drawparallels(np.arange(-20,61,20),labels=[0,1,0,0],fontsize=20, color='grey')
    
    ##m.drawparallels(np.arange(-80,81,40))#, labels=[0,1,0,0])
    ##m.drawmeridians(np.arange(-160,160,20), labels=[0,0,0,1])#,fontsize=8,linewidth=0,labelstyle="+/-")
    
    ##m.pcolormesh( LON, LAT, local_interped + non_local, latlon=True, cmap=cmap_bwr_jwi,\
        ##norm=norm,vmin=-plotmax, vmax=plotmax, rasterized=True, zorder=1)    
    
    ##im = m.pcolormesh( LON, LAT, local_interped, latlon=True, cmap=cmap_bwr_jwi,\
        ##norm=norm,vmin=-plotmax, vmax=plotmax, rasterized=True, zorder=1)    
    
    #im = m.pcolormesh(LON, LAT, CHESSBOARD, cmap='Wistia', latlon=True)# rasterized=True, zorder=1
    #im = m.pcolormesh(LON, LAT, differencemap, cmap='Wistia', latlon=True)# rasterized=True, zorder=1
    #im = m.pcolormesh(LON, LAT, non_local, cmap='Wistia', latlon=True)#, zorder=1)#, rasterized=True,
    #im = m.pcolormesh(LON, LAT, local, cmap='Wistia', latlon=True)#, zorder=1)#, rasterized=True, 
    #im = m.pcolormesh(LON, LAT, local_interped, cmap='Wistia', latlon=True)#, rasterized=True, zorder=1)
    
    #plt.show()
    #im = m.pcolormesh(ii - 0.05, jj - 0.05, interpolated_grid, cmap='bwr', vmin= -5, vmax= 5)
    #cb = plt.colorbar(im, orientation='horizontal')
    #plt.show()
