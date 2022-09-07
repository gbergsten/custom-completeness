import numpy as np
import os

# Creating the period-radius grid used for the detection efficiency maps
def get_Period_Radius_grid():
    nP= 21
    nR= 25
    Pgrid = np.geomspace(0.5, 730, nP)
    Rpgrid = np.geomspace(0.3, 20, nR)
    
    return Pgrid, Rpgrid
     
def make_Detection_Efficiency_Map(stars):
    
    Pgrid, Rpgrid = get_Period_Radius_grid()
        
    # Reading in files of detection metrics from KeplerPORTS
    # Split into two files to get around GitHub file size limitations
    part1 = np.load('files/KeplerPorts_DetectionEfficiency_part1.npz')
    part2 = np.load('files/KeplerPorts_DetectionEfficiency_part2.npz')
    
    idlist = np.concatenate([part1['KID'], part2['KID']])
    det3D = np.concatenate([part1['det3D'], part2['det3D']])
    
    _, use_deteff, use_stars = np.intersect1d(idlist, stars['KIC'], return_indices=True)
        
    ### Calculate the detection efficiency
    fsnr= np.nanmean(det3D[use_deteff,:,:],axis=0)
    
    ### Setting up file for EPOS
    mean_Mstar = np.mean(stars[use_stars]['Mass'])
    mean_Rstar = np.mean(stars[use_stars]['Rad'])
    
    return Pgrid, Rpgrid, fsnr.T, mean_Mstar, mean_Rstar, len(stars)