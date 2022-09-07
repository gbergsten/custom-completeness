import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy import constants as c

from detection import *
from vetting import *

def get_completeness(ranges={}, score_cut=0.0, dwarfcut=True, Verbose=False, return_components=False):
    
    # Load up data files
    # For stellar properties, use Berger et al. 2020a Gaia-Kepler Steller Catalog
    stars = Table.read('files/Berger2020a_GaiaKeplerCatalog_Table2.fits');
    # Kepler detection metrics
    DetectionMetrics = np.load('files/metricslist.npz')
    
    # Get the subset of stars with both detection metrics and Gaia-updated properties
    _, _, use = np.intersect1d(DetectionMetrics['KID'], stars['KIC'], return_indices=True)
    stars = stars[use]
    
    # Apply a cut to select only dwarf stars
    if dwarfcut:
        isdwarf = stars['logg'] > 1./4.671 * np.arctan((stars['Teff']-6300.)/-67.172)+ 3.876
        stars = stars[isdwarf]
    
    # Filter the stellar catalog by any parameter ranges specified
    # Format for 'ranges' entries must use same name as Berger catalog
    use = np.full(len(stars),True)
    for key in ranges.keys(): # for each parameter specified,
        # Find entries that fall within the given range
        inRange = (stars[key] >= ranges[key][0]) & (stars[key] < ranges[key][1])
        # Keep a running list of which stars have met all properties
        use = np.logical_and(use, inRange)
    # Update table to only have desired stars.
    stars = stars[use] 
    
    # Calculate the various components of completeness (following the methods of epos - Mulders et al. 2018).
    
    # Create detection efficiency map using this sample
    Pgrid, Rpgrid, det_eff, mean_Mstar, mean_Rstar, n_stars = make_Detection_Efficiency_Map(stars)
    
    # Fit the vetting efficiency parameters for the given score cut
    vetpars = make_Vetting_Efficiency_Parameters(stars, score_cut, Verbose)
    # Make a grid using those parameters
    X,Y = np.meshgrid(Pgrid, Rpgrid, indexing='ij')
    vet_eff = fbpl2d( (X,Y), *vetpars)
    
    # Compute the geometric transit probability map
    fourpi2_GM = 4. * np.pi**2. / (c.G * mean_Mstar * u.Msun)
    fgeo_prefac = mean_Rstar * u.Rsun * fourpi2_GM**(1./3.)
    fgeo = (fgeo_prefac * (X * u.day)**(-2/3)).si
    
    # Calculate the total completeness
    completeness = det_eff * vet_eff * fgeo
    
    # Can return the three components of total completeness if specified.
    if return_components:
        return Pgrid, Rpgrid, completeness, n_stars, [det_eff, vet_eff, fgeo]
    
    return Pgrid, Rpgrid, completeness, n_stars
    
    
    
    
