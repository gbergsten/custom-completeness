import numpy as np
import os
from astropy.io import ascii
from scipy.optimize import curve_fit

from detection import *

# Functional form to fit the vetting efficiency grid
def fbpl2d(xy, a, b, c, d, e, f, g):
    x,y= xy
    bpl= a* (x/b)**np.where(x<b, c, d) * (y/e)**np.where(y<e, f,g)
    return np.maximum(0.2, np.minimum(bpl, 1.))

def make_Vetting_Efficiency_Parameters(stars, score_cut, Verbose):
    
    # Load the Injected TCEs file
    injTces = ascii.read("files/kplr_dr25_lc/kplr_dr25_inj1_tces.txt")
    if Verbose:
        print("num injected TCEs: " + str(np.size(injTces)))
    
    # TCE Drop List
    tcelist = "files/kplr_dr25_lc/DR25-Injected-Recovered-OnTarget-Planet-TCEs-1-1-Prat.txt"
    tces = np.genfromtxt(tcelist, dtype='str')

    tceKepids = np.zeros(len(tces));
    for i in range(len(tces)):
        s = tces[i].split('-');
        tceKepids[i] = int(s[0]);
        
    if Verbose:
        print("num injected/recovered TCEs: " + str(np.size(tceKepids)))
    
    injTces = injTces[np.in1d(injTces['TCE_ID'],tces)]
    if Verbose:
        print("num injected TCEs after trimming to injected/recovered: " + str(np.size(injTces)))

    # Select only those TCEs that are in this stellar population
    injTces = injTces[np.in1d(injTces['KIC'],stars['KIC'])]
    if Verbose:
        print("after: " + str(np.size(injTces)))

    # Do some basic stats
    if Verbose:
        print("# of injected TCEs: " + str(len(injTces)))
        print("# of injected PCs: " + str(len(injTces[injTces['Disp']=='PC'])))
        print("# of injected FPs: " + str(len(injTces[injTces['Disp']=='FP'])))
        print(' ')
    
    # Finding all the TCE KICs with Kepler-Gaia info
    kepIDs, tce_idx, star_idx = np.intersect1d(injTces['KIC'], stars['KIC'], return_indices=True)

    # Making a new pair of arrays that have the same KIC-index order
    stars = stars[star_idx]
    injTces = injTces[tce_idx]
    
    isPC = injTces['Disp']=='PC'
    isFP = injTces['Disp']=='FP'

    # Calculate Gaia-updated TCE radius based on updated stellar radii
    Rearth = 6371 #km
    Rsun = 696340 #km
    new_planet_radii = (injTces['Rp/Rs'] * stars['Rad'] * Rsun) / Rearth
    
    # Making a referential array to store key values in
    vet_star = {}
    vet_star['ID'] = np.asarray(injTces['KIC'])
    vet_star['isPC'] = np.asarray(isPC)
    vet_star['score']= np.asarray(injTces['Score'])
    vet_star['Teff']= np.asarray(stars['Teff'])
    vet_star['logg']= np.asarray(stars['logg'])
    vet_star['Rs']= np.asarray(stars['Rad'])
    vet_star['Mstar']= np.asarray(stars['Mass'])
    vet_star['Rp']= np.asarray(new_planet_radii)
    vet_star['P']= np.asarray(injTces['period'])
    
    # Remove injections outside grid boundaries 
    Pgrid, Rpgrid = get_Period_Radius_grid()
    ingrid = (Rpgrid[0] <= vet_star['Rp']) & (vet_star['Rp'] < Rpgrid[-1]) \
                & (Pgrid[0] <= vet_star['P']) & (vet_star['P'] < Pgrid[-1])
    for key in vet_star:
        vet_star[key]= vet_star[key][ingrid]
    grid2d = (Pgrid, Rpgrid)
        
    # Fitting the Vetting Params
    Rst = np.average(vet_star['Rs'])
    Rst_median = np.median(vet_star['Rs'])
    nst = len(stars)
    if Verbose:
        print('{}, Rst = {:.2f}, median = {:.2f}'.format(nst, Rst, Rst_median))
    
    isPC = vet_star['isPC']
    isrel = vet_star['isPC'] & (vet_star['score'] > score_cut)
    if Verbose:
        print('{} Injections'.format(isPC.size))
        print('{} Candidates, {:.0%}'.format(isPC.sum(), 1.*isPC.sum()/isPC.size))
        print('{} w/ score > {:.2f}, {:.0%}'.format(isrel.sum(), score_cut, 1.*isrel.sum()/isPC.size))
    
    # 2D Histogram
    all_2d, _, _ = np.histogram2d(vet_star['P'],vet_star['Rp'], bins=grid2d)
    cand_2d, _, _ = np.histogram2d(vet_star['P'][isPC],vet_star['Rp'][isPC], bins=grid2d)
    score_2d, _, _ = np.histogram2d(vet_star['P'][isrel],vet_star['Rp'][isrel], bins=grid2d)
    
    # 1D Histogram
    score_P, _ = np.histogram(vet_star['P'][isrel], bins=Pgrid)
    score_R, _ = np.histogram(vet_star['Rp'][isrel], bins=Rpgrid)
    all_P, _ = np.histogram(vet_star['P'], bins=Pgrid)
    all_R, _ = np.histogram(vet_star['Rp'], bins=Rpgrid)
    
    if Verbose:
        print('{} x {} = {} (apparently)'.format(Pgrid.size, Rpgrid.size, all_2d.shape))
    Pc= np.sqrt(Pgrid[:-1]*Pgrid[1:])
    Rc= np.sqrt(Rpgrid[:-1]*Rpgrid[1:])
    
    with np.errstate(invalid='ignore'): 
        fvet= cand_2d / all_2d
        fscore=  score_2d / all_2d
    with np.errstate(invalid='ignore', divide='ignore'): 
        sqrtn= np.sqrt(score_2d)
        ferr= np.where(sqrtn<1,1./all_2d,fscore/sqrtn)

    with np.errstate(divide='ignore'): 
        fscore_R= 1.* score_R / all_R
        fscore_P= 1.* score_P / all_P
        
    Xc,Yc = np.meshgrid(Pc, Rc, indexing='ij')
    
    score= np.ravel(fscore)
    ncell= np.ravel(all_2d)
    Pravel= np.ravel(Xc)
    Rravel= np.ravel(Yc)
            
    nonzero= np.isfinite(score)
    xx_obs=(Pravel[nonzero], Rravel[nonzero])
    yy_obs= score[nonzero]
    n_obs= ncell[nonzero]	
                    
    # fit, all bins equal weight (??)
    fitargs= {}
    fitargs['sigma']= np.ravel(ferr)[nonzero]
    
    # Some initial guesses to help with fitting
    if score_cut==0.9:
        p0 = [0.9, 50, -0.07, -0.4, 5.7, 0.1, -2.7]
    else:
        p0 = [0.9, 100, 0, -0.2, 5.7, 0.1, -2.7]
    
    popt, pcov = curve_fit(fbpl2d, xx_obs, yy_obs, p0= p0, **fitargs)

    perr = np.sqrt(np.diag(pcov))
    if Verbose:
        print('\nFitted parameters:')
    
    parstring = 'vetpars= '
    this_fit = []
    
    for i, (_popt, _perr) in enumerate(zip(popt, perr)):
        if i==6:
            parstring = parstring + '{:.2g}, '.format(popt[5])
            this_fit.append(popt[5])
        else:
            parstring = parstring + '{:.2g}, '.format(_popt)
            this_fit.append(_popt)
    
    if Verbose:
        print(parstring)
    VettingParams = this_fit

    if Verbose:
        print('Vetting Pars: {}'.format(VettingParams))
        print('---')
    
    return VettingParams