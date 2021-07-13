import astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyWarning

import os,sys
import numpy as np
import pylab as pl
import json

import warnings
warnings.simplefilter('ignore', category=AstropyWarning)

from PyContourBin.contour_binning import *

pl.rcParams['figure.figsize'] = [16,16]
pl.rcParams['figure.dpi'] = 150

# -----------------------------------------------------------
# -----------------------------------------------------------

def main(config_dict):

    # read information from config file:
    
    obs_image      = config_dict['data']['obs_image']
    exptime        = config_dict['data']['exptime']
    t              = config_dict['binning']['threshold']
    maxbin         = config_dict['binning']['maxbin']
    constraint     = config_dict['binning']['constraint']
    xmmsas         = config_dict['sas']['issas']
    plots          = config_dict['outputs']['plots']
    plotdir        = config_dict['outputs']['plotdir']
    regdir         = config_dict['outputs']['regdir']
    
    if plotdir==None: plotdir = './'
    if regdir==None: regdir = './'
    
    if not os.path.exists(plotdir): os.mkdir(plotdir)
    if not os.path.exists(regdir): os.mkdir(regdir)
    
    if config_dict['binning']['maxrad']!=None:
        maxrad = config_dict['binning']['maxrad']
        try:
            xraycen = config_dict['binning']['centre']
        except:
            print("Maximum radius for binning requires that the Xray centre is specified. Please specify this in the config file or set 'maxrad: None'.")
            exit()
            
    # read adaptively smoothed image:
    # should be in units of counts/s/deg2

    obs_hdu = fits.open(obs_image)[0]
    obs_wcs = WCS(obs_hdu.header)

    data = obs_hdu.data
    
    pixsize = obs_hdu.header['CDELT2'] # degrees
    pixarea = pixsize**2

    censky = SkyCoord(xraycen.split()[0], xraycen.split()[1], frame='icrs')
    cenpix = obs_wcs.world_to_pixel(censky)
        
    data = obs_hdu.data
    if maxrad>0.:
        maxradpix = maxrad/(pixsize*60.)
        h, w = data.shape
        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - cenpix[0])**2 + (Y-cenpix[1])**2)
        data[np.where(dist_from_center>maxradpix)] = 0

    # plot the input data:
    
    if plots:
        pl.subplot(111, projection=obs_wcs)
        pl.imshow(data, origin='lower', vmin=0, vmax=30.0, cmap='cubehelix')
        pl.grid(color='white', ls='solid')
        pl.xlabel('Right Ascension [J2000]')
        pl.ylabel('Declination [J2000]')
        pl.savefig(plotdir+'input_masked.png')
     
    # run the contour binning initial pass:

    threshold = t/(exptime*pixarea)
    mask, bins = contour_binning(data, maxbin=maxbin, threshold=threshold, constraint=constraint, maxrad=maxradpix, xraycen=cenpix, verbose=False)

    # plot the initial bins:
    
    if plots:
        pl.subplot(111, projection=obs_wcs)
        pl.imshow(mask, origin='lower', cmap='tab20c')
        pl.grid(color='white', ls='solid')
        pl.xlabel('Right Ascension [J2000]')
        pl.ylabel('Declination [J2000]')
        pl.savefig(plotdir+'initial_bins.png')

    # run the cleaning pass:

    initial_mask = mask.copy()
    initial_bins = bins.copy()
    mask, bins = clean_bins(data, mask, bins, threshold=threshold)

    # plot the cleaned bins:
    
    if plots:
        pl.subplot(111, projection=obs_wcs)
        pl.imshow(mask, origin='lower', cmap='tab20c')
        pl.grid(color='white', ls='solid')
        pl.xlabel('Right Ascension [J2000]')
        pl.ylabel('Declination [J2000]')
        pl.savefig(plotdir+'cleaned_bins.png')

    # reorder the bins:

    cleaned_mask = mask.copy()
    cleaned_bins = bins.copy()
    mask, bins = reorder_bins(mask, bins)

    # create the regions in sky coordinates:

    for i in tqdm(range(len(bins))):
        edge = find_edge(mask,region=i+1)
        if np.count_nonzero(edge)>0:
            poly_pix = np.argwhere(edge==1)
            poly_sky = obs_wcs.pixel_to_world(poly_pix[:,1], poly_pix[:,0])
            
            # JSON output:
            polygon = {'edge_'+str(i+1) : poly_sky.to_string()}
            with open(regdir+'edge_'+str(i+1)+'.json', 'w') as f:
                json.dump(polygon, f)
                
            # DS9 output:
            filename = regdir+'edge_'+str(i+1)+'.reg'
            ds9_region(filename, poly_sky)
            
            # SAS output:
            #if sas:
            #    make_sas_region(polygon)
                

    return
    
    
# -----------------------------------------------------------
# -----------------------------------------------------------

if __name__ == "__main__":

    vars = parse_args()
    config_dict, config = parse_config(vars['config'])

    main(config_dict)
