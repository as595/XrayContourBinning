#
#  contour_binning.py
#
#
#  Created by Anna Scaife on 02/06/2021.
# -------------------------------------------------------------------------------------------------------


import numpy as np
import pylab as pl
import warnings
import copy
from tqdm import tqdm
import scipy.ndimage as ndimage

import argparse
import configparser as ConfigParser
import ast

#verbose=True

# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------

def get_centroid(data, mask, region):
    
    tmp = np.zeros(mask.shape)
    tmp[np.where(mask!=region)] = 0
    tmp[np.where(mask==region)] = 1
    
    pix = np.argwhere(tmp>0)
    flux = data[np.where(tmp>0)]
    
    pix_x = pix[:,0]
    pix_y = pix[:,1]
    
    xcen = np.sum(pix_x*flux)/np.sum(flux)
    ycen = np.sum(pix_y*flux)/np.sum(flux)
    
    return xcen,ycen

# -------------------------------------------------------------------------------------------------------

def contour_binning(data, maxbin=50, threshold=70, constraint=2, maxrad=None, xraycen=None, verbose=False):

    """
    Parameters:
    
    data       : image data array in units of counts/s/deg2
    maxbin     : maximum number of countour bins
    threshold  : minimum counts per contour bin
    constraint : geometry constraint d/R < C
    maxrad     : maximum radius within which to contour bin [pixels]
    xraycen    : xray centroid for maxrad constraint [pixels]
    verbose    : flag on print output
    """

    mask = np.zeros(data.shape)
    bins=[]
    
    if maxrad!=None and maxrad>0.:
        if xraycen==None:
            print("Maximum radius for binning requires that the Xray centre is specified. Please specify this in the config file or set 'maxrad: None'.")
            exit()
    
    for i in tqdm(range(maxbin)):
        
        tmp = np.zeros(data.shape)
        tmp[:,:] = data[:,:]
        tmp[np.where(mask>0)] = 0
        
        i0 = np.unravel_index(np.argmax(tmp, axis=None), data.shape)
        
        # find first pixel in bin:
        pix_0 = data[i0]
        
        bin_val = pix_0
        bin_g   = (1 + np.sqrt(pix_0+0.75))**2
        
        mask[i0] = i+1
        
        while True:
            
            nomorepix=False
            
            # find offset pixels:
            ofset = np.array([[1,0],[0,1]])
            idx = np.argwhere(mask==i+1)
            pls = np.array([ix+ofset for ix in idx[:]]).reshape(2*idx.shape[0],2)
            mns = np.array([ix-ofset for ix in idx]).reshape(2*idx.shape[0],2)
            sides = np.squeeze(np.vstack((pls,mns)))
            
            # check if outside fov:
            rm=[]
            for j in range(sides.shape[0]):
                if sides[j,0]>=data.shape[0] or sides[j,1]>=data.shape[1] or sides[j,0]<0 or sides[j,1]<0:
                    rm.append(j)
            
            sides = np.delete(sides,(rm), axis=0)
            
            # check if already masked:
            rm=[]
            for j in range(sides.shape[0]):
                if mask[sides[j,0],sides[j,1]]>0:
                    rm.append(j)
            
            sides = np.delete(sides,(rm), axis=0)
            
            if sides.shape[0]>0:
            
                # get values of offsets:
                tmp = data[sides[:,0],sides[:,1]]
                
                pix_i = sides[np.argsort(np.abs(tmp-pix_0))]
                
                j=0
                while True:
                    
                    if j>=len(pix_i):
                        nomorepix = True
                        break
                    else:
                        tmp_p = pix_i[j]
                
                    # check length to width ratio:
                    r_parm = np.sqrt(np.count_nonzero(mask==i+1)/np.pi)
                    xcen, ycen = get_centroid(data,mask,region=i+1)
                    d_x = tmp_p[0] - xcen
                    d_y = tmp_p[1] - ycen
                    d1 = np.sqrt(d_x**2 + d_y**2)
                    
                    # check maximum radius constraint:
                    if maxrad!=None and maxrad>0.:
                        d_x = tmp_p[0] - xraycen[0]
                        d_y = tmp_p[1] - xraycen[1]
                        d2 = np.sqrt(d_x**2 + d_y**2)
                        
                    if d1/r_parm < constraint and d2 < maxrad:
                        new_pix = pix_i[j]
                        break
                    else:
                        j+=1
                        
                
                # update mask:
                mask[new_pix[0],new_pix[1]] = i+1

                # update bin value:
                bin_val += data[new_pix[0],new_pix[1]]
                bin_g += (1 + np.sqrt(data[new_pix[0],new_pix[1]]+0.75))**2

                # check threshold:
                #if bin_val/np.sqrt(bin_g) >= threshold:
                if bin_val >= threshold:
                    if verbose:
                        tqdm.write("Bin {} sum: {}".format(i+1,bin_val))
                    bins.append(bin_val)
                    break

            else:
                if verbose:
                    tqdm.write("Bin {} sum: {}".format(i+1,bin_val))
                bins.append(bin_val)
                break
                
            if nomorepix:
                if verbose:
                    tqdm.write("Bin {} sum: {}".format(i+1,bin_val))
                break
            
                
    return mask, bins

# -------------------------------------------------------------------------------------------------------

def sum_func(values):
    return values.sum()
    
    
def find_edge(mask,region):

    tmp = np.zeros(mask.shape)
    tmp[np.where(mask!=region)] = 0
    tmp[np.where(mask==region)] = 1
    
    footprint = np.array([[0,1,0],[1,1,1],[0,1,0]])
    edge = ndimage.generic_filter(tmp, sum_func, footprint=footprint, mode='constant', cval=0.0)
    
    edge[np.where(tmp==0)]  = 0
    edge[np.where(edge==5)] = 0
    edge[np.where(edge>0)]  = 1
    
    return edge

# -------------------------------------------------------------------------------------------------------

def get_diflist(data, mask, edge, region):

    # find offset pixels:
    ofset = np.array([[1,0],[0,1]])
    idx = np.argwhere(edge==1)
    
    dlist=[]
    for ix in idx:
    
        flux = data[ix[0],ix[1]]

        pls = (ix+ofset)
        mns = (ix-ofset)
        sides = np.squeeze(np.vstack((pls,mns)))
           
        # check if outside fov:
        rm=[]
        for j in range(sides.shape[0]):
            if sides[j,0]>=data.shape[0] or sides[j,1]>=data.shape[1] or sides[j,0]<0 or sides[j,1]<0:
                rm.append(j)
                
        sides = np.delete(sides,(rm), axis=0)
                
        # check if inside region
        rm=[]
        for j in range(sides.shape[0]):
            if mask[sides[j,0],sides[j,1]]==region:
                rm.append(j)
                
        sides = np.delete(sides,(rm), axis=0)
                
        for side in sides:
            diff = np.abs(data[side[0],side[1]] - flux)
            dlist.append([ix,side,diff])

    return np.array(dlist, dtype=object)
    
# -------------------------------------------------------------------------------------------------------

def remake_bins(data,mask):

    nbin = int(np.max(mask))
    bins = np.zeros(nbin)
        
    for i in range(0,nbin):
    
        # update bin value:
        bin_val = np.sum(data[np.where(mask==i+1)])
        bin_g = np.sqrt(np.sum((1 + np.sqrt(data[np.where(mask==i+1)]+0.75))**2))
        if bin_val>0:
            bins[i] = bin_val/bin_g
        else:
            bins[i] = 0
            
    return bins
    
# -------------------------------------------------------------------------------------------------------

def clean_bins(data, mask, bins, threshold=70.):

    idx = np.argsort(bins)
    
    for ix in idx:
    
        if bins[ix]> 0 and bins[ix]<threshold:
        
            npix = np.count_nonzero(mask==ix+1)
        
            while npix>0:
                # get bin edge
                edge = find_edge(mask,region=ix+1)
                # get differences by edge:
                dlist= get_diflist(data, mask, edge, region=ix+1)
                imin = np.argmin(dlist[:,2])
                i0 = dlist[imin,0]
                i1 = dlist[imin,1]
                # update mask
                mask[i0[0],i0[1]] = mask[i1[0],i1[1]]
                npix-=1
                
    bins = remake_bins(data, mask)

    return mask, bins

# -------------------------------------------------------------------------------------------------------

def reorder_bins(mask, bins):
    
    nbin = len(bins)
    for i in range(nbin):
        ib = nbin-i
        ix = nbin-i-1
        if (ix-1)>=0 and bins[ix-1]==0.:
            mask[np.where(mask==ib)] -= 1
            bins[ix-1] = bins[ix]
            bins[ix] = 0.
        
    return mask, bins

# -------------------------------------------------------------------------------------------------------

def parse_args():
    """
        Parse the command line arguments
        """
    parser = argparse.ArgumentParser()
    parser.add_argument('-C','--config', default="myconfig.txt", required=True, help='Name of the input config file')
    
    args, __ = parser.parse_known_args()
    
    return vars(args)

# -------------------------------------------------------------------------------------------------------

def parse_config(filename):
    
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.read(filename)
    
    # Build a nested dictionary with tasknames at the top level
    # and parameter values one level down.
    taskvals = dict()
    for section in config.sections():
        
        if section not in taskvals:
            taskvals[section] = dict()
        
        for option in config.options(section):
            # Evaluate to the right type()
            try:
                taskvals[section][option] = ast.literal_eval(config.get(section, option))
            except (ValueError,SyntaxError):
                err = "Cannot format field '{0}' in config file '{1}'".format(option,filename)
                err += ", which is currently set to {0}. Ensure strings are in 'quotes'.".format(config.get(section, option))
                raise ValueError(err)

    return taskvals, config

# -------------------------------------------------------------------------------------------------------

def order_edge(edge):
    
    poly_pix = np.argwhere(edge==1)
    nvert = poly_pix.shape[0]
    mask_pix = np.ones(nvert,dtype=bool)
    idx=[]
    i=0; ix=0
    while True:
        pix = poly_pix[ix,:]
        seps = np.sum((poly_pix - pix)**2, axis=1)
        mask_pix[ix] = False
        seps[~mask_pix] = 1e10
        ix = np.argmin(seps)
        if seps[ix]>10:
            #print("Finished polygon")
            break
        idx.append(ix)
        i+= 1
        if (i>=nvert): break

    pts = poly_pix[idx,:]
    
    return pts

# -------------------------------------------------------------------------------------------------------

def ds9_region(filename, poly_sky):

    polygon = 'polygon('+','.join(' '.join(poly_sky.to_string()).split())+')'
    with open(filename, 'w') as f:
        f.write('# Region file format: DS9 version 4.1 \n')
        f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        f.write('fk5 \n')
        f.write(polygon)

    return
    
# -------------------------------------------------------------------------------------------------------

