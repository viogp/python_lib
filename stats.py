#! /usr/bin/env python

"""
Some useful functions
  percentiles(xbins,xarray,yarray,per): obtains percentiles of yarray in xbins

  convert_to_stdev(grid): normalised a grid to cumulative standard deviations.

NOTE: this module requires the numpy and scipy libraries to be
      available for import!

"""
import sys
import numpy as np

def percentiles(val,data,weights=None):
    if (val <0 or val >1):
        sys.exit('STOP percentiles: 0<val<1')

    if (weights is None):
        ws = np.zeros(shape=(len(data))) ; ws.fill(1.)
    else:
        ws = weights

    data = np.array(data) ; ws = np.array(ws)
    ind_sorted = np.argsort(data)  # Median calculation from wquantiles
    sorted_data = data[ind_sorted] ; sorted_weights = ws[ind_sorted]
    
    num = np.cumsum(sorted_weights) - 0.5*sorted_weights 
    den = np.sum(sorted_weights) 
    if (den!=0): 
        pn = num/den   
        percentiles = np.interp(val, pn, sorted_data)  
    else:
        sys.exit('STOP percentiles: problem with weights')
    return percentiles

def perc_2arrays(xbins,xarray,yarray,weights,nmin,val):
    """ Returns percentiles of yarray over xbins"""
    xlen = len(xbins)-1
    perc_2arrays = np.zeros(shape=(xlen)) ; perc_2arrays.fill(-999.)

    if len(xarray) != len(yarray):
        sys.exit('ERROR @ perc_2arrays: The lenght of the input arrays should be equal.')

    for i in range(xlen):
        ind = np.where((xarray >= xbins[i]) & (xarray < xbins[i+1]))
        # We require at least nmin points per bin
        if (np.shape(ind)[1] > nmin): 
            data = yarray[ind] ; ws = weights[ind]
            perc_2arrays[i] = percentiles(val,data,weights=ws)

    return perc_2arrays


def av_2arrays(xbins,xarray,yarray,weights,nmin):
    """ Returns average of yarray over xbins"""
    xlen = len(xbins)-1
    av_2arrays = np.zeros(shape=(xlen)) ; av_2arrays.fill(-999.)

    if len(xarray) != len(yarray):
        sys.exit('ERROR @ perc_2arrays: The lenght of the input arrays should be equal.')

    for i in range(xlen):
        ind = np.where((xarray >= xbins[i]) & (xarray < xbins[i+1]))
        # We require at least nmin points per bin
        num = np.shape(ind)[1]
        if (num > nmin): 
            data = yarray[ind] ; ws = weights[ind]
            dw = ws*data
            av_2arrays[i] = np.sum(dw)/np.sum(ws)

    return av_2arrays

def convert_to_stdev(grid):
    """
    From Coleman Krawczyk
    Based on https://pypi.python.org/simple/astroml/
    Given a grid of values, convert them to cumulative standard deviation.
    This is useful for drawing contours with standard deviations as the levels.
    """
    shape = grid.shape
    grid = grid.ravel()
    # Obtain indices to sort and unsort the flattened array
    i_sort = np.argsort(grid)[::-1]
    i_unsort = np.argsort(i_sort)
    grid_cumsum = grid[i_sort].cumsum()
    grid_cumsum /= grid_cumsum[-1]

    return grid_cumsum[i_unsort].reshape(shape)
