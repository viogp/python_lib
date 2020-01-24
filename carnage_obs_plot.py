#! /usr/bin/env python

import numpy as np

def plot_smf(ax,obsdir,zz,h0):
    file = obsdir+'smf_z'+zz+'.dat'
    ml,mh,bp3,bdp3 = np.loadtxt(file,unpack=True)
    lm = 0.14 + (ml+mh)/2. - 2.*np.log10(h0)
    mmin = 9.31 + 0.14 - 2.*np.log10(h0)
    mmax = 10.8152 + 0.14 - 2.*np.log10(h0)

    mass = 10.**lm
    # Check again the following lines!!!!!!!!!
    p3 = bp3*h0**3.
    dp3 = bdp3*h0**3.
    err = dp3/(p3*np.log(p3))

    lobs = "Obs, z="+zz
    # Plot all the points
    xobs = lm

    yobs = xobs*0. - 999.  ; indx = np.where( p3 > 0)
    yobs[indx] = np.log10(p3[indx]) 

    #herr = yobs*0. + 999.  ; indx = np.where( (p3+dp3) > 0)
    #herr[indx]  = np.log10(p3[indx] + dp3[indx])

    ax.errorbar(xobs, yobs, yerr=err, ls='None', mfc='None', ecolor = 'LightGrey', mec='LightGrey',marker='o')

# Plot only the points to be used for calibration clearly
    ind = np.where((lm>=mmin)&(lm<=mmax))
    xa = lm[ind]
    ya = yobs[ind] ; erra = err[ind]

    ax.errorbar(xa, ya, yerr=erra, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label=lobs)


def plot_coldgas(ax,obsdir,zz):
    file = obsdir+'cgmf_z'+zz+'.dat'
    xobs,y1,err1,y2,err2 = np.loadtxt(file,unpack=True)
    lobs = "Bosselli+14, z="+zz

    ax.errorbar(xobs, y1, yerr=err1, ls='None', mfc='None', ecolor = 'grey', mec='LightGrey',marker='o',label=lobs)
    ax.errorbar(xobs, y2, yerr=err2, ls='None', mfc='None', ecolor = 'grey', mec='Grey',marker='o')


def plot_sfrf(ax,obsdir,zz):
    file = obsdir+'sfrf_z'+zz+'.dat'
    ml,mh,p3,dp3 = np.loadtxt(file,unpack=True)
    lm = (ml+mh)/2.
    lobs = "Gruppionni+15, z=[0.0,0.3]"
    # Plot all the points
    xobs = lm
    yobs = p3
    err = dp3

    ax.errorbar(xobs, yobs, yerr=err, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label=lobs)


def plot_bhbulge(ax,obsdir,zz):
    file = obsdir+'bh_bulge_av_z'+zz+'.dat'
    ml,mh,p3,dp3 = np.loadtxt(file,unpack=True)
    lm = (ml+mh)/2.
    lobs = "McConnell+13, Kormendy+13"
    # Plot all the points
    xobs = lm
    yobs = p3
    err = dp3

    ax.errorbar(xobs, yobs, yerr=err, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label=lobs)


