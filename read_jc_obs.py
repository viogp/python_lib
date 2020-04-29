#! /usr/bin/env python

import numpy as np
import sys,glob

def set_h0(h0=None):
    global obsh0
    if(h0 is None):
        obsh0 = 0.677
    else:
        obsh0 = h0
    return  

def set_dz(dz=None):
    global obsdz
    if(dz is None):
        obsdz = 0.2
    else:
        obsdz = dz
    return


def read_jc_lf(obs_dir,zz,dz=None,h0=None,infile=None,survey=None,limit=None,inlog=True):
    set_h0(h0) ; set_dz(dz) ; firstpass=True

    files = [ ]
    if(infile is None):
        files = obs_dir
    else:
        files.append(obs_dir+infile)
   
    for i in range(len(files)):
        fil = files[i]
        oz,ozlow,ozhigh,oL,oPhi,oPlow,oPhigh = np.loadtxt(fil,usecols=(0,1,2,3,6,7,8),comments='#',unpack=True)

        ind = np.where((abs(oz-zz)<obsdz) & (ozlow<zz) & (ozhigh>zz) & (oPhi>0.))
    
        # Johan's units: L(erg/s), Phi(Mpc^-3/dLogL)
        # Output: L(h-2 erg/s), Phi(Mpc^-3h^3/dLogL)
        if (np.shape(ind)[1]>0):
            if inlog:
                if firstpass:
                    firstpass = False
                    L = np.log10(oL[ind]) + 2*np.log10(obsh0) 
                    Phi = np.log10(oPhi[ind]) - 3*np.log10(obsh0)
                    el = np.log10(oPhi[ind])  - np.log10(oPlow[ind])
                    eh = np.log10(oPhigh[ind])- np.log10(oPhi[ind])
                else:
                    L = np.append(L,np.log10(oL[ind]) + 2*np.log10(obsh0))
                    Phi = np.append(Phi,\
                              np.log10(oPhi[ind])- 3*np.log10(obsh0))
                    el = np.append(el,np.log10(oPhi[ind])\
                              - np.log10(oPlow[ind]))
                    eh = np.append(eh,np.log10(oPhigh[ind])\
                              - np.log10(oPhi[ind]))
            else:
                if firstpass:
                    firstpass = False
                    L = oL[ind]*obsh0**2.
                    Phi = oPhi[ind]/(obsh0**3)
                    el = (oPhi[ind]-oPlow[ind])/(obsh0**3)
                    eh = (oPhigh[ind]-oPhi[ind])/(obsh0**3)
                else:
                    L = np.append(L,oL[ind]*obsh0**2.)
                    Phi = np.append(Phi,oPhi[ind]/(obsh0**3))
                    el = np.append(el,(oPhi[ind]-oPlow[ind])/(obsh0**3))
                    eh = np.append(eh,(oPhigh[ind]-oPhi[ind])/(obsh0**3))
        else:
            L, Phi, el, eh = [ -999. for i in range(4)]
            print('No adecuate data for this redshift')

    return L, Phi, el, eh


def read_jc_indlf(obs_dir,zz,dz=None,h0=None,survey=None,line=None,band=None):
    set_h0(h0) ; set_dz(dz) ; firstpass=True
    
    if(survey is None or band is None):
        L, Phi, el, eh = [ -999. for i in range(4)]
        print('No band or observational survey given')
    else:
        if (line is None):
            line = 'O2_3728'

        files = glob.glob(obs_dir+line+'-'+survey+band+"*.txt")
        lf = len(files) ; zf = np.zeros(shape=(lf))
        for i,ifile in enumerate(files):
            zf[i] = float(ifile.split('z')[1].split('.txt')[0])
        if (np.min(abs(zf-zz)) < 0.2):
            i = np.argmin(abs(zf-zz))
            fil = obs_dir+line+'-'+survey+band+'-z'+str(zf[i])+'.txt'
            print('Obs. {}'.format(fil))

            oz,ozlow,ozhigh,oL,oPhi,oPlow,oPhigh = np.loadtxt(fil,usecols=(0,1,2,3,6,7,8),comments='#',unpack=True)

            ind = np.where(oPhi>0.)

            # Johan's units: L(erg/s), Phi(Mpc^-3/dLogL)
            # Output: L(h-2 erg/s), Phi(Mpc^-3h^3/dLogL)
            if (np.shape(ind)[1]>0):
                if firstpass:
                    firstpass = False
                    L = np.log10(oL[ind]) + 2*np.log10(obsh0) 
                    Phi = np.log10(oPhi[ind]) - 3*np.log10(obsh0)
                    el = np.log10(oPhi[ind])  - np.log10(oPlow[ind])
                    eh = np.log10(oPhigh[ind])- np.log10(oPhi[ind])
                else:
                    L = np.append(L,np.log10(oL[ind]) + 2*np.log10(obsh0))
                    Phi = np.append(Phi,np.log10(oPhi[ind])- 3*np.log10(obsh0))
                    el = np.append(el,np.log10(oPhi[ind])\
                                       - np.log10(oPlow[ind]))
                    eh = np.append(eh,np.log10(oPhigh[ind])
                                   - np.log10(oPhi[ind]))
            else:
                L, Phi, el, eh = [ -999. for i in range(4)]
                print('No adecuate data for this redshift')
        else:
            L, Phi, el, eh = [ -999. for i in range(4)]
            print('No adecuate data for this redshift')

    return L, Phi, el, eh

