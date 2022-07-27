"""
Astronomical functions
"""

import numpy as np
from astropy.constants import c,M_sun,L_sun

c_ms = c.value
Msun_kg = M_sun.value

J2erg = 10**7
s_in_year = 365.*24.*3600.

er_mcc16 = 0.1
ef_mcc16 = 0.15

def lbol_coeff(x):
    match x:
        case 1:
            return er_mcc16*(1-ef_mcc16)
        case _:
            return 0   # 0 is the default case if x is not found

def get_lbol(lmdotbh,mdot2SI=True,eq=1,units='SI'):
    """
    Obtain the bolometric luminosity given the BH accretion rate

    Parameters
    ----------
    lmdotbh : numpy array of floats
      log10(Mdot_BH/Msun/yr) or log10(Mdot_BH/kg/s)
    mdot2SI : boolean
      True if the input given is log10(Mdot_BH/Msun/yr)
    eq : integer
      Specifies the equation to be used
      1 corresponds to that from McCarthy+2016
    unitsOUT : string
      SI for W; cgs for erg/s, sun for L_sun

    Returns
    -------
    lbol : numpy array of floats
        log10(Lbol) with units depending on flags
    """

    if (mdot2SI):
        lmdotbh = lmdotbh + np.log10(Msun_kg) - np.log10(s_in_year)

    case of     
    lLEdd = np.log10(1.26) + 46. + lmbh - 8.
    return lLEdd
    

def get_lLEdd(lmbh):
    """
    Obtain the Eddington luminosity limit, following Eq. 3 in Griffin+2020

    Parameters
    ----------
    lmbh : numpy array of floats
        log10(M_BH/Msun)

    Returns
    -------
    lLedd : numpy array of floats
        log10(LEdd/erg/s)
    """

    lLEdd = np.log10(1.26) + 46. + lmbh - 8.
    return lLEdd


def get_lmdotEdd(lmbh,SI=True,h0=None):
    """
    Calculate the Eddington BH accretion limit, following Eq. 4 in Griffin+2020.

    Parameters
    ----------
    lmbh : numpy array of floats
       BH mass as log10(M_BH).
    SI : boolean
       If true, the output log10(Mdot_Edd/kg/s), 
       otherwise log10(Mdot_Edd/Msun/year).        
    h0 : float
       Hubble constant unit, H=100h. 
       If None, the input is assumed log10(M_BH/Msun),
       otherwise log10(M_BH/Msun/h0).

    Returns
    -------
    lmdotedd : numpy array of floats
       Eddington limit for the mass accretion, 
       log10(Mdot_Edd/kg/s) if SI True
       log10(Mdot_Edd/Msun/yr) if SI False
       log10(Mdot_Edd/Msun/h0/yr) if SI False and given h0
    """

    if (h0 is not None):
        lmbh = lmbh - np.log10(h0)

    ledd = get_lLEdd(lmbh) - 7. #W
    lmdotedd = ledd + 1 - 2.*np.log10(c_ms)

    if (not SI):
        lmdotedd = lmdotedd - np.log10(Msun_kg) + np.log10(s_in_year)
        if (h0 is not None):
            lmdotedd = lmdotedd + np.log10(h0)
        
    return lmdotedd


if __name__ == "__main__":
    lmbh = np.array([7.,8.])
    print(get_lLEdd(lmbh))
    print(get_lmdotEdd(lmbh,SI=True,h0=0.73))
