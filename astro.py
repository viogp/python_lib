"""
Astronomical functions
"""

import numpy as np
from scipy.constants import c,pi

J2erg = 10**7

def get_LEdd(mbh):
    """
    Following Eq. 3 in Griffin+2020, obtain the Eddington BH accretion mass 

    Parameters
    ----------
    Ledd : numpy array of floats
    Eddington luminosity ()

    Returns
    -------
    mdotedd : numpy array of floats

    """

    LEdd = mbh
    
    return LEdd

def get_mdotEdd(mbh):
    """
    Following Eq. 4 in Griffin+2020, obtain the Eddington BH accretion mass 

    Parameters
    ----------
    Ledd : numpy array of floats
    Eddington luminosity ()

    Returns
    -------
    mdotedd : numpy array of floats

    """

    #Can I get a conditional output(return mdotedd and ledd if flag true)?
    print(c)
    mdotedd = Ledd/(0.1*c*c)
    return mdotedd


if __name__ == "__main__":
    lmbh = np.array([7])
    print(lmbh,type(lmbh),lmbh[0])
    #print(get_lmdotEdd(lmbh))
