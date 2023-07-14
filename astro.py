"""
Astronomical functions
"""

import numpy as np
from astropy.constants import c,M_sun,L_sun,pc

c_ms    = c.value
Msun_kg = M_sun.value
Lsun_W  = L_sun.value

m_in_pc = pc.value

J2erg = 10**7
s_in_year = 365.*24.*3600.


def get_r(x1,y1,z1,x2,y2,z2):
    """
    Calculate the distance between 2 objects
    
    Parameters
    ----------
    x1,y1,z1 : floats
       Coordinattes of object 1 (subhalo or satellite)
    x2,y2,z2 : floats
       Coordinattes of object 2 (halo or central)

    Returns
    -------
    r : float
       Distnace between the objects
    """

    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    return r


def get_vr(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2):
    """
    Calculate the relative radial velocity of an object
    
    Parameters
    ----------
    x1,y1,z1 : floats
       Coordinattes of object 1 (halo or central)
    x2,y2,z2 : floats
       Coordinattes of object 2 (subhalo or satellite)

    Returns
    -------
    vr : float
       Radial relative velocity
    """

    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    dvx = vx2-vx1
    dvy = vy2-vy1
    dvz = vz2-vz1

    r = get_r(x1,y1,z1,x2,y2,z2)
    
    vr = (dx*dvx + dy*dvy + dz*dvz)/r
    return vr


def get_vt(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2):
    """
    Calculate the relative radial velocity of an object
    
    Parameters
    ----------
    x1,y1,z1 : floats
       Coordinattes of object 1 (halo or central)
    x2,y2,z2 : floats
       Coordinattes of object 2 (subhalo or satellite)

    Returns
    -------
    vt : float
       Tangential relative velocity in the x-y plane (v theta)
    """

    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    dvx = vx2-vx1
    dvy = vy2-vy1
    dvz = vz2-vz1

    vt =  (dx*dvy - dy*dvx)/np.sqrt(dx*dx + dy*dy) 
    return vt



def get_vg(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2):
    """
    Calculate the relative radial velocity of an object
    
    Parameters
    ----------
    x1,y1,z1 : floats
       Coordinattes of object 1 (halo or central)
    x2,y2,z2 : floats
       Coordinattes of object 2 (subhalo or satellite)

    Returns
    -------
    vg : float
       Tangential relative velocity perpendicular to the x-y plane (v gamma)
    """

    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    dvx = vx2-vx1
    dvy = vy2-vy1
    dvz = vz2-vz1

    r = get_r(x1,y1,z1,x2,y2,z2)

    den = r*r*np.sqrt(dx*dx + dy*dy) 
    num = dz*(dx*dvx + dy*dvy) - dvz*(dx*dx + dy*dy)
    vg =  num/den
    return vg



def lbol_coeff(nomeq):
    """
    Return the adequate coefficient to obtain the bolometric luminosity
    as Lbol = coeff*Mdot_BH*c**2
    
    Parameters
    ----------
    nomeq : string
      Indicates the allowed options
      mcc16 for McCarthy+2016 (Table 1 and Sec 7) coefficient 

    Returns
    -------
    coeff : float
       Coefficient in Lbol = coeff*Mdot_BH*c**2
    """
    if (nomeq == 'mcc16'):
        er_mcc16 = 0.1
        ef_mcc16 = 0.15
        return er_mcc16*(1.-ef_mcc16)
    else:
        return None   # default case if neq is not found

    
def get_lbol(lmdotbh,mdot2SI=True,eq='mcc16',units='SI'):
    """
    Obtain the bolometric luminosity given the BH accretion rate

    Parameters
    ----------
    lmdotbh : numpy array of floats
      log10(Mdot_BH/Msun/yr) or log10(Mdot_BH/kg/s)
    mdot2SI : boolean
      True if the input given is log10(Mdot_BH/Msun/yr)
    eq : string
      Specifies the equation to be used (see options in lbol_coeff)
    units : string
      SI for W; cgs for erg/s, sun for L_sun

    Returns
    -------
    lbol : numpy array of floats
        log10(Lbol) with units depending on flags
    """

    if (mdot2SI):
        lmdotbh = lmdotbh + np.log10(Msun_kg) - np.log10(s_in_year)

    coeff = lbol_coeff(eq)
    if (coeff is None): return None
    
    lbol = np.log10(coeff) + lmdotbh + 2*np.log10(c_ms)
    if (units == 'cgs'):
        lbol = lbol + 7.
    elif (units == 'sun'):
        lbol = lbol - np.log10(Lsun_W)
        
    return lbol
    

def get_lLEdd(lmbh):
    """
    Obtain the Eddington luminosity limit, 
    following Eq. 9 and 10 in Griffin+2020 (The evolution ...)

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
       BH mass as log10(M_BH/Msun) or log10(M_BH/Msun/h0).
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
    print(get_r(0,0,1,0,1,0))
    print(get_vg(1,0,0,0,1,0,0,0,0,0,0,0))
    exit()
    lmbh = np.array([7.,8.])
    lmdotbh = lmbh - 2.
    print(lbol_coeff(1),lbol_coeff(34))
    print(get_lbol(lmdotbh),get_lbol(lmdotbh,mdot2SI=False),
          get_lbol(lmdotbh,units='cgs'),get_lbol(lmdotbh,units='sun'))
    print(get_lLEdd(lmbh))
    print(get_lmdotEdd(lmbh,SI=True,h0=0.73))
    
