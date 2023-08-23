"""
Generic functions dealing with cosmological simulations
"""

import numpy as np


def boundary_correction(xin,box,groups=False):
    """
    Correct an array for periodic boundary conditions

    Parameters
    ----------
    xin : numpy array of floats
       Array with 1-D coordinates.
    box: float
       Size of the boundary to be used.
    groups: boolean
       True for distances within haloes

    Returns
    -------
    xout : numpy array of floats
       Array with the corrected 1-D coordinates.

    Example
    -------
    >>> from cosmosim import boundary_correction; import numpy as np
    >>> boundary_correction(np.array([110.,2,-0.5,100.]),100.)
    >>> [10.   2.  99.5  0. ]
    >>> boundary_correction(np.array([110.,2,-0.5,100.,-49,-51,49,51]),100.,groups=True)
    >>> [10.   2.  -0.5  0. -49 49 49 -49]
    """

    xout = xin.copy()

    if groups:
        lbox2 = box/2.
        xout[xin < -lbox2] = xout[xin < -lbox2] + box
        xout[xin >= lbox2] = xout[xin >= lbox2] - box
    else:
        xout[xin<0] = xout[xin<0] + box
        xout[xin>=box] = xout[xin>=box] - box

    return xout


def get_diffpos(x1,y1,z1,x2,y2,z2,box=None):
    """
    Calculate the relative positions
    
    Parameters
    ----------
    x1,y1,z1 : numpy array of floats
       If groups, coordinates of subhalo or satellites
    x2,y2,z2 : numpy array of floats
       If groups, oordinattes of haloes or centrals
    box : float
       If a simulation, side of the simulation box

    Returns
    -------
    dx,dy,dz : numpy arrays of floats
       Relative coordinates

    Example
    -------
    >>> from cosmosim import get_diffpos; import numpy as np
    >>> get_diffpos(np.array([2.]),np.array([95.]),np.array([2.]),
                    np.array([5.]),np.array([5.]),np.array([80.]),100.)
    >>> (array([-3.]), array([-10.]), array([22.]))
    """

    if box:
        dx = boundary_correction(x1-x2, box, groups=True)
        dy = boundary_correction(y1-y2, box, groups=True)
        dz = boundary_correction(z1-z2, box, groups=True)
    else:
        dx = x1-x2
        dy = y1-y2
        dz = z1-z2

    return dx,dy,dz


def get_r(x1,y1,z1,x2,y2,z2,box=None):
    """
    Calculate the distance between objects
    
    Parameters
    ----------
    x1,y1,z1 : array of floats
       If groups, coordinates of subhalo or satellites
    x2,y2,z2 : array of floats
       If groups, oordinattes of haloes or centrals
    box : float
       If a simulation, side of the simulation box

    Returns
    -------
    r : float
       Distnace between the objects


    Example
    -------
    >>> from cosmosim import get_r; import numpy as np
    >>> get_r(np.array([2.]),np.array([95.]),np.array([2.]),np.array([5.]),np.array([5.]),np.array([80.]),100.,groups=True)
    >>> (array([]))
    """

    dx,dy,dz = get_diffpos(x1,y1,z1,x2,y2,z2,box)
    return np.sqrt(dx*dx + dy*dy + dz*dz)


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




if __name__ == "__main__":
    print(boundary_correction(np.array([110.,2,-0.5,100.]),100.))
    print(boundary_correction(np.array([110.,2,-0.5,100.,-49,-51,49,51]),100.,groups=True))
    print(get_diffpos(np.array([2.]),np.array([95.]),np.array([2.]),np.array([5.]),np.array([5.]),np.array([80.]))) 
    print(get_diffpos(np.array([2.]),np.array([95.]),np.array([2.]),np.array([5.]),np.array([5.]),np.array([80.]),100.))
    print(get_r(np.array([2.]),np.array([95.]),np.array([2.]),np.array([5.]),np.array([5.]),np.array([80.]),100.)) 
    print(get_r(0,0,1,0,1,0),100.)
    #print(get_vg(1,0,0,0,1,0,0,0,0,0,0,0))

    
