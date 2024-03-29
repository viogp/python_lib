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
    
    if hasattr(xin, "__len__"):
        xout = xin.copy()

        if groups:
            lbox2 = box/2.
            xout[xin < -lbox2] = xout[xin < -lbox2] + box
            xout[xin >= lbox2] = xout[xin >= lbox2] - box
        else:
            xout[xin<0] = xout[xin<0] + box
            xout[xin>=box] = xout[xin>=box] - box
    else:
        if groups:
            lbox2 = box/2.
            if (xin < -lbox2):
                xout = xin + box
            elif (xin >= lbox2):
                xout = xin - box
            else:
                xout = xin
        else:
            if (xin < 0):
                xout = xin + box
            elif (xin >= box):
                xout = xin - box
            else:
                xout = xin
            
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
    x1,y1,z1 : (array of) floats
       Coordinattes of object 1 (halo or central). Array for simulations.
    x2,y2,z2 : (array of) floats
       Coordinattes of object 2 (subhalo or satellite). Array for simulations.
    box : float
       If a simulation, side of the simulation box

    Returns
    -------
    r : float
       Distnace between the objects


    Example
    -------
    >>> from cosmosim import get_r; import numpy as np
    >>> get_r(np.array([2.,3.]),np.array([95.,5.]),np.array([2.,0.]),np.array([5.,-5.]),np.array([5.,50.]),np.array([80.,90.]),100.)
    >>> array([24.35159132, 46.78675026])
    """

    dx,dy,dz = get_diffpos(x1,y1,z1,x2,y2,z2,box)

    dr = np.sqrt(dx*dx + dy*dy + dz*dz)
    return dr


def get_vr(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,box=None):
    """
    Calculate the relative radial velocity of an object
    
    Parameters
    ----------
    x1,y1,z1 : (array of) floats
       Coordinattes of object 1 (halo or central). Array for simulations.
    x2,y2,z2 : (array of) floats
       Coordinattes of object 2 (subhalo or satellite). Array for simulations.
    vx1,vy1,vz1 : (array of) floats
       Velocity of object 1 (halo or central). Array for simulations.
    vx2,vy2,vz2 : (array of) floats
       Velocity of object 2 (subhalo or satellite). Array for simulations.
    box : float
       If a simulation, side of the simulation box

    Returns
    -------
    vr : float
       Radial relative velocity
    """

    r = get_r(x1,y1,z1,x2,y2,z2,box)

    dx,dy,dz = get_diffpos(x1,y1,z1,x2,y2,z2,box)
    dvx,dvy,dvz = get_diffpos(vx1,vy1,vz1,vx2,vy2,vz2)

    if hasattr(r, "__len__"):
        vr = np.empty(len(r)); vr[:] = np.nan    
        ind = np.where(r>0)
        if (np.shape(ind)[1]>0):
            vr[ind] = (dx[ind]*dvx[ind] + dy[ind]*dvy[ind] + dz[ind]*dvz[ind])/r[ind]
    else:
        vr = np.nan
        if (abs(r)>0): vr = (dx*dvx + dy*dvy + dz*dvz)/r

    return vr


def get_vtheta(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,box=None):
    """
    Calculate the relative tangential velocity in the x-y plane (v theta).
    
    Parameters
    ----------
    x1,y1,z1 : (array of) floats
       Coordinattes of object 1 (halo or central). Array for simulations.
    x2,y2,z2 : (array of) floats
       Coordinattes of object 2 (subhalo or satellite). Array for simulations.
    vx1,vy1,vz1 : (array of) floats
       Velocity of object 1 (halo or central). Array for simulations.
    vx2,vy2,vz2 : (array of) floats
       Velocity of object 2 (subhalo or satellite). Array for simulations.
    box : float
       If a simulation, side of the simulation box

    Returns
    -------
    vtheta : float
       Relative tangential velocity in the x-y plane (v theta)
    """

    dx,dy,dz = get_diffpos(x1,y1,z1,x2,y2,z2,box)
    dvx,dvy,dvz = get_diffpos(vx1,vy1,vz1,vx2,vy2,vz2)

    den = np.sqrt(dx*dx + dy*dy) 

    if hasattr(den, "__len__"):
        vtheta = np.empty(len(den)); vtheta[:] = np.nan    
        ind = np.where(den>0)
        if (np.shape(ind)[1]>0):
            vtheta[ind] =  (dx[ind]*dvy[ind] - dy[ind]*dvx[ind])/den[ind]
    else:
        vtheta = np.nan
        if (abs(den)>0): vtheta =  (dx*dvy - dy*dvx)/den
        
    return vtheta


def get_vphi(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,box=None):
    """
    Calculate the relative tangential velocity (perpendicular to the x-y plane) of an object
    
    Parameters
    ----------
    x1,y1,z1 : (array of) floats
       Coordinattes of object 1 (halo or central). Array for simulations.
    x2,y2,z2 : (array of) floats
       Coordinattes of object 2 (subhalo or satellite). Array for simulations.
    vx1,vy1,vz1 : (array of) floats
       Velocity of object 1 (halo or central). Array for simulations.
    vx2,vy2,vz2 : (array of) floats
       Velocity of object 2 (subhalo or satellite). Array for simulations.
    box : float
       If a simulation, side of the simulation box

    Returns
    -------
    vphi : float
       Relative tangential velocity perpendicular to the x-y plane (v phi)
    """
    r = get_r(x1,y1,z1,x2,y2,z2,box)

    dx,dy,dz = get_diffpos(x1,y1,z1,x2,y2,z2,box)
    dvx,dvy,dvz = get_diffpos(vx1,vy1,vz1,vx2,vy2,vz2)

    num = dz*(dx*dvx + dy*dvy) - dvz*(dx*dx + dy*dy)
    den = r*r*np.sqrt(dx*dx + dy*dy)

    if hasattr(r, "__len__"):
        vphi = np.empty(len(r)); vphi[:] = np.nan    
        ind = np.where(den>0)
        if (np.shape(ind)[1]>0):
            vphi[ind] = num[ind]/den[ind]
    else:
        vphi = np.nan
        if (abs(den)>0): vphi = num/den
        
    return vphi



def get_vlos(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,box=None):
    """
    Calculate the relative velocity on the line of sight (assumed z-axis)
    
    Parameters
    ----------
    x1,y1,z1 : (array of) floats
       Coordinattes of object 1 (halo or central). Array for simulations.
    x2,y2,z2 : (array of) floats
       Coordinattes of object 2 (subhalo or satellite). Array for simulations.
    vx1,vy1,vz1 : (array of) floats
       Velocity of object 1 (halo or central). Array for simulations.
    vx2,vy2,vz2 : (array of) floats
       Velocity of object 2 (subhalo or satellite). Array for simulations.
    box : float
       If a simulation, side of the simulation box

    Returns
    -------
    vlos : float
       Relative velocity on the line of sight (z-axis)
    """

    vr   =   get_vr(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,box)
    vphi = get_vphi(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,box)

    dx,dy,dz = get_diffpos(x1,y1,z1,x2,y2,z2,box)    
    r = get_r(x1,y1,z1,x2,y2,z2,box)

    if hasattr(r, "__len__"):                                                                                
        vlos = np.empty(len(r)); vlos[:] = np.nan
        cosphi = np.empty(len(r)); cosphi[:] = np.nan
        ind = np.where(r>0)
        if (np.shape(ind)[1]>0):
            cosphi[ind] = dz[ind]/r[ind]

            sinphi = np.sqrt(1-cosphi*cosphi)

            vlos = vr*cosphi - vphi*sinphi
    else:                                                                                                    
        vlos = np.nan
        if (abs(r)>0):
            cosphi = dz/r
            sinphi = np.sqrt(1-cosphi*cosphi)
            vlos = vr*cosphi - vphi*sinphi

    return vlos


if __name__ == "__main__":
    x1 = np.array([2.])
    y1 = np.array([95.])
    z1 = np.array([2.])

    x2 = np.array([5.])
    y2 = np.array([5.])
    z2 = np.array([80.])

    vx1 = np.array([-3.])
    vy1 = np.array([-10.])
    vz1 = np.array([22.])

    vx2 = np.array([0.])
    vy2 = np.array([0.])
    vz2 = np.array([0.])
        
    print(boundary_correction(110.,100.))
    print(boundary_correction(2,100.,groups=True))
    print(boundary_correction(np.array([110.,2,-0.5,100.]),100.))
    print(boundary_correction(np.array([110.,2,-0.5,100.,-49,-51,49,51]),100.,groups=True))
    print(get_diffpos(x1,y1,z1,x2,y2,z2)) 
    print(get_diffpos(x1,y1,z1,x2,y2,z2,100.))
    print('r(0,0,1,0,1,0)={}'.format(get_r(0,0,1,0,1,0)))
    print('vr(1,0,0,0,1,0,0,0,0,0,0,0)={}'.format(get_vr(1,0,0,0,1,0,0,0,0,0,0,0)))
    print('vtheta(5,3,9,2,0,9,2,3,9,0,1,9)={}'.format(get_vtheta(5,3,9,2,0,9,2,3,9,0,1,9)))
    print('vphi(1,2,0,3,4,0,1,2,0,3,4,0)={}'.format(get_vphi(1,2,0,3,4,0,1,2,0,3,4,0)))
    print('vlos(1,2,9,3,4,8,1,2,9,3,4,9)={}'.format(get_vlos(1,2,9,3,4,8,1,2,9,3,4,9)))
    print('r(box)={}'.format(get_r(x1,y1,z1,x2,y2,z2,100.)))
    print('vr(box)={}'.format(get_vr(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,100.)))
    print('vtheta(box)={}'.format(get_vtheta(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,100.)))
    print('vphi(box)={}'.format(get_vphi(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,100.)))
    print('vlos(box)={}'.format(get_vlos(x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,100.)))
