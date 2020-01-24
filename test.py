#! /usr/bin/env python

import Cosmology as c
import numpy as np
from astropy.cosmology import Planck15 as cosmo

print cosmo

area = 100

#--------------------
z=1.03
c.set_Planck15() 
l = 10.**39
lum = c.h*c.h*l/(10.**40.) #; print lum
print ("F(erg/s/cm^2 for L=%e erg/s, z=%f)=%e" % (l,z,c.emission_line_flux(lum,z)))

c.set_MR7()
print ("F(erg/s/cm^2 for L=%e erg/s, z=%f)=%e" % (l,z,c.emission_line_flux(lum,z)))

c.set_Planck15() 
flux_data = 2.7*1./10.**17.
print ("L(erg/s for F=%e erg/s/cm^2, z=%f)=%e" % (flux_data,z,c.emission_line_luminosity(flux_data,z)*10.**40./(c.h*c.h)))



#----------------------------------------
h0 =0.704 
sarea = 4.*np.pi*(180./np.pi)**2
vol = c.comoving_volume(1.6)*h0**3*area/sarea

#print vol,np.power(vol,1./3.)

z = 0.
bias = [1., 1., 1.]
print c.kaiser_factor(z,bias) 

z=1.
om = c.omegam(z)
sigma8 = 0.52*om**(-0.52+0.13*om) # Eke, Cole & Frenk 1996
bias = 1.7*0.81/sigma8
print om,sigma8
print bias,0.53+0.29*(1.+z)**2


