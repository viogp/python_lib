# Dust
import numpy as np

Lsolar_g = 0.02
extinction_curve_g = '/cosma5/data/durham/violeta/gaea/dust_attenuation/al_avcalz.dat'

def gaea_calzetti_lum(llum,leff,av):
    '''
    Function that obtains the dust attenuated log10(luminosity)
    at an effective wavelength, leff, assuming Av,
    using the extinction curves provided for GAEA
    by Michaela Hirschmann.

    Arguments:
    llum : float, intrinsic log10(luminosity)
    leff : float, wavelength at which to calculate the attenuation
    av : float, Av

    Return:
    llum_att : float, dust attenuated log10(luminosity)
    '''

    # Read the tables with the normalised extinction curves:
    # wavelenght (\AA), ec (wavelength)
    wl, ec_wl = np.loadtxt(extinction_curve_g,unpack=True)

    # Interpolate the extinction curve at 5500 \AA)
    ec5500 = np.interp(5500.,wl, ec_wl)

    # Interpolate the extinction curve at wl
    ecl = np.interp(leff, wl, ec_wl)

    # log10(Lobs) = log10(Lint) - 0.4*Av*ec(interpolated at leff)/ec(interpolated at 5500))
    llum_att = llum - 0.4*av*ecl/ec5500

    return llum_att

def gaea_calzetti_mag(mag,leff,av):
    '''
    Function that obtains the dust attenuated magnitud
    at an effective wavelength, leff, assuming Av,
    using the extinction curves provided for GAEA
    by Michaela Hirschmann.

    Arguments:
    mag : float, intrinsic magnitude
    wl : float, wavelength at which to calculate the attenuation
    av: float, Av

    Return:
    mag_att: float, dust attenuated magnitude
    '''

    # Read the tables with the normalised extinction curves:
    # wavelenght (\AA), ec (wavelength)
    wl, ec_wl = np.loadtxt(extinction_curve_g,unpack=True)

    # Interpolate the extinction curve at 5500 \AA)
    ec5500 = np.interp(5500.,wl, ec_wl)

    # Interpolate the extinction curve at wl
    ecl = np.interp(leff, wl, ec_wl)

    # mag_obs = mag_int +  Av * ec(interpolated at leff)/ec(interpolated at 5500))
    mag_att = mag + av*ecl/ec5500

    return mag_att

if __name__ == '__main__':
    leff = 4000.
    av = 1.

    mag_att = gaea_calzetti_mag(24.,leff,av)
    print(mag_att, type(mag_att))
    llum_att = gaea_calzetti_lum(42.,leff,av)
    print(llum_att, type(llum_att))
