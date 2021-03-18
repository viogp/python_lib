import os.path
import numpy as np
import h5py
import glob
import subprocess
from astropy import constants as const
from iotools import stop_if_no_file, is_sorted
#print('\n \n')

ptypes = ['gas','DM','bp1','bp2','star','BH']

dirbahamasarilega = '/hpcdata0/simulations/BAHAMAS/'
dirbahamasari = '/hpcdata0/arivgonz/BAHAMAS/'
dirobsari = '/hpcdata0/Obs_Data/'

dirbahamascosmalega = '/cosma6/data/dp004/Eagle/jsTestRuns/BAHAMAS_XL/'
dirbahamascosma = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'
dirobscosma = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/Obs_Data/'

tblz = 'snap_z.txt'

defaultdz = 0.25

n0 = 3


def get_zminmaxs(zz,dz=None):
    """
    Get the previous (min) and next (max) values
    given an array. If the input is a single value 
    zmaxs = zz+dz, zmins = zz-dz

    Parameters
    -----------
    zz : list of floats
        Redshift 
    dz : float
        Optional parameter with an interval

    Returns
    -----
    zmins : list of floats
        Minimum redshifts
    zmaxs : list of floats
        Maximum redshifts

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_zminmaxs([0.,1.])
    >>> [-1.0, 0.0], [1.0, 2.0]
    """

    # Check that the input list is sorted
    if(zz != sorted(zz)):
        print('STOP (get_zminmaxs): Sort the redshift list')
        exit()
        return -999.,-999.
    
    if (dz is not None):
        zmins = [iz - dz for iz in zz]
        zmaxs = [iz + dz for iz in zz]
    elif(len(zz)<2):
        print('WARNING (get_zminmaxs): Setting to default dz={}'.format(defaultdz))
        dz = defaultdz
        zmins = [iz - dz for iz in zz]
        zmaxs = [iz + dz for iz in zz]
    else:
        zmins = zz[:-1]
        zmins.insert(0,2*zz[0]-zz[1])

        zmaxs = zz[1:]
        zmaxs.append(2*zz[-1]-zz[-2]) 
        
    return zmins,zmaxs



def mb2msun(massb,h0):
    """
    Provide bahamas masses as M/Msun

    Parameters
    -----------
    massb : float
        Mass in 10^10Msun/h units
    h0 : float
        Hubble constant

    Returns
    -----
    mass : float
        M/Msun

    Examples
    ---------
    >>> import bahamas as b
    >>> b.mb2msun(0.0048,0.7)
    >>> 68571428.57142857
    """

    if(massb <= 0):
        print('WARNING (get_mb2msun): Input mass <= 0, returning -999.')
        return -999.
    else:
        mass = massb*np.power(10.,10.)/h0
        return mass


def mb2lmsun(massb,h0):
    """
    Provide bahamas masses as log10(M/Msun)

    Parameters
    -----------
    massb : float
        Mass in 10^10Msun/h units
    h0 : float
        Hubble constant

    Returns
    -----
    lmass : float
        log10(M/Msun)

    Examples
    ---------
    >>> import bahamas as b
    >>> b.mb2lmsun(0.0048,0.7)
    >>> 7.836143197361331
    """

    if(massb <= 0):
        print('WARNING (get_mb2lmsun): Input mass <= 0, returning -999.')
        return -999.
    else:
        lmass = np.log10(massb) + 10. - np.log10(h0)
        return lmass

    
def get_dirb(env):
    """
    Get the Bahamas directory given the environment

    Parameters
    -----------
    env : string
        cosma, cosmalega, ari or arilega, to use the adecuate paths
 
    Returns
    -----
    dirb : string
       Bahamas directory

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_dirb('cosma')
    """

    dirb = None
    if (env == 'ari'):
        dirb = dirbahamasari
    elif (env == 'arilega'):
        dirb = dirbahamasarilega
    elif (env == 'cosma'):
        dirb = dirbahamascosma
    elif (env == 'cosmalega'):
        dirb = dirbahamascosmalega
    else:
        print('WARNING (b.get_dirb): environment name not cosma, cosmalega, ari or arilega')

    return dirb

def get_dirobs(env):
    """
    Get the directory with observations given the environment

    Parameters
    -----------
    env : string
        cosma, cosmalega, ari or arilega, to use the adecuate paths
 
    Returns
    -----
    dirb : string
       Bahamas directory

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_dirobs('cosma')
    """

    dirobs = None
    if (env == 'ari' or env == 'arilega'):
        dirobs = dirobsari
    elif (env == 'cosma' or env == 'cosmalega'):
        dirobs = dirobscosma

    return dirobs


def get_path2data(sim,env):
    """
    Get the path to the directory with the data

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    env : string
        cosma, cosmalega, ari or arilega, to use the adecuate paths
 
    Returns
    -----
    path2data : string
       Path to data

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_path2data('L050N256/WMAP9/Sims/ex','cosma')
    >>> b.get_path2data('AGN_TUNED_nu0_L100N256_WMAP9','arilega')
    """

    # Simulation input
    if (env == 'ari'):
        path2data = dirbahamasari+sim+'/data/'
    elif (env == 'arilega'):
        path2data = dirbahamasarilega+sim+'/Data/EagleSubGroups_5r200/'
    elif (env == 'cosma'):
        path2data = dirbahamascosma+sim+'/data/'
    elif (env == 'cosmalega'):
        path2data = dirbahamascosmalega+sim+'/data/'
    else:
        exit('get_path2data set to handle env=cosma, cosmalega, ari or arilega')

    return path2data    


def get_particle_files(snap,sim,env):
    """
    Get the particle files

    Parameters
    -----------
    snap : integer
        Snapshot number
    sims : list of strings
        Array with the names of the simulation
    env : string
        ari(lega) or cosma(lega), to use the adecuate paths
 
    Returns
    -----
    files : array of string
       Subfind files with full paths

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_particle_files(8,'L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path1 = get_path2data(sim,env)+'particledata_'+str(snap).zfill(n0)

    # Get path to subfind files
    paths = glob.glob(path1+'*/') 
    if (len(paths) == 1):
        path = paths[0]
    else:
        print('STOP(~/python_lib/bahamas): more than one or none directories with root {}'.format(path1+'*/'))
        exit()

    root = path+'eagle_subfind_particles_'+str(snap).zfill(n0) 
    files = glob.glob(root+'*.hdf5')

    return files     


def get_subfind_files(snap,sim,env):
    """
    Get the subfind files

    Parameters
    -----------
    snap : integer
        Snapshot number
    sims : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
 
    Returns
    -----
    files : array of string
       Subfind files with full paths

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_subfind_files(8,'L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path1 = get_path2data(sim,env)+'groups_'+str(snap).zfill(n0)

    # Get path to subfind files
    paths = glob.glob(path1+'*/')
    if (len(paths) == 1):
        path = paths[0]
    else:
        print('STOP(~/python_lib/bahamas): more than one or none directories with root {}'.format(path1+'*/'))
        exit()

    root = path+'eagle_subfind_tab_'+str(snap).zfill(n0)
    files = glob.glob(root+'*.hdf5')

    return files     


def get_particle_files(snap,sim,env): 
    """
    Get the particle files

    Parameters
    -----------
    snap : integer
        Snapshot number
    sims : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
 
    Returns
    -----
    files : array of string
       Subfind files with full paths

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_particle_files(8,'L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path1 = get_path2data(sim,env)+'particledata_'+str(snap).zfill(n0)

    # Get path to the particle files
    paths = glob.glob(path1+'*/')
    if (len(paths) == 1):
        path = paths[0]
    else:
        print('STOP(~/python_lib/bahamas): more than one or none directories with root {}'.format(path1+'*/'))
        exit()

    root = path+'eagle_subfind_particles_'+str(snap).zfill(n0)
    files = glob.glob(root+'*.hdf5')

    return files     


def get_cosmology(sim,env):
    """
    Get the cosmology for a simulation

    Parameters
    -----------
    sim : string
        Name of the Bahamas directory.
    env : string
        ari or cosma, to use the adecuate paths.

    Returns
    -----
    omega0, omegab, lambda0, h0 : floats
        Cosmological parameters

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_cosmology('AGN_TUNED_nu0_L100N256_WMAP9','ari')
    >>> b.get_cosmology('L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path = get_path2data(sim,env)

    # Initialize arrays for z and sn
    files = glob.glob(path+'groups_*/group_tab*')
    infile = files[0]
    f = h5py.File(infile, 'r')
    header = f['Header']
    #print(list(header.attrs.items()))

    omega0 = header.attrs['Omega0']
    omegab = header.attrs['OmegaBaryon']
    lambda0 = header.attrs['OmegaLambda']
    h0 = header.attrs['HubbleParam']

    return omega0, omegab, lambda0, h0


def table_z_sn(sim,env,dirz=None):
    """
    Produce a table with redshifts and snapshot numbers
    for a simulation.

    Parameters
    -----------
    sim : string
        Name of the Bahamas directory.
    env : string
        ari or cosma, to use the adecuate paths.
    dirz : string
        Alternative path to table with z and snapshot.

    Returns
    -----
    tablez : string
        Full path to file containing the list of 
        redshifts and snapshots for this simulation.

    Examples
    ---------
    >>> import bahamas as b
    >>> b.table_z_sn('AGN_TUNED_nu0_L100N256_WMAP9','ari')
    >>> b.table_z_sn('L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path = get_path2data(sim,env)

    # Output file
    if (dirz == None):
        tablez = path+tblz
    else:
        dirz = dirz+sim+'/'
        if (not os.path.exists(dirz)):
            os.makedirs(dirz)
        tablez = dirz+tblz

    # Initialize arrays for z and sn
    dirs = glob.glob(path+'groups_0*')
    if (len(dirs) < 1) :
        print('STOP (bahamas.table_z_sn): {} not containing expected files'.format(path))
        exit()
        
    zzs = np.zeros(shape=len(dirs)) ; zzs.fill(-999.)
    sns = np.zeros(shape=len(dirs), dtype=int)

    for ii, dir in enumerate(dirs):
        ending = dir.split('groups_')[1]
        if ('_z' in ending):
            snap = ending.split('_z')[0]
        else:
            snap = ending
        sns[ii] = int(snap)

        infile = dir+'/group_tab_'+ending+'.0.hdf5'
        if (not os.path.isfile(infile)):
            print('WARNING: Files missing in {}'.format(dir))
            continue

        f = h5py.File(infile, 'r')
        header = f['Header']

        zzs[ii] = header.attrs['Redshift']
        #zz = '{0:.2f}'.format(header.attrs['Redshift'])
        f.close()

    # Sort by redshift
    ind = np.argsort(zzs)
    zz = zzs[ind]
    sn = sns[ind]

    # Write list to file
    with open(tablez, 'w') as f:
        f.write('# Redshift Snapshot \n')
        for ii in range(len(zz)):
            tofile = '{:.2f} {:d}'.format(zz[ii],sn[ii])
            f.write("%s\n" % tofile)

    return tablez


def get_z(snap,sim,env,dirz=None):
    """
    Get the redshift for a snapshot number and
    a simulation name.

    Parameters
    -----------
    snap : int
        Snapshot number
    sim : string
        Name of the simulation
    env : string
        ari or cosma, to use the adecuate paths
    dirz : string
        Alternative path to table with z and snapshot.

    Returns
    -----
    snap : string
        Snapshot value

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_z(26,'AGN_TUNED_nu0_L100N256_WMAP9','ari',dirz='/hpcdata3/arivgonz/bahamas/')
    >>> b.get_z(8,'L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path = get_path2data(sim,env)

    # Table with increasing redshifts and corresponding snapshots
    if (dirz == None):
        tablez = path+tblz
    else:
        tablez = dirz+sim+'/'+tblz

    if (not os.path.isfile(tablez)):
        # Generate the table if it doesn't exist
        tablez = table_z_sn(sim,env,dirz=dirz)

    # Read the table:
    zzs, snsf = np.loadtxt(tablez, unpack=True)

    if (np.isscalar(zzs)):
        # Case of scalars
        sns = int(snsf)
        if (sns == snap):
            return zzs
        else:
            print('WARNING: no snapshot number {} found in {}'.format(snap,tablez))
            return -999.            
    else:
        sns = np.asarray(snsf,dtype=int)
        # Find the index of the snapshot number
        ind = np.where(sns == snap)
        if (np.shape(ind)[1]<1):
            print('WARNING: no snapshot number {} found in {}'.format(snap,tablez))
            return -999.
        elif (np.shape(ind)[1]>1):
            print('WARNING: several snapshot number {} found in {}'.format(snap,tablez))
            return -999.
        else:
            zz = zzs[ind][0]
            return zz


def get_snap(zz,zmin,zmax,sim,env,dirz=None):
    """
    Get the closest snapshot given a redshift and
    a simulation name, within some limits

    Parameters
    -----------
    zz : float
        Redshift
    zmin : float
        Minimum redsfhit to look for the sanpshot
    zmax : float
        Maximum redsfhit to look for the sanpshot
    sim : string
        Name of the Bahamas simulation
    dirz : string
        Alternative directory where the table with z and snapshots is.
    env : string
        ari, arilega or cosma, cosmalega to use the adecuate paths

    Returns
    -----
    snap : int
        Snapshot value
    z_snap : float
        Redshift for that snapshot

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_snap(3.2,2.8,3.6,'AGN_TUNED_nu0_L100N256_WMAP9','ari')
    >>> b.get_snap(12.3,12.0,12.6,'L050N256/WMAP9/Sims/ex','cosma')
    >>> (2, 12.5)
    >>> snap, z_snap = b.get_snap(99.,20.,150.,'L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path = get_path2data(sim,env)

    # Table with increasing redshifts and corresponding snapshots
    if (dirz == None):
        tablez = path+tblz
    else:
        tablez = dirz+sim+'/'+tblz

    if (not os.path.isfile(tablez)):
        # Generate the table if it doesn't exist
        tablez = table_z_sn(sim,env,dirz=dirz)

    # Read the table:
    zzs, lsns = np.loadtxt(tablez, unpack=True)        
            
    # Case of having a single snapshot:
    if (not hasattr(zzs, "__len__")):
        if (zmin < zzs and zzs < zmax):
            return int(lsns),float(zzs)
        else:
            print('WARNING: {} far from single {} found in {}'.format(zz,zzs,tablez))
            return -999,-999.
    sns = np.asarray(lsns,dtype=int)

    # Find the closest redshift
    if ((zz < zzs[0] and zmax < zzs[0]) or (zz<0.)):
        print('WARNING: {}<{}, min. z found in {}'.format(zz,zzs[0],tablez))
        return -999,-999.
    elif(zz > zzs[-1] and zmin > zzs[-1]):
        print('WARNING: {}>{}, max. z found in {}'.format(zz,zzs[-1],tablez))
        return -999.,-999.
    else:
        idx = (np.abs(zzs - zz)).argmin()
        if (zmin < zzs[idx] and zzs[idx] < zmax):
            return sns[idx],zzs[idx]
        else:
            print('WARNING: z={} outside range {}<z<{}, {}'.format(zzs[idx],zmin,zmax,tablez))
            return -999,-999.


def cenids(snap,sim,env):
    """
    Get the list of indexes for central galaxies

    Parameters
    -----------
    snap : integer
        Snapshot 
    sim : string
        Simulation name
    env : string
        ari or cosma, to use the adecuate paths

    Returns
    -----
    cenids : numpy array int
        Global indexes for central galaxies

    Examples
    ---------
    >>> import bahamas as b
    >>> b.cenids(31,'HIRES/AGN_TUNED_nu0_L050N256_WMAP9','ari')
    >>> b.cenids(8,'L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    files = get_subfind_files(snap,sim,env)

    # Cycle through the files
    lenf = len(files)
    if (lenf<1):
        print('STOP(~/python_lib/bahamas): Make sure you can see the path {}'.format(path))
        exit()

    for ii,ff in enumerate(files):
        stop_if_no_file(ff)
        f = h5py.File(ff, 'r')
        haloes = f['FOF']
        if (ii == 0):
            cenids = haloes['FirstSubhaloID'][:]
        else:
            np.append(cenids, haloes['FirstSubhaloID'][:]) 

    # Check that the array is sorted
    sorted = is_sorted(cenids)
    if (not sorted):
        print('STOP(bahamas/cenids): central indexes not sorted')
        exit()
    
    return cenids


def resolution(sim,env,zz=0.,verbose=True):
    """
    Get the mass resolution of a simulation

    Parameters
    -----------
    sim : string
        Simulation name
    env : string
        ari or cosma, to use the adecuate paths
    zz : float
        Redshift to look (should be the same for all)
    verbose : boolean
        True to print the resolution

    Returns
    -----
    mdm, mgas : float
        Mass resolution (Msun) for DM and gas particles

    Examples
    ---------
    >>> import bahamas as b
    >>> b.resolution('HIRES/AGN_TUNED_nu0_L050N256_WMAP9','arilega')
    >>> b.resolution('L050N256/WMAP9/Sims/ex','cosma')
    """

    snap, z_snap = get_snap(zz,-999.,999.,sim,env)
    
    # Simulation input
    files = get_particle_files(snap,sim,env)

    f= h5py.File(files[0],'r')
    header=f['Header']
    masstable = header.attrs['MassTable']

    itype = ptypes.index('DM') 
    mdm = masstable[itype]
    if (mdm<0):
        print('WARNING: negative or 0 input mass, returning -999.')
        return -999., -999.
    else:
        omega0, omegab, lambda0, h0 = get_cosmology(sim,env)

        mdm = mb2msun(mdm,h0)

        mgas = mdm*omegab/(omega0-omegab)

        if (verbose):
            print('Particle resolution of sim.: {}'.format(sim))
            print('mdm (Msun) = {:.2e}, mgas (Msun)=mb= {:.2e}'.format(mdm,mgas))
        
        return mdm,mgas


def get_min_sfr(sim,env,zz=0.,A=1.515*1e-4,gamma=5/3,fg=1.,n=1.4,verbose=True):
    '''
    Get the minimum star formation rate, following eq 1 in Schaye et al. 2015

    Parameters
    -----------
    sim : string
        Simulation name
    env : string
        ari or cosma, to use the adecuate paths
    zz : float
        Redshift to read files
    A : float
        Kennicutt-Schmidt's law constant in Msun yr^-1 kpc^-2 units
    gamma : float
        Ratio of specific heats
    fg : float
        Mass fraction in gas
    n : float
        Index of the Kennicutt-Schmidt's law
    verbose : boolean
        True to print the resolution

    Returns
    -----
    minsfr : float
        Minimum theoretical star formation rate in Msun/yr units

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_min_sfr('HIRES/AGN_TUNED_nu0_L050N256_WMAP9','arilega')
    '''

    # Typical values at the outskirts of haloes, where SF happens
    nHsf   = 0.001  # cm^-3  
    Tsf = 100000 # K      

    # Get the particle resolution
    mdm, mgas = resolution(sim,env,verbose=False)

    apc = A/1e6 # Msun yr^-1 pc^-2

    nHm = nHsf * 100.**3 # m^-3
    P = nHm*const.k_B.value*Tsf
    
    gp_si = np.sqrt(gamma*fg*P/const.G.value)
    gp = gp_si*const.pc.value**2/const.M_sun.value

    minsfr = mgas*apc*np.power(gp,n-1)*10**9
    if verbose:
        print('  Min. theoretical log10(SFR (Msun/Gyr)) = {:2f}'.format(np.log10(minsfr)))
    
    return minsfr


def get_nh(zz,massdef,sim,env,mmin=9.,mmax=16.,dm=0.1,outdir=None,Testing=True):
    '''
    Calculate the number of haloes per mass bin and write this into a file

    Parameters
    -----------
    zz : float
        Redshift to get the number of haloes
    massdef : string
        Name of the mass definition to be used
    sim : string
        Name of the simulation
    env : string
        ari, arilega or cosma, to use the adecuate paths
    mmin : float
        Mimimum mass to be considered
    mmax : float
        Maximum mass to be considered
    dm : float
        Intervale step for the halo mass
    outdir : string
        Path to output file
    Testing : boolean
        Calculations on part or all the simulation

    Returns
    -----
    outfile : string
        Path to file with: 

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.get_nh(31,'Group_M_Mean200',sim,'arilega',outdir='/hpcdata0/arivgonz/BAHAMAS/')
    '''

    # Get snapshot
    zmin,zmax = get_zminmaxs([zz])
    snap, z_snap = get_snap(zz,zmin,zmax,sim,env)

    # Output file
    if outdir:
        path = outdir+sim
    else:
        path = get_dirb(env)+sim
    if (not os.path.exists(path)):
        os.makedirs(path)

    if Testing:
        outfil = path+'/nh_'+massdef+'_sn'+str(snap)+'_dm'+str(dm)+'_test.txt'
    else:
        outfil = path+'/nh_'+massdef+'_sn'+str(snap)+'_dm'+str(dm)+'.txt'
    if (os.path.isfile(outfil)):
        nlines = subprocess.call(['wc', '-l', outfil])
        if (nlines > 0):
            return outfil

    # The subfiles to loop over
    nvols = 'All'
    if Testing: nvols = 2

    # Bins in halo mass
    edges = np.array(np.arange(mmin,mmax,dm))
    mhist = edges[1:]-0.5*dm
    nh  = np.zeros(shape=(len(mhist)))

    elow  = edges[:-1]
    ehigh = edges[1:]
    
    # Get subfind files
    files = get_subfind_files(snap,sim,env)
    if (len(files)<1):
        print('WARNING (b.get_nh): no subfind files at snap={}, {} '.format(snap,sim))
        return None
    
    # Loop over the files
    volume = 0.
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r')
        # Read volume in first iteration
        if (iff == 0):
            header = f['Header']
            volume = np.power(header.attrs['BoxSize'],3.)

        haloes = f['FOF']
        mh = haloes[massdef][:]  #10^10Msun/h
        ind = np.where(mh > 0.)
        lmh = np.log10(mh[ind]) + 10. # log10(M/Msun/h)

        # Number of haloes
        H, bins_edges = np.histogram(lmh,bins=edges)
        nh[:] = nh[:] + H

        if (nvols != 'All'):
            if (iff>nvols): break

    # Output only bins with haloes
    indh = np.where(nh > 0) 
    tofile = np.array([mhist[indh], elow[indh], ehigh[indh], nh[indh]])

    # Write to output file
    with open(outfil, 'w') as outf:
        # Header
        outf.write('# '+sim+', snapshot='+str(snap)+' \n')
        outf.write('# Volume='+str(volume)+', dm='+str(dm)+' \n' )
        outf.write('# log(Mh/Msun/h)_midpoint, log(Mh)_low, log(Mh)_high, Number of haloes \n')

        # Data
        np.savetxt(outf,tofile.T,fmt='%.5f %.5f %.5f %.0f')

    return outfil


if __name__== "__main__":
    env = 'ari'
    env = 'cosmalega'

    if (env == 'cosmalega'):
        sim = 'L400N1024/WMAP9/Sims/BAHAMAS'
        dirz = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/')
        print(table_z_sn(sim,'cosmalega',dirz=dirz)
        print(resolution(sim,'cosmalega'))

    if (env == 'ari'):
        print(resolution('AGN_TUNED_nu0_L400N1024_WMAP9','arilega'))
        print(resolution('HIRES/AGN_RECAL_nu0_L100N512_WMAP9','arilega'))

        sim = 'L050N256/WMAP9/Sims/ws_324_23_mu_7_05_dT_8_35_n_75_BH_beta_1_68_msfof_1_93e11'

        #print(get_z(-1,sim,dirz,env))
        #print(get_z(26,sim,dirz,env))
        #
        #snap,zsnap = get_snap(3.2,2.8,3.8,sim,dirz,env)
        #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
        #print(get_snap(-100.,-200.,-5.,sim,dirz,env))
        #print(get_snap(0.28,0.26,0.3,sim,dirz,env))
        #print(resolution(sim,env))
        #print(np.log10(get_min_sfr(sim,env))+9.)

    #print(get_zminmaxs([0.]))
    #print(get_zminmaxs([0.,1.],dz=0.5))
    #print(get_zminmaxs([0.,1.]))
    #print(get_zminmaxs([0.,1.,0.5]))

