import sys,os.path
import numpy as np
import h5py
import glob
from iotools import stop_if_no_file, is_sorted
#print('\n \n')

ptypes = ['gas','DM','bp1','bp2','star','BH']

dirbahamasari = '/hpcdata0/simulations/BAHAMAS/'
dirobsari = '/hpcdata0/Obs_Data/'

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
        sys.exit()
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


def get_dirb(env):
    """
    Get the Bahamas directory given the environment

    Parameters
    -----------
    env : string
        ari or cosma, to use the adecuate paths
 
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
    elif (env == 'cosma'):
        dirb = dirbahamascosma

    return dirb


def get_path2data(sim,env):
    """
    Get the path to the directory with the data

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
 
    Returns
    -----
    path2data : string
       Path to data

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_path2data('L050N256/WMAP9/Sims/ex','cosma')
    >>> b.get_path2data('AGN_TUNED_nu0_L100N256_WMAP9','ari')
    """

    # Simulation input
    if (env == 'ari'):
        path2data = dirbahamasari+sim+'/Data/EagleSubGroups_5r200/'
    elif (env == 'cosma'):
        path2data = dirbahamascosma+sim+'/data/'

    return path2data    


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
        sys.exit()

    root = path+'eagle_subfind_tab_'+str(snap).zfill(n0)
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
        ari or cosma, to use the adecuate paths

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
        sys.exit()

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
        sys.exit()
    
    return cenids
    

if __name__== "__main__":
    env = 'ari'

    if (env == 'ari'):
        sim = 'AGN_TUNED_nu0_L100N256_WMAP9'
        dirz = '/hpcdata3/arivgonz/bahamas/'

        print(get_z(-1,sim,dirz,env))
        print(get_z(26,sim,dirz,env))

        snap,zsnap = get_snap(3.2,2.8,3.8,sim,dirz,env)
        print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
        print(get_snap(-100.,-200.,-5.,sim,dirz,env))
        print(get_snap(0.28,0.26,0.3,sim,dirz,env))

    print(get_zminmaxs([0.]))
    print(get_zminmaxs([0.,1.],dz=0.5))
    print(get_zminmaxs([0.,1.]))
    print(get_zminmaxs([0.,1.,0.5]))

