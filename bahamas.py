import os.path
import numpy as np
import h5py
import glob
import subprocess
from astropy import constants as const
from iotools import stop_if_no_file, is_sorted, create_dir
#print('\n \n')

ptypes = ['gas','DM','bp1','bp2','star','BH']

unitdefault = {
    'volume': '(Mpc/h)^3',
    'mass': 'Msun/h',
    'sfr': 'to be filled',
    'ssfr': 'to be filled'
}

#dirbahamasarilega = '/hpcdata0/simulations/BAHAMAS/'
#dirbahamasari = '/hpcdata0/arivgonz/BAHAMAS/'
#dirobsari = '/hpcdata0/Obs_Data/'

dirbahamasarilega = '/beegfs2/hpcdata0_backup/simulations/BAHAMAS/' # Output: '/beegfs2/arivgonz/BAHAMAS/'
dirbahamasari = '/enc1/hpcddn/hpcdata3/arivgonz/BAHAMAS/' #from havok
dirobsari = '/beegfs2/Obs_Data/'

dirbahamascosmalega = '/cosma6/data/dp004/Eagle/jsTestRuns/BAHAMAS_XL/'
dirbahamascosma = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'
dirobscosma = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/Obs_Data/'

tblz = 'snap_z.txt'

defaultdz = 0.25

n0 = 3

def print_h5attributes(infile,inhead='Header'):
    """
    Print out the group attributes of a hdf5 file

    Parameters
    ----------
    infile : string
      Name of input file (this should be a hdf5 file)
    inhead : string
      Name of the group to read the attributes from

    Example
    -------
    >>> import bahamas as b
    >>> infile = '/hpcdata0/simulations/BAHAMAS/AGN_TUNED_nu0_L100N256_WMAP9/Data/Snapshots/snapshot_026/snap_026.27.hdf5'
    >>> b.print_header5(infile)
    """
    try:
        f = h5py.File(infile, 'r')
    except:
        print('Check that the file provided is correct')
        return
    
    header = f[inhead]
    for hitem in list(header.attrs.items()): 
        print(hitem)
    f.close()
    return


def get_zminmaxs(zz,dz=None,verbose=False):
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
    verbose : bool
        If True write out warning messages

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
        if (verbose):
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

    
def get_simlabels(sims,labels=None):
    """
    Check that the arrays sims and inlabels have the same lenght,
    if not, get labels from the simulation names

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    labels : list of strings
        Array with the labels to be used

    Returns
    -----
    outlabels : list of strings
        Labels for simulations

    Examples
    ---------
    >>> import bahamasplot as bp
    >>> bp.get_simlabels(['HIRES/AGN_RECAL_nu0_L100N512_WMAP9','L400N1024/WMAP9/Sims/BAHAMAS'])
    """ 

    outlabels = labels
    # Check that the size of the arrays for the simulations and labels is the same
    if (labels == None or len(labels) != len(sims)):
        # Generate labels
        labels0 = [x.split('/')[0] for x in sims]        
        labels1 = [x.split('/')[-1] for x in sims]        

        outlabels = labels1 ; newlabel = ''
        for ii,label in enumerate(labels1):
            if (label == 'BAHAMAS'):
                newlabel = labels0[ii]

            val = '_msfof'
            if val in label:
                label = label.split(val)[0]

            val = '_beta_'
            if val in label:
                l1 = label.split(val)[1].replace('_','.')
                newlabel=',$\\beta=$'+l1
                label = label.split(val)[0]

            val = '_BH'
            if val in label:
                label = label.split(val)[0]
            
            val = '_n_'
            if val in label:
                l1 = label.split(val)[1].replace('_','.')
                newlabel=',$n=$'+l1+newlabel
                label = label.split(val)[0]

            val = '_dT_'
            if val in label:
                l1 = label.split(val)[1].replace('_','.')
                newlabel=',$\\Delta T=$'+l1+newlabel
                label = label.split(val)[0]

            val = '_mu_'
            if val in label:
                l1 = label.split(val)[1].replace('_','.')
                newlabel=',$\\mu=$'+l1+newlabel
                label = label.split(val)[0]

            val = 'ws_'
            if val in label:
                l1 = label.split(val)[1].replace('_','.')
                newlabel='$v_{\\rm wind}=$'+l1+newlabel
                label = label.split(val)[0]

            if (newlabel != ''):    
                outlabels[ii] = newlabel
        
    return outlabels

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


def get_outdirs(env,dirz=None,outdir=None,sim_label=None):
    """
    Get the output directories given the environment

    Parameters
    -----------
    env : string
        cosma, cosmalega, ari or arilega, to use the adecuate paths
    sim_label : string
        Simulation label to be used in path, otherwise *'_models'

    Returns
    -----
    outdir : string
       Path to output files
    dirz :string
       Path to list with z, snapshot
    dirplots : string
       Path to output plots

    Examples
    ---------
    >>> import bahamas as b
    >>> outdir, dirz, dirplots = b.get_outdirs('ari')
    """

    dirplots = None
    
    # Outdir
    if not outdir:
        if (env == 'ari'):
            outdir = dirbahamasari
        if (env == 'arilega'):
            outdir = dirbahamasarilega
        elif (env == 'cosma' or env == 'cosmalega'): 
            outdir = dirbahamascosma
        else:
            print('WARNING (b.get_outdirs): environment name not cosma(lega) or ari(lega)')
            return None,None,None

    new_dir = create_dir(outdir) 
    if not new_dir: return None,None,None

    # Dirz
    if not dirz:
        dirz = outdir 

    # Dirplots
    if not sim_label:
        if ('lega' in env):
            dirplots = outdir+'plots/published_models/'
        else:
            dirplots = outdir+'plots/compared_models/'
    else:
        dirplots = outdir+'plots/'+sim_label+'/'

    new_dir = create_dir(dirplots) 
    if not new_dir: return None,None,None

    return outdir,dirz,dirplots


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


def get_path2part(sim,env):
    """
    Get the path to the directory with all the particle data

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    env : string
        cosma, cosmalega, ari or arilega, to use the adecuate paths
 
    Returns
    -----
    path2part : string
       Path to data

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_path2part('L050N256/WMAP9/Sims/ex','cosma')
    >>> b.get_path2part('AGN_TUNED_nu0_L100N256_WMAP9','arilega')
    """

    # Simulation input
    if (env == 'ari'):
        path2part = dirbahamasari+sim+'/data/'
    elif (env == 'arilega'):
        path2part = dirbahamasarilega+sim+'/Data/Snapshots/'
    elif (env == 'cosma'):
        path2part = dirbahamascosma+sim+'/data/'
    elif (env == 'cosmalega'):
        path2part = dirbahamascosmalega+sim+'/data/'
    else:
        exit('get_path2data set to handle env=cosma, cosmalega, ari or arilega')

    return path2part    


def get_path2data(sim,env):
    """
    Get the path to the directory with the subhalo properties

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


def get_allparticle_files(snap,sim,env):
    """
    Get the particle files (in and out haloes)

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
    allfiles : boolean
       True if all files encountered given the numbers in the files names

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_allparticle_files(8,'L050N256/WMAP9/Sims/ex','cosma')
    >>> b.get_allparticle_files(27,'L400N1024/WMAP9/Sims/BAHAMAS','cosmalega')
    """

    allfiles = True

    # Simulation input
    path1 = get_path2part(sim,env)+'snapshot_'+str(snap).zfill(n0)

    # Get path to particle files
    paths = glob.glob(path1+'*/') 
    if (len(paths) == 1):
        path = paths[0]
    else:
        print('WARNING(b.get_allparticle_files): '+
              'more than one or none directories with root {}'.format(path1+'*/'))
        return None, False

    root = path+'snap_'+str(snap).zfill(n0)
    files = glob.glob(root+'*.hdf5')
    if (len(files)<1):
        print('WARNING (b.get_allparticle_files): no files in path {}'.format(path1+'*/'))
        return None, False

    numff =[int(ff.split('.')[1]) for ff in files]
    ind = np.argsort(numff)
    outff = [None]*len(files)
    for ii,iind in enumerate(ind):
        outff[ii]=files[iind]

        # Check if there are missing files
        if (ii != numff[iind]):
            print('WARNING (b.get_allparticle_files): missing file {}.{}.hdf5'.format(root,ii))
            allfiles = False

    return outff, allfiles


def get_particle_files(snap,sim,env):
    """
    Get the halo particle files

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
    allfiles : boolean
       True if all files encountered given the numbers in the files names

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_particle_files(8,'L050N256/WMAP9/Sims/ex','cosma')
    >>> b.get_particle_files(27,'L400N1024/WMAP9/Sims/BAHAMAS','cosmalega')
    """

    allfiles = True

    # Simulation input
    path1 = get_path2data(sim,env)+'particledata_'+str(snap).zfill(n0)

    # Get path to particle files
    paths = glob.glob(path1+'*/') 
    if (len(paths) == 1):
        path = paths[0]
    else:
        print('WARNING(b.get_particle_files): more than one or none directories with root {}'.format(path1+'*/'))
        return None, False

    root = path+'eagle_subfind_particles_'+str(snap).zfill(n0) 
    files = glob.glob(root+'*.hdf5')
    if (len(files)<1):
        print('WARNING (b.get_particle_files): no files in path {}'.format(path1+'*/'))
        return None, False

    numff =[int(ff.split('.')[1]) for ff in files]
    ind = np.argsort(numff)
    outff = [None]*len(files)
    for ii,iind in enumerate(ind):
        outff[ii]=files[iind]

        # Check if there are missing files
        if (ii != numff[iind]):
            print('WARNING (b.get_particle_files): missing file {}.{}.hdf5'.format(root,ii))
            allfiles = False

    return outff, allfiles


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
    allfiles : boolean
       True if all files encountered given the numbers in the files names

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_subfind_files(8,'L050N256/WMAP9/Sims/ex','cosma')
    >>> files, allfiles = b.get_subfind_files(27,'L400N1024/WMAP9/Sims/BAHAMAS','cosmalega')
    """

    allfiles = True

    # Simulation input
    path1 = get_path2data(sim,env)+'groups_'+str(snap).zfill(n0)

    # Get path to subfind files
    paths = glob.glob(path1+'*/')
    if (len(paths) == 1):
        path = paths[0]
    else:
        print('WARNING (b.get_subfind_files): more than one or none directories with root {}'.format(path1+'*/'))
        return None, False

    root = path+'eagle_subfind_tab_'+str(snap).zfill(n0)
    files = glob.glob(root+'*.hdf5')
    if (len(files)<1):
        print('WARNING (b.get_subfind_files): no files in path {}'.format(path1+'*/'))
        return None, False

    numff =[int(ff.split('.')[1]) for ff in files]
    ind = np.argsort(numff)
    outff = [None]*len(files)
    for ii,iind in enumerate(ind):
        outff[ii]=files[iind]

        # Check if there are missing files
        if (ii != numff[iind]):
            print('WARNING (b.get_subfind_files): missing file {}.{}.hdf5'.format(root,ii))
            allfiles = False

    return outff, allfiles


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
    omega0, omegab, lambda0, h0, volume : floats
        Cosmological parameters and volume

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

    volume = header.attrs['BoxSize']**3 #(Mpc/h)^3
    
    f.close()
    return omega0, omegab, lambda0, h0, volume


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
        print('WARNING (b.table_z_sn): {} not containing expected files'.format(path))
        return None
        
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
        if tablez is None: return -999.

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


def get_snap(zz,sim,env,zmin=None,zmax=None,dirz=None):
    """
    Get the closest snapshot given a redshift and
    a simulation name, within some limits

    Parameters
    -----------
    zz : float
        Redshift
    sim : string
        Name of the Bahamas simulation
    zmin : float
        Minimum redsfhit to look for the sanpshot
    zmax : float
        Maximum redsfhit to look for the sanpshot
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

    # Check if the redshift range has been provided
    if (zmin == None and zmax == None):
        zmin,zmax = get_zminmaxs([zz])
    elif (zmin == None):
        zmin,dum = get_zminmaxs([zz])
    elif (zmax == None):
        dum,zmax = get_zminmaxs([zz])

    # Table with increasing redshifts and corresponding snapshots
    if (dirz == None):
        tablez = path+tblz
    else:
        tablez = dirz+sim+'/'+tblz

    if (not os.path.isfile(tablez)):
        # Generate the table if it doesn't exist
        tablez = table_z_sn(sim,env,dirz=dirz)
        if tablez is None: return -999.,-999.

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


def get_cenids(snap,sim,env,Testing=False,nfiles=2):
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
    Testing: boolean
        True or False
    nfiles : integer
        Number of files to be considered for testing

    Returns
    -----
    cenids : numpy array int
        Global indexes for central galaxies

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_cenids(31,'HIRES/AGN_TUNED_nu0_L050N256_WMAP9','ari')
    >>> b.get_cenids(27,'L400N1024/WMAP9/Sims/BAHAMAS','cosmalega',Testing=True)
    """

    # Simulation input
    files, allfiles = get_subfind_files(snap,sim,env)
    if allfiles is False: return -999.

    # Cycle through the files
    for ii,ff in enumerate(files):
        if (Testing and ii>=nfiles): break
        stop_if_no_file(ff)

        f = h5py.File(ff, 'r')
        haloes = f['FOF/FirstSubhaloID'][:]
        if (ii == 0):
            cenids = np.unique(haloes)
        else:
            cenids = np.append(cenids, np.unique(haloes))  

    if (not is_sorted(cenids)):
        print('WARNING (b.get_cenids): Not ordered indeces {}'.format(path))
        return -999.

    return cenids



def get_prop(snap,sim,env,propdef,Testing=False,nfiles=2):
    """
    Get an array with a given property

    Parameters
    -----------
    snap : integer
        Snapshot 
    sim : string
        Simulation name
    env : string
        ari or cosma, to use the adecuate paths
    propdef : string
        Name of the property, including path within hdf5 file
    Testing: boolean
        True or False
    nfiles : integer
        Number of files to be considered for testing

    Returns
    -----
    fofhmass : numpy array float
        Property for FOF groups, 10^10Msun/h 

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_prop(27,'L400N1024/WMAP9/Sims/BAHAMAS','cosmalega','FOF/Group_M_Crit200',Testing=True)
    """

    # Simulation input
    files, allfiles = get_subfind_files(snap,sim,env)
    if allfiles is False: return -999.

    # Cycle through the files
    for ii,ff in enumerate(files):
        if (Testing and ii>=nfiles): break
        stop_if_no_file(ff)

        f = h5py.File(ff, 'r')

        if (ii == 0):
            try:
                prop = f[propdef][:]
            except:
                print('WARNING (bahamas): no {} found in {}'.format(propdef,ff))
                return None
        else:
            prop = np.append(prop, f[propdef][:], axis=0)

    return prop


def resolution(sim,env,zz=0.,msunh=True,dirz=None,verbose=True):
    """
    Get the mass resolution of a simulation

    Parameters
    -----------
    sim : string
        Simulation name
    env : string
        ari, arilega or cosma, to use the adecuate paths
    zz : float
        Redshift to look (should be the same for all)
    msunh : boolean
        If True, mass in Msun/h units
    dirz : string
        Alternative path to table with z and snapshot.
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

    snap, z_snap = get_snap(zz,sim,env,dirz=dirz)

    # Simulation input
    files, allfiles = get_particle_files(snap,sim,env)
    if (files is None): return -999., -999.

    f= h5py.File(files[0],'r') #;print(files[0])
    header=f['Header']
    masstable = header.attrs['MassTable']

    itype = ptypes.index('DM') 
    mdm = masstable[itype]
    if (mdm<0):
        print('WARNING: negative or 0 input mass, returning -999.')
        return -999., -999.
    else:
        omega0, omegab, lambda0, h0, volume = get_cosmology(sim,env)

        if (msunh):
            mdm = mdm*10**10
        else:
            mdm = mb2msun(mdm,h0)
            
        mgas = mdm*omegab/(omega0-omegab)

        if (verbose):
            print('Particle resolution of sim.: {}'.format(sim))
            if msunh:
                print('mdm (Msun/h) = {:.2e}, mgas (Msun/h)=mb= {:.2e}'.format(mdm,mgas))
            else:
                print('mdm (Msun) = {:.2e}, mgas (Msun)=mb= {:.2e}'.format(mdm,mgas))
        
        return mdm,mgas


def get_min_sfr(sim,env,zz=0.,A=1.515*1e-4,gamma=5/3,fg=1.,n=1.4,dirz=None,verbose=True):
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
    dirz : string
        Alternative path to table with z and snapshot.
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
    mdm, mgas = resolution(sim,env,dirz=dirz,verbose=False)

    apc = A/1e6 # Msun yr^-1 pc^-2

    nHm = nHsf * 100.**3 # m^-3
    P = nHm*const.k_B.value*Tsf
    
    gp_si = np.sqrt(gamma*fg*P/const.G.value)
    gp = gp_si*const.pc.value**2/const.M_sun.value

    minsfr = mgas*apc*np.power(gp,n-1)
    if verbose:
        print('  Min. theoretical log10(SFR (Msun/Gyr)) = {:2f}'.format(np.log10(minsfr)+9))
    
    return minsfr


def get_nh(zz,massdef,sim,env,mmin=9.,mmax=16.,dm=0.1,
           dirz=None,outdir=None,Testing=True,verbose=False):
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
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    Testing : boolean
        Calculations on part or all the simulation
    verbose : boolean
        True to print the resolution

    Returns
    -----
    outfile : string
        Path to file with: 

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.get_nh(0.5,'FOF/Group_M_Mean200',sim,'arilega',outdir='/hpcdata4/arivgonz/BAHAMAS/')
    '''

    # Convert / to _ in massdef name
    mnom = massdef.replace('/','_')
    
    # Get snapshot
    snap, z_snap = get_snap(zz,sim,env,dirz=dirz)

    # Output file
    if outdir:
        path = outdir+sim
    else:
        path = get_dirb(env)+sim
    if (not os.path.exists(path)):
        os.makedirs(path)

    if Testing:
        outfil = path+'/nh_'+mnom+'_sn'+str(snap)+'_dm'+str(dm)+'_test.txt'
    else:
        outfil = path+'/nh_'+mnom+'_sn'+str(snap)+'_dm'+str(dm)+'.txt'
        
    if (os.path.isfile(outfil)):
        #data = subprocess.run(['wc', '-l', outfil], capture_output=True) #python3.7
        process = subprocess.run(['wc', '-l', outfil], universal_newlines=True, 
                                 stdout=subprocess.PIPE).stdout
        nlines = int(process.split(' ')[0])
        if (nlines > 0):
            if verbose: print('b.get_nh, file already exists: ',outfil)
            return outfil

    # Get the halo mass
    mh = get_prop(snap,sim,env,massdef,Testing=Testing)
    ind = np.where(mh > 0.)
    lmh = np.log10(mh[ind]) + 10. # log10(M/Msun/h)

    # Bins in halo mass
    edges = np.array(np.arange(mmin,mmax,dm))
    mhist = edges[1:]-0.5*dm
    nh  = np.zeros(shape=(len(mhist)))

    elow  = edges[:-1]
    ehigh = edges[1:]

    # Number of haloes
    H, bins_edges = np.histogram(lmh,bins=edges)
    nh[:] = nh[:] + H

    # Output only bins with haloes
    indh = np.where(nh > 0) 
    tofile = np.array([mhist[indh], elow[indh], ehigh[indh], nh[indh]])

    # Get cosmology and volume
    omega0, omegab, lambda0, h0, volume = get_cosmology(sim,env)

    # Write to output file
    with open(outfil, 'w') as outf:
        # Header
        outf.write('# '+sim+', snapshot='+str(snap)+' \n')
        outf.write('# Volume ((Mpc/h)^3)='+str(volume)+', dm='+str(dm)+' \n' )
        outf.write('# omega0='+str(omega0)+', omegab='+str(omegab)+
                   ', lambda0='+str(lambda0)+', h0='+str(h0)+' \n')
        outf.write('# log(Mh/Msun/h)_midpoint, log(Mh)_low, log(Mh)_high, Number of haloes \n')

        # Data
        np.savetxt(outf,tofile.T,fmt='%.5f %.5f %.5f %.0f')

    return outfil


def get_propfunc(zz,propdefs,proplabel,sim,env,ptype=['star'],mmin=9.,mmax=16.,dm=0.1,
           dirz=None,outdir=None,Testing=True):
    '''
    Calculate the functions of variables (normalized histograms in log scales)

    Parameters
    -----------
    zz : float
        Redshift to get the number of haloes
    propdef : array of string
        Name of the properties (needs to include the FOF/ or Subhalo/ part)
    proplabel : string
        Short name to identify the properties. E.g. 'mass'
    sim : string
        Name of the simulation
    env : string
        ari, arilega or cosma, to use the adecuate paths
    ptype : array of string
        array containing one of the allowed ptypes
    mmin : float
        Mimimum mass to be considered
    mmax : float
        Maximum mass to be considered
    dm : float
        Intervale step for the halo mass
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    Testing : boolean
        Calculations on part or all the simulation

    Returns
    -----
    outfile : string
        Path to file with the funtion of 'propdefs'

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.get_propfunc(31,['FOF/Group_M_Mean200'],'mass',sim,'arilega',outdir='/hpcdata4/arivgonz/Junk/')
    '''

    itype = ptypes.index(ptype)
    
    keys = [key for key, value in unitdefault.items()]
    if (proplabel not in keys):
        print('WARNING (get_propfunc): '+proplabel+' not in ('+','.join(keys)+')')
        return None

    # Get snapshot
    zmin,zmax = get_zminmaxs([zz])
    snap, z_snap = get_snap(zz,zmin,zmax,sim,env,dirz=dirz)
    
    # Output file
    outdir, dirz, plotdir = get_outdirs(env,dirz=dirz,outdir=outdir)
    path = outdir+sim ; create_dir(path)
    outfil = path+'/'+proplabel+'F_z'+str(z_snap).replace('.','_')+ \
             '_min'+str(mmin).replace('.','_')+'_max'+str(mmax).replace('.','_')+ \
             '_dm'+str(dm).replace('.','_')+'.hdf5'

    # Bins in halo mass
    edges = np.array(np.arange(mmin,mmax,dm))
    mhist = edges[1:]-0.5*dm
    elow  = edges[:-1]
    ehigh = edges[1:]

    # Get cosmology and volume
    omega0, omegab, lambda0, h0, volume = get_cosmology(sim,env)

    # Check which properties are already in the file
    if (os.path.isfile(outfil)):
        print('File in place:',outfil,propdefs)
        #check if the property is already there
        # if so, continue  (#here do I want this?)

    # Separate the FOF and Subhalo properties
    pfof = [p for p in propdefs if 'FOF' in p]
    psub = [p for p in propdefs if 'Subhalo' in p]

    # Initialize matrices for arrays
    propf = np.zeros((len(propdefs),len(mhist)))
    iprop = -1
    
    if (len(pfof) > 0):
        for propdef in pfof:
            iprop += 1
            pdef = propdef.split('/')[1]
            prop = get_fofprop(snap,sim,env,pdef,Testing=Testing)
            if (prop is None): continue
            
            ind = np.where(prop > 0.)
            if (np.shape(ind)[1]<1): continue
            
            if (proplabel == 'mass'):
                lprop = np.log10(prop[ind]) + 10. # log10(M/Msun/h)
            else:
                lprop = np.log10(prop[ind])
            
            # Numbers of haloes
            H, bins_edges = np.histogram(lprop,bins=edges)
            propf[iprop][:] = propf[iprop][:] + H
        
    if (len(psub) > 0):
        #here work through subvolumes
        print(propdef)

    # Take logs
    ind = np.where(propf > 0)
    propf[ind] = np.log10(propf[ind]/volume/dm)

    # Output
    hf = h5py.File(outfil, 'w') #here if appending columns here it'll be the place

    # Output header
    head = hf.create_dataset('header',(100,))
    head.attrs[u'sim']          = sim
    head.attrs[u'snapshot']     = snap
    head.attrs[u'redshift']     = z_snap
    head.attrs[u'h0']           = h0

    # Output data
    units = unitdefault[proplabel]
    hfdat = hf.create_group('data')

    pn = 'midpoint'
    hfdat.create_dataset(pn,data=mhist)
    hfdat[pn].dims[0].label = 'log10('+proplabel+'_'+pn+' / '+units+')'

    pn = 'low'
    hfdat.create_dataset(pn,data=elow)
    hfdat[pn].dims[0].label = 'log10('+proplabel+'_'+pn+' / '+units+')'

    pn = 'high'
    hfdat.create_dataset(pn,data=ehigh)
    hfdat[pn].dims[0].label = 'log10('+proplabel+'_'+pn+' / '+units+')'

    for ii, propdef in enumerate(propdefs):
        propii = propf[ii][:]
        hfdat.create_dataset(propdef,data=ehigh)
        hfdat[propdef].dims[0].label = 'log10(N/dlog'+propdef.replace('/','_')+'/'+unitdefault['volume']+')'

    hf.close()    
    return outfil


if __name__== "__main__":
    dirz = None ; outdir = None
    snap = 31
    zz = 3.

    env = 'arilega'
    #env = 'cosmalega'

    if (env == 'cosmalega'):
        sim = 'L400N1024/WMAP9/Sims/BAHAMAS'
        dirz = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'
        outdir = '/cosma6/data/dp004/dc-gonz3/Junk/'
    if (env == 'arilega'):
        #sim = 'HIRES/AGN_RECAL_nu0_L100N512_WMAP9'
        sim = 'AGN_TUNED_nu0_L400N1024_WMAP9'
        dirz = '/hpcdata4/arivgonz/BAHAMAS/'
        outdir = '/hpcdata4/arivgonz/Junk/'
    if (env == 'ari'):
        sim = 'L050N256/WMAP9/Sims/ws_324_23_mu_7_05_dT_8_35_n_75_BH_beta_1_68_msfof_1_93e11'

    #print(get_zminmaxs([0.,1.],dz=0.5))
    #print(get_simlabels(['AGN_TUNED_nu0_L100N256_WMAP9',
    #               'HIRES/AGN_RECAL_nu0_L100N512_WMAP9',
    #               'L400N1024/WMAP9/Sims/BAHAMAS']))
    #print(get_outdirs(env,dirz=dirz,outdir=outdir))
    #print(table_z_sn(sim,env,dirz=dirz))
    #print(get_z(27,sim,env,dirz=dirz))
    #print(get_z(-1,sim,env,dirz=dirz))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz)
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz,zmax=[3.8])
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz,zmin=[2.8])
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz,zmin=[2.8],zmax=[3.8])
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #print(get_allparticle_files(snap,sim,env))
    #print(get_cenids(snap,sim,env))
    #print(get_prop(snap,sim,env,'FOF/Group_M_Crit200'))
    print(resolution(sim,env,dirz=dirz))
    #print('log10(SFR (Msun/Gyr)) = {:2f}'.format(np.log10(get_min_sfr(sim,env,dirz=dirz))+9))
    #print(get_nh(zz,'FOF/Group_M_Mean200',sim,env,dirz=dirz,outdir=outdir))
    #print(get_propfunc(zz,['FOF/Group_M_Mean200','FOF/m2'],
    #                   'mass',sim,env,ptype='DM',dirz=dirz,outdir=outdir))

    #infile = '/hpcdata0/simulations/BAHAMAS/AGN_TUNED_nu0_L100N256_WMAP9/Data/Snapshots/snapshot_026/snap_026.27.hdf5'
    #infile = '/hpcdata0/simulations/BAHAMAS/AGN_TUNED_nu0_L100N256_WMAP9/Data/EagleSubGroups_5r200/groups_026/eagle_subfind_tab_026.0.hdf5'
    #print(print_h5attributes(infile,'Constants'))
