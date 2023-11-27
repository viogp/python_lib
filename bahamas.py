import os.path
import numpy as np
import h5py
import glob
import subprocess
import pandas as pd
import astro as ast
from cosmosim import get_r, get_vr, get_vphi, get_vlos
from astropy import constants as const
import iotools as io

ptypes = ['gas','DM','bp1','bp2','star','BH']

nogroup = 1073741824.

# Default units in Bahamas
unitdefault = { 
    'volume': '(Mpc/h)^3',
    'mass': '10^10 Msun/h',
    'sfr': 'Msun/h/yr',     ##here not sure
    'ssfr': 'to be filled'
}

BH_seed_mass = 7.63*10**5 #Msun/h

#dirbahamasarilega = '/hpcdata0/simulations/BAHAMAS/'
#dirbahamasari = '/hpcdata0/arivgonz/BAHAMAS/'
#dirobsari = '/hpcdata0/Obs_Data/'

dirbahamasarilega = '/beegfs1/backup/beegfs2/hpcdata0_backup/simulations/BAHAMAS/'
# '/mnt/autofs/beegfs2/hpcdata0_backup/simulations/'
# '/beegfs2/hpcdata0_backup/simulations/BAHAMAS/' # Output: '/beegfs2/arivgonz/BAHAMAS/'
dirbahamasari = '/enc1/hpcddn/hpcdata3/arivgonz/BAHAMAS/' #from havok
dirobsari = '/users/arivgonz/trunk/Obs_Data/'
dirplotari = '/users/arivgonz/output/'

dirbahamascosmalega = '/cosma7/data/dp004/Eagle/jsTestRuns/BAHAMAS/'
#'/cosma7/data/dp004/fkgm22/ANTILLES/L100N256/WMAP9/' 
dirbahamascosma = '/cosma7/data/dp004/dc-gonz3/BAHAMAS/'
dirobscosma = '/cosma7/data/dp004/dc-gonz3/BAHAMAS/Obs_Data/'

dirbahamaslap = '/home/violeta/soil/BAHAMAS/'
dirobslap = '/home/violeta/soil/Obs_Data/'
dirplotlap = '/home/violeta/buds/'

tblz = 'snap_z.txt'

defaultdz = 0.25

n0 = 3

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
            print('\n WARNING (get_zminmaxs): Setting to default dz={}'.format(defaultdz))
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
    >>> import bahamas as b
    >>> b.get_simlabels(['HIRES/AGN_RECAL_nu0_L100N512_WMAP9','L400N1024/WMAP9/Sims/BAHAMAS'])
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

def mb2msun(massb,h0,verbose=False):
    """
    Provide bahamas masses as M/Msun

    Parameters
    -----------
    massb : float
        Mass in 10^10Msun/h units
    h0 : float
        Hubble constant
    verbose : bool
        True for informative messages

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

    mass = -999.
    if(massb <= 0):
        if (verbose): print('\n WARNING (get_mb2msun): Input mass <= 0, returning -999.')
    else:
        mass = massb*np.power(10.,10.)/h0
        
    return mass


def mb2lmsun(massb,h0,verbose=False):
    """
    Provide bahamas masses as log10(M/Msun)

    Parameters
    -----------
    massb : float
        Mass in 10^10Msun/h units
    h0 : float
        Hubble constant
    verbose : bool
        True for informative messages

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

    lmass = -999.
    if(massb <= 0):
        if (verbose): print('\n WARNING (get_mb2lmsun): Input mass <= 0, returning -999.')
    else:
        lmass = np.log10(massb) + 10. - np.log10(h0)
        
    return lmass


def get_dirb(env):
    """
    Get the Bahamas directory given the environment

    Parameters
    -----------
    env : string
        cosma, cosmalega, ari, arilega or lap, to use the adecuate paths
 
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
    elif (env == 'lap'):
        dirb = dirbahamaslap
    else:
        print('\n WARNING (b.get_dirb): environment name not cosma, cosmalega, ari or arilega')

    return dirb


def get_outdirs(env,dirz=None,outdir=None,sim_label=None):
    """
    Get the output directories given the environment

    Parameters
    -----------
    env : string
        cosma, cosmalega, ari, arilega or lap, to use the adecuate paths
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
        elif (env == 'arilega'):
            outdir = dirbahamasarilega
        elif (env == 'cosma' or env == 'cosmalega'): 
            outdir = dirbahamascosma
        elif (env == 'lap'):
            outdir = dirbahamaslap
        else:
            print('\n WARNING (b.get_outdirs): environment name not cosma(lega) or ari(lega)')
            return None,None,None

    new_dir = io.create_dir(outdir) 
    if not new_dir: return None,None,None

    # Dirz
    if not dirz:
        dirz = outdir 

    # Dirplots
    outdir1 = outdir
    if ('ari' in env):
        outdir1 = dirplotari

    if not sim_label:
        if ('lega' in env):
            dirplots = outdir1+'plots/published_models/'
        else:
            dirplots = outdir1+'plots/compared_models/'
    else:
        dirplots = outdir1+'plots/'+sim_label+'/'

    new_dir = io.create_dir(dirplots) 
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
    elif (env == 'lap'):
        dirobs = dirobslap

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
        if (sim==None):
            path2part = dirbahamascosmalega+'data/'
        else:
            path2part = dirbahamascosmalega+sim+'/data/'
    else:
        exit('b.get_path2part set to handle env=cosma, cosmalega, ari or arilega')

    return path2part    


def get_path2data(sim,env):
    """
    Get the path to the directory with the subhalo properties

    Parameters
    -----------
    sims : string
        Name of the simulation
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
        if (sim==None):
            path2data = dirbahamascosmalega+'data/'
        else:
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
        print('\n WARNING(b.get_allparticle_files): '+
              'more than one or none directories with root {}'.format(path1+'*/'))
        return None, False

    root = path+'snap_'+str(snap).zfill(n0)
    files = glob.glob(root+'*.hdf5')
    if (len(files)<1):
        print('\n WARNING (b.get_allparticle_files): no files in path {}'.format(path1+'*/'))
        return None, False

    numff =[int(ff.split('.')[1]) for ff in files]
    ind = np.argsort(numff)
    outff = [None]*len(files)
    for ii,iind in enumerate(ind):
        outff[ii]=files[iind]

        # Check if there are missing files
        if (ii != numff[iind]):
            print('\n WARNING (b.get_allparticle_files): missing file {}.{}.hdf5'.format(root,ii))
            allfiles = False

    return outff, allfiles


def get_particle_files(snap,sim,env,subfind=True):
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
    subfind : boolean
        True for Subfind data, False for snapshot

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
    if subfind:
        path1 = get_path2data(sim,env)+'particledata_'+str(snap).zfill(n0)
    else:
        path1 = get_path2part(sim,env)+'snapshot_'+str(snap).zfill(n0)

    # Get path to particle files
    paths = glob.glob(path1+'*/') 
    if (len(paths) == 1):
        path = paths[0]
    else:
        print('\n WARNING(b.get_particle_files): more than one or none directories with root {}'.format(path1+'*/'))
        return None, False

    if subfind:
        root = path+'eagle_subfind_particles_'+str(snap).zfill(n0)
    else:
        root = path+'snap_'+str(snap).zfill(n0)

    files = glob.glob(root+'*.hdf5')
    if (len(files)<1):
        print('\n WARNING (b.get_particle_files): no files in path {}'.format(path1+'*/'))
        return None, False

    numff =[int(ff.split('.')[1]) for ff in files]
    ind = np.argsort(numff)
    outff = [None]*len(files)
    for ii,iind in enumerate(ind):
        outff[ii]=files[iind]

        # Check if there are missing files
        if (ii != numff[iind]):
            print('\n WARNING (b.get_particle_files): missing file {}.{}.hdf5'.format(root,ii))
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
        print('\n WARNING (b.get_subfind_files): more than one or none directories with root {}'.format(path1+'*/'))
        return None, False

    root = path+'eagle_subfind_tab_'+str(snap).zfill(n0)
    files = glob.glob(root+'*.hdf5')
    if (len(files)<1):
        print('\n WARNING (b.get_subfind_files): no files in path {}'.format(path1+'*/'))
        return None, False
    
    numff =[int(ff.split('.')[1]) for ff in files]
    ind = np.argsort(numff)
    outff = [None]*len(files)
    for ii,iind in enumerate(ind):
        outff[ii]=files[iind]

        # Check if there are missing files
        if (ii != numff[iind]):
            print('\n WARNING (b.get_subfind_files): missing file {}.{}.hdf5'.format(root,ii))
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
    omega0, omegab, lambda0, h0, boxside : floats
        Cosmological parameters and boxside

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_cosmology('AGN_TUNED_nu0_L100N256_WMAP9','arilega')
    >>> b.get_cosmology('L050N256/WMAP9/Sims/ex','cosma')
    """

    # Simulation input
    path = get_path2data(sim,env)

    # Initialize arrays for z and sn
    files = glob.glob(path+'groups_*/group_tab*')
    infile = files[0] #; print(infile) 
    f = h5py.File(infile, 'r')
    header = f['Header']
    #print(list(header.attrs.items()))
    
    omega0 = header.attrs['Omega0']
    omegab = header.attrs['OmegaBaryon']
    lambda0 = header.attrs['OmegaLambda']
    h0 = header.attrs['HubbleParam']

    boxside = header.attrs['BoxSize'] # Mpc/h
    
    f.close()
    return omega0, omegab, lambda0, h0, boxside


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
        if (sim==None):
            dirz = dirz
        else:
            dirz = dirz+sim+'/'
        if (not os.path.exists(dirz)):
            os.makedirs(dirz)
        tablez = dirz+tblz

    # Initialize arrays for z and sn
    dirs = glob.glob(path+'groups_0*')
    if (len(dirs) < 1) :
        print('\n WARNING (b.table_z_sn): {} not containing expected files'.format(path))
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
            print('\n WARNING: Files missing in {}'.format(dir))
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
        if (sim==None):
            tablez = dirz+tblz
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
            print('\n WARNING: no snapshot number {} found in {}'.format(snap,tablez))
            return -999.            
    else:
        sns = np.asarray(snsf,dtype=int)
        # Find the index of the snapshot number
        ind = np.where(sns == snap)
        if (np.shape(ind)[1]<1):
            print('\n WARNING: no snapshot number {} found in {}'.format(snap,tablez))
            return -999.
        elif (np.shape(ind)[1]>1):
            print('\n WARNING: several snapshot number {} found in {}'.format(snap,tablez))
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
    >>> b.get_snap(2.8,'AGN_TUNED_nu0_L100N256_WMAP9','ari')
    >>> b.get_snap(12.6,'L050N256/WMAP9/Sims/ex','cosma')
    >>> (2, 12.5)
    >>> snap, z_snap = b.get_snap(20.,'L050N256/WMAP9/Sims/ex','cosma')
    """

    # Check if the redshift range has been provided
    if (zmin == None and zmax == None):
        zmin,zmax = get_zminmaxs([zz])
    elif (zmin == None):
        zmin,dum = get_zminmaxs([zz])
    elif (zmax == None):
        dum,zmax = get_zminmaxs([zz])

    # Table with increasing redshifts and corresponding snapshots
    if (dirz == None):
        # Simulation input
        path = get_path2data(sim,env)
        tablez = path+tblz
    else:
        if (sim==None):
            tablez = dirz+tblz
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
            print('\n WARNING: {} far from single {} found in {}'.format(zz,zzs,tablez))
            return -999,-999.
    sns = np.asarray(lsns,dtype=int)

    # Find the closest redshift
    if ((zz < zzs[0] and zmax < zzs[0]) or (zz<0.)):
        print('\n WARNING: {}<{}, min. z found in {}'.format(zz,zzs[0],tablez))
        return -999,-999.
    elif(zz > zzs[-1] and zmin > zzs[-1]):
        print('\n WARNING: {}>{}, max. z found in {}'.format(zz,zzs[-1],tablez))
        return -999.,-999.
    else:
        idx = (np.abs(zzs - zz)).argmin()
        if (zmin < zzs[idx] and zzs[idx] < zmax):
            return sns[idx],zzs[idx]
        else:
            print('\n WARNING: z={} outside range {}<z<{}, {}'.format(zzs[idx],zmin,zmax,tablez))
            return -999,-999.


def get_cenids(snap,sim,env,Testing=False,nfiles=2):
    """
    Get the list of indexes for central galaxies
    using FOF/FirstSubhaloID (index of first subhalos, or centrals)

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
    >>> b.get_cenids(31,'HIRES/AGN_TUNED_nu0_L050N256_WMAP9','arilega')
    >>> b.get_cenids(27,'L400N1024/WMAP9/Sims/BAHAMAS','cosmalega',Testing=True)
    """

    # Simulation input
    files, allfiles = get_subfind_files(snap,sim,env)
    if allfiles is False: return -999.

    # Cycle through the files
    ii = 0; cenids = [];
    for ff in files:
        io.stop_if_no_file(ff)

        f = h5py.File(ff, 'r')

        # Index of first subhalo in SubHalo list (starts at 0) 
        first = f['FOF/FirstSubhaloID'][:]

        # Haloes with no substructure are not well defined
        nsubh = f['FOF/NumOfSubhalos'][:]
        ind = np.where(nsubh>0)
        if (np.shape(ind)[1]>0):
            cenh = first[ind]            
            if (abs(len(cenh) - len(np.unique(first))) > 1):
                print('\n WARNING (b.get_cenids): Not unique central IDs {}'.format(ff))
                return None
        
            if (ii == 0):
                cenids = cenh
            else:
                cenids = np.append(cenids, cenh)  

            if (not io.is_sorted(cenids)):
                print('\n WARNING (b.get_cenids): Not ordered indeces {}'.format(ff))
                return None

            ii += 1
            if (Testing and ii>=nfiles): break
            
        f.close()    

    return cenids



def get_subfind_prop(snap,sim,env,propdef,proptype=None,Testing=False,nfiles=2,verbose=False):
    """
    Get an array with a given property from the Subfind output

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
    proptype : string
        'DM', 'star', 'gas', 'BH', etc. for relevant properties
    Testing: boolean
        True or False
    nfiles : integer
        Number of files to be considered for testing
    verbose : boolean
        True to write first Subfind file out

    Returns
    -----
    prop : numpy array float
        Property within Subfind files

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_subfind_prop(27,'L400N1024/WMAP9/Sims/BAHAMAS','cosmalega',
                           'FOF/Group_M_Crit200',Testing=True)
    """

    # Simulation input
    files, allfiles = get_subfind_files(snap,sim,env)
    if allfiles is False: return -999.
    if verbose: print('get_subfind_prop: First Subfind file is {}'.format(files[0]))

    if (proptype is not None):
        itype = ptypes.index(proptype)
        
    # Cycle through the files
    for ii,ff in enumerate(files):
        if (Testing and ii>=nfiles): break
        io.stop_if_no_file(ff)

        f = h5py.File(ff, 'r')
        if (ii == 0):
            if (proptype is None):
                try:
                    prop = f[propdef][:]
                except:
                    print('\n WARNING (bahamas): no {} found in {}'.format(propdef,ff))
                    return None
            else:
                try:
                    prop = f[propdef][:,itype]
                except:
                    print('\n WARNING (bahamas): no {} found in {}'.format(propdef,ff))
                    return None
        else:
            if (proptype is None):
                prop = np.append(prop, f[propdef][:], axis=0)
            else:
                prop = np.append(prop, f[propdef][:,itype], axis=0)
    return prop


def resolution(sim,env,zz=0.,msunh=True,dirz=None,verbose=False):
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
        print('\n WARNING: negative or 0 input mass, returning -999.')
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
        if (sim==None):
            path = outdir
        else:
            path = outdir+sim
    else:
        if (sim==None):
            path = get_dirb(env)
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
        print('\n WARNING (get_propfunc): '+proplabel+' not in ('+','.join(keys)+')')
        return None

    # Get snapshot
    zmin,zmax = get_zminmaxs([zz])
    snap, z_snap = get_snap(zz,zmin,zmax,sim,env,dirz=dirz)
    
    # Output file
    outdir, dirz, plotdir = get_outdirs(env,dirz=dirz,outdir=outdir)
    if (sim==None):
        path = outdir
    else:
        path = outdir+sim
    io.create_dir(path)
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


def get_m500_file(outdir,sim,snap):
    '''
    Get the name and existance check of the map_m500 file

    Parameters
    ----------
    outdir: string
       Directory to write or find the file
    sim: string
       Simulation name or path
    snap: string
       Snapshot of the simulation
    
    Returns
    -------
    outfile: string
       Name of the map_m500 file
    file_exists: boolean
       True if file exists
    '''
    
    if (sim==None):
        outdir2 = outdir
    else:
        outdir2 = outdir+sim+'/'
    dir_exists = io.create_dir(outdir2)
    outfile = outdir2+'m500_snap'+str(snap)+'.hdf5'
    file_exists = io.check_file(outfile)

    return outfile, file_exists


def map_m500(snap,sim,env,ptype='gas',overwrite=False,mlim=0.,dirz=None,outdir=None,
             Testing=True, verbose=False):
    '''
    Map the mass of a given type of (subfind) particles into FOF's R500 

    Parameters
    -----------
    snap : int
        Snapshot number
    sim : string
        Name of the simulation
    env : string
        ari, arilega or cosma, to use the adecuate paths
    ptype : string
        Name of one of the allowed ptypes, 0:gas, 4: stars, 5:BH
    overwrite : boolean
        If True the output file will be overwritten
    mlim : float
        mass limit for M500 [Msun/h]
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    Testing : boolean
        Calculations on part or all the simulation
    verbose : boolean
        To output extra information

    Returns
    -----
    prop : float array
        Mapped property

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.get_propfunc(31,sim,'arilega')
    '''

    if (ptype == 'DM'):
        print('\n WARNING (bahamas.map_m500): For DM, use directly the halo masses')
        return None

    # Type of particles to be read
    itype = ptypes.index(ptype) # 0:gas, 1:DM, 4: stars, 5:BH
    inptype = 'PartType'+str(itype)
    if verbose: print('\n {}: {} \n'.format(ptype,inptype))
    
    # Output file
    outfile, file_exists = get_m500_file(outdir,sim,snap)
    if(overwrite): file_exists = False 

    # Check if the dataset already exists
    nompartmass = 'm500_'+ptype
    if (file_exists):
        f = h5py.File(outfile, 'r')
        e = 'data/'+nompartmass in f
        f.close()
        if (e):
            print('\n WARNING (b.map_m500): {} already in file.'.format(nompartmass))
            return outfile

    # Get particle files
    files, allfiles = get_particle_files(snap,sim,env)

    if (not allfiles):
        print('\n WARNING (b.map_m500): no adequate particle files found, {}, {}'.
              format(snap,env))
        return None
    if (Testing): files = [files[0]]
    if verbose: print('\n Particles: {} \n'.format(files[0]))
    
    # Loop over the particle files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff,inptype)
        p0 = f[inptype]  

        # Read particle information
        if (iff == 0):
            groupnum = p0['GroupNumber'][:] # FoF group number particle is in
            # Negative values: particles within r200 but not part of the halo
            subgroupnum = p0['SubGroupNumber'][:]
            partmass = p0['Mass'][:]            # 1e10 Msun/h
            partx = p0['Coordinates'][:,0]      # Mpc/h
            party = p0['Coordinates'][:,1]
            partz = p0['Coordinates'][:,2] 
        else:
            groupnum    = np.append(groupnum,p0['GroupNumber'][:])
            subgroupnum = np.append(subgroupnum,p0['SubGroupNumber'][:])
            partmass    = np.append(partmass,p0['Mass'][:])
            partx       = np.append(partx,p0['Coordinates'][:,0])
            party       = np.append(party,p0['Coordinates'][:,0])
            partz       = np.append(partz,p0['Coordinates'][:,0])

    # If all groupnum are less than 0, take abs()
    allgneg = False
    ind = np.where(groupnum<0)
    if(np.shape(ind)[1] == len(groupnum)):
        allgneg = True
        groupnum = abs(groupnum)-1
        if verbose: print('All {}/GroupNumber < 0'.format(inptype))
        
    # Get particle information into a pandas dataset to facilitate merging options
    #here: This operation changes groupnum and subgroupnum into floats, but doesn't seem to matter
    df_part = pd.DataFrame(data=np.vstack([groupnum,subgroupnum,partmass,partx,party,partz]).T,
                           columns=['groupnum','subgroupnum','partmass','partx','party','partz'])
    groupnum,subgroupnum,partmass,partx,party,partz=[[] for i in range(6)] #Empty individual arrays
    df_part.sort_values(by=['groupnum', 'subgroupnum'], inplace=True)
    df_part.reset_index(inplace=True, drop=True)  

    # Get FOF/Subfind files
    files, allfiles = get_subfind_files(snap,sim,env)
    if (not allfiles):
        print('\n WARNING (bahamas.map_m500): no adequate Subfind files found, {}, {}'.
              format(snap,env))
        return None
    if (Testing): files = [files[0],files[1]]
    if verbose: print('\n Haloes (FoF): {} \n'.format(files[0]))
    
    # Loop over the FOF/Subfind files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff)
        fof = f['FOF']

        # Read halo information
        if (iff == 0):
            m500  = fof['Group_M_Crit500'][:]           #1e10Msun/h
            r500  = fof['Group_R_Crit500'][:]           #cMpc/h
            cop_x = fof['GroupCentreOfPotential'][:,0]  #cMpc/h
            cop_y = fof['GroupCentreOfPotential'][:,1]  #cMpc/h
            cop_z = fof['GroupCentreOfPotential'][:,2]  #cMpc/h
        else:
            m500  = np.append(m500,fof['Group_M_Crit500'][:])
            r500  = np.append(r500,fof['Group_R_Crit500'][:])
            cop_x  = np.append(cop_x,fof['GroupCentreOfPotential'][:,0])
            cop_y  = np.append(cop_y,fof['GroupCentreOfPotential'][:,1])
            cop_z  = np.append(cop_z,fof['GroupCentreOfPotential'][:,2])

    # Get the FOF data into a pandas dataset
    df_fof = pd.DataFrame(data=np.vstack([m500,r500,cop_x,cop_y,cop_z]).T,
                          columns=['m500','r500','cop_x','cop_y','cop_z'])
    m500,r500,cop_x,cop_y,cop_z=[[] for i in range(5)] #Empty individual arrays

    # Generate a groupnum column from 1 to the last halo
    df_fof.index += 1
    df_fof.index.names = ['groupnum'] 
    df_fof.reset_index(inplace=True)

    # Join the particle and FoF information
    merge = pd.merge(df_part, df_fof, on=['groupnum'])

    # Get the boxsize
    omega0, omegab, lambda0, h0, boxsize = get_cosmology(sim,env)
    lbox2 = boxsize/2.

    # Position of particles relative to the center of the group
    merge['partx'] = merge.partx - merge.cop_x
    merge['party'] = merge.party - merge.cop_y
    merge['partz'] = merge.partz - merge.cop_z

    # Correct for periodic boundary conditions (for gal. in groups)
    merge.partx.loc[merge.partx < -lbox2] = merge.partx.loc[merge.partx < -lbox2] + boxsize
    merge.party.loc[merge.party < -lbox2] = merge.party.loc[merge.party < -lbox2] + boxsize
    merge.partz.loc[merge.partz < -lbox2] = merge.partz.loc[merge.partz < -lbox2] + boxsize

    merge.partx.loc[merge.partx >= lbox2] = merge.partx.loc[merge.partx >= lbox2] - boxsize
    merge.party.loc[merge.party >= lbox2] = merge.party.loc[merge.party >= lbox2] - boxsize
    merge.partz.loc[merge.partz >= lbox2] = merge.partz.loc[merge.partz >= lbox2] - boxsize

    # Distances to selected particles
    merge = merge.loc[merge.m500 > mlim]
    merge['distance'] = (merge.partx**2 +     
                         merge.party**2 +
                         merge.partz**2) ** 0.5

    # Mass of those particles enclosed in R500
    merge['inside_r500'] = merge.distance <= merge.r500
    merge = merge.loc[merge.inside_r500 == True]

    groups = merge.groupby(['groupnum'], as_index=False)
    massinr500 = groups.partmass.sum() # partmass now = particle mass (1e10 Msun/h)
    final = pd.merge(massinr500, df_fof, on=['groupnum'])
    final.partmass = np.log10(final.partmass) + 10. #log10(M/Msun/h)
    final.m500     = np.log10(final.m500) + 10.     #log10(M/Msun/h)

    # Write properties to output file
    if (not file_exists):
        # Generate the file
        hf = h5py.File(outfile, 'w') #here if appending columns here it'll be the place
        
        # Output header
        headnom = 'header'
        head = hf.create_dataset(headnom,(100,))
        head.attrs[u'sim']          = sim
        head.attrs[u'snapshot']     = snap
        head.attrs[u'redshift']     = get_z(snap,sim,env,dirz=dirz)
        head.attrs[u'omega0']       = omega0
        head.attrs[u'omegab']       = omegab
        head.attrs[u'lambda0']      = lambda0        
        head.attrs[u'h0']           = h0
        head.attrs[u'boxsize']      = boxsize

        # Output data with units
        hfdat = hf.create_group('data')
        
        prop = final[['cop_x', 'cop_y', 'cop_z']].to_numpy()
        hfdat.create_dataset('pos',data=prop); prop = []
        hfdat['pos'].dims[0].label = 'x,y,z (Mpc/h)'

        prop = final[['groupnum']].to_numpy()
        hfdat.create_dataset('groupnum',data=prop); prop = []
        hfdat['groupnum'].dims[0].label = 'FoF group number' 

        prop = final[['m500']].to_numpy()
        hfdat.create_dataset('m500',data=prop); prop = []
        hfdat['m500'].dims[0].label = 'log10(M/Msun/h)' 

        prop = final[['r500']].to_numpy()
        hfdat.create_dataset('r500',data=prop); prop = []
        hfdat['r500'].dims[0].label = 'cMpc/h'
        
        prop = final[['partmass']].to_numpy()
        hfdat.create_dataset(nompartmass,data=prop); prop = []
        hfdat[nompartmass].dims[0].label = 'log10(M/Msun/h)' 

        hf.close()
    elif (file_exists):
        # Open output file to append dataset
        hf = h5py.File(outfile, 'a') 
        prop = final[['partmass']].to_numpy()
        hf.create_dataset('data/'+nompartmass,data=prop); prop = []
        hf['data/'+nompartmass].dims[0].label = 'log10(M/Msun/h)' 
        hf.close()

    # Retrurn name of file with output
    return outfile


def get_mHMRmap_file(outdir,sim,snap,nhmr=2.,cop=False):
    '''
    Get the name and existance check of the map_HMR file

    Parameters
    ----------
    outdir: string
       Directory to write or find the file
    sim: string
       Simulation name or path
    snap: string
       Snapshot of the simulation
    nhrm: float
       Times the HalfMassRadius is considered
    cop: boolean
       If True, using CentreOfPotential, otherwise CentreOfMass 
    
    Returns
    -------
    outfile: string
       Name of the map_HMR file
    file_exists: boolean
       True if file exists
    '''
    
    if (sim==None):
        outdir2 = outdir
    else:
        outdir2 = outdir+sim+'/'
    dir_exists = io.create_dir(outdir2)
    snhmr = ('%f' % nhmr).rstrip('0').rstrip('.').replace('.','_')
    if cop:
        outfile = outdir2+'m'+snhmr+'HMRmap_cop_snap'+str(snap)+'.hdf5'
    else:
        outfile = outdir2+'m'+snhmr+'HMRmap_snap'+str(snap)+'.hdf5'
    file_exists = io.check_file(outfile)

    return outfile, file_exists


def map_mHMR(snap,sim,env,ptype='BH',mlim=0.,nhmr=2.,cop=False,
             dirz=None,outdir=None,Testing=True,verbose=False):
    '''
    Map particle mass into the half mass radius (HMR) of (central) subhaloes
    Within 'arilega' there is not enough information to map on satellite subhaloes.

    Parameters
    -----------
    snap : int
        Snapshot number
    sim : string
        Name of the simulation
    env : string
        ari, arilega or cosma, to use the adecuate paths
    ptype : string
        Name of one of the allowed ptypes, 0:gas, 4: stars, 5:BH
    mlim : float
        mass limit for subhaloes to be considered, M_30kp [Msun/h] 
    nhmr : float
        Enclosure radius = nhmr*HalfMassRadius(DM)
    cop  : boolean
        True to use the CentreOfPotential, False for CentreOfMass
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    Testing : boolean
        Calculations on part or all the simulation
    verbose : boolean
        To output extra information

    Returns
    -----
    prop : float array
        Mapped property

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.map_mHMR(31,sim,'arilega',ptype='BH')
    '''

    # Stop for environments different to arilega
    if (env != 'arilega'):
        print('STOP: Function bahamas.map_mHMR developed for env=arilega.')
        return None

    # Type of particles to be read
    itype = ptypes.index(ptype) # 0:gas, 1:DM, 4: stars, 5:BH
    inptype = 'PartType'+str(itype)
    nompartmass = 'mHMR_'+ptype
    
    # Output file
    outfile, file_exists = get_mHMRmap_file(outdir,sim,snap,nhmr,cop)
    
    # Get particle files
    files, allfiles = get_particle_files(snap,sim,env)
    if (not allfiles):
        print('\n WARNING (bahamas.map_mHMR): no adequate particle files found, {}, {}'.
              format(snap,env))
        return None
    if Testing: files = [files[0]]

    # Loop over the particle files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff,inptype)
        p0 = f[inptype]  

        # Read particle information
        if (iff == 0):
            groupnum = p0['GroupNumber'][:] # FoF group number particle is in
            # Negative values: particles within r200 but not part of the halo
            subgroupnum = p0['SubGroupNumber'][:]
            partmass = p0['Mass'][:]            # 1e10 Msun/h
            partx = p0['Coordinates'][:,0]      # Mpc/h
            party = p0['Coordinates'][:,1]
            partz = p0['Coordinates'][:,2] 
        else:
            groupnum    = np.append(groupnum,p0['GroupNumber'][:])
            subgroupnum = np.append(subgroupnum,p0['SubGroupNumber'][:])
            partmass    = np.append(partmass,p0['Mass'][:])
            partx       = np.append(partx,p0['Coordinates'][:,0])
            party       = np.append(party,p0['Coordinates'][:,0])
            partz       = np.append(partz,p0['Coordinates'][:,0])

    # If all groupnum are less than 0, take abs()
    allgneg = False
    ind = np.where(groupnum<0)
    if(np.shape(ind)[1] == len(groupnum)):
        allgneg = True
        groupnum = abs(groupnum)-1

    # Get particle information into a pandas dataset to facilitate merging options
    #here: This operation changes groupnum and subgroupnum into floats, but doesn't seem to matter
    df_part = pd.DataFrame(data=np.vstack([groupnum,subgroupnum,partmass,partx,party,partz]).T,
                           columns=['groupnum','subgroupnum','partmass','partx','party','partz'])
    groupnum,subgroupnum,partmass,partx,party,partz=[[] for i in range(6)] #Empty individual arrays
    df_part.sort_values(by=['groupnum', 'subgroupnum'], inplace=True)
    df_part.reset_index(inplace=True, drop=True)  

    # Get FOF&Subfind files
    files, allfiles = get_subfind_files(snap,sim,env)
    if (not allfiles):
        print('\n WARNING (bahamas.map_mHMR): no adequate Subfind files found, {}, {}'.
              format(snap,env))
        return None
    if Testing: files = [files[0],files[1]]

    # Prop index
    stype = ptypes.index('star')
    dmtype = ptypes.index('DM')

    # Loop over the FOF&Subfind files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff)
        sh = f['Subhalo']

        # Read halo information
        if (iff == 0):
            groupnum  = sh['GroupNumber'][:]      #FOF GroupNumber
            ms30  = sh['Mass_030kpc'][:,stype]    #1e10Msun/h
            HMRdm = sh['HalfMassRad'][:,dmtype]   #cMpc/h
            if cop:
                cop_x = sh['CentreOfPotential'][:,0]  #cMpc/h
                cop_y = sh['CentreOfPotential'][:,1]  #cMpc/h
                cop_z = sh['CentreOfPotential'][:,2]  #cMpc/h
            else:
                cop_x = sh['CentreOfMass'][:,0]  #cMpc/h
                cop_y = sh['CentreOfMass'][:,1]  #cMpc/h
                cop_z = sh['CentreOfMass'][:,2]  #cMpc/h
        else:
            groupnum  = np.append(groupnum,sh['GroupNumber'][:])
            ms30  = np.append(ms30,sh['Mass_030kpc'][:,stype])
            HMRdm = np.append(HMRdm,sh['HalfMassRad'][:,dmtype])
            if cop:
                cop_x = np.append(cop_x,sh['CentreOfPotential'][:,0])
                cop_y = np.append(cop_y,sh['CentreOfPotential'][:,1])
                cop_z = np.append(cop_z,sh['CentreOfPotential'][:,2])
            else:
                cop_x = np.append(cop_x,sh['CentreOfMass'][:,0])
                cop_y = np.append(cop_y,sh['CentreOfMass'][:,1])
                cop_z = np.append(cop_z,sh['CentreOfMass'][:,2])

    if verbose: print('All read galaxies = {:d}'.format(len(groupnum)))
            
    # Get indexes for centrals
    cind = get_cenids(snap,sim,env)
    if (max(cind) > len(groupnum) and Testing):
        ind = np.where(cind < len(groupnum)-1) 
        if(np.shape(ind)[1]<1):
            print('STOP (bahamas.map_mHMR): no centrals in Subfind file.')
            return None
        cind = cind[ind]
    elif (max(cind) > len(groupnum) and not Testing):
        print('STOP (bahamas.map_mHMR): problem with centrals indexes.')
        return None
    
    if verbose:
        print('Number of central galaxies = {:d}'.format(len(cind)))
        print('Min. HMR = {:.3f} Mpc/h; Max. = {:.3f} Mpc/h'.format(
            min(HMRdm[cind]),max(HMRdm[cind])))
        print('Min. n*HMR = {:.3f} Mpc/h; Max. = {:.3f} Mpc/h'.format(
            min(nhmr*HMRdm[cind]),max(nhmr*HMRdm[cind])))
    
    # Mapping for central galaxies with stellar mass
    data = np.vstack([groupnum[cind],ms30[cind],HMRdm[cind],cop_x[cind],cop_y[cind],cop_z[cind]]).T
    df_sh = pd.DataFrame(data=data,
                         columns=['groupnum','ms30','HMRdm','cop_x','cop_y','cop_z'])
    data,groupnum,ms30,HMRdm,cop_x,cop_y,cop_z=[[] for i in range(7)] #Empty individual arrays
    df_sh = df_sh.loc[df_sh.ms30 > mlim] # With stellar mass
    if df_sh.empty:
        print('STOP (bahamas.map_mHMR): no centrals with stellar mass.')
        return None
    df_sh.ms30 = np.log10(df_sh.ms30) + 10.    #log10(M/Msun/h)
    
    # Join the particle and FoF information
    merge = pd.merge(df_part, df_sh, on=['groupnum'])
    del df_part

    # Get the boxsize
    omega0, omegab, lambda0, h0, boxsize = get_cosmology(sim,env)
    lbox2 = boxsize/2.

    # Position of particles relative to the center of the group
    merge['partx'] = merge.partx - merge.cop_x
    merge['party'] = merge.party - merge.cop_y
    merge['partz'] = merge.partz - merge.cop_z

    # Correct for periodic boundary conditions (for gal. in groups)
    merge.partx.loc[merge.partx < -lbox2] = merge.partx.loc[merge.partx < -lbox2] + boxsize
    merge.party.loc[merge.party < -lbox2] = merge.party.loc[merge.party < -lbox2] + boxsize
    merge.partz.loc[merge.partz < -lbox2] = merge.partz.loc[merge.partz < -lbox2] + boxsize

    merge.partx.loc[merge.partx >= lbox2] = merge.partx.loc[merge.partx >= lbox2] - boxsize
    merge.party.loc[merge.party >= lbox2] = merge.party.loc[merge.party >= lbox2] - boxsize
    merge.partz.loc[merge.partz >= lbox2] = merge.partz.loc[merge.partz >= lbox2] - boxsize

    # Distances to selected particles
    merge['distance'] = (merge.partx**2 +     
                         merge.party**2 +
                         merge.partz**2) ** 0.5
    if verbose: print('Min. distance to centre = {:.3f} Mpc/h; Max. = {:.3f} Mpc/h'.format(
            merge['distance'].min(),merge['distance'].max()))
    
    # Mass of those particles enclosed in radius
    radius = nhmr*merge.HMRdm
    merge['inside_HMRdm'] = merge.distance <= radius
    merge = merge.loc[merge.inside_HMRdm == True]
    if merge.empty:
        print('STOP (bahamas.map_mHMR): no particles within DM HMR.')
        return None

    groups = merge.groupby(['groupnum'], as_index=False)
    massinHMRdm = groups.partmass.sum() # partmass now = particle mass (1e10 Msun/h)

    final = pd.merge(massinHMRdm, df_sh, on=['groupnum'])
    final.partmass = np.log10(final.partmass) + 10. #log10(M/Msun/h)
    if verbose: print(final)
    
    # Write properties to output file        
    hf = h5py.File(outfile, 'w') # Generate the file
    
    # Output header
    headnom = 'header'
    head = hf.create_dataset(headnom,(100,))
    head.attrs[u'sim']          = sim
    head.attrs[u'snapshot']     = snap
    head.attrs[u'redshift']     = get_z(snap,sim,env,dirz=dirz)
    head.attrs[u'omega0']       = omega0
    head.attrs[u'omegab']       = omegab
    head.attrs[u'lambda0']      = lambda0        
    head.attrs[u'h0']           = h0
    head.attrs[u'boxsize']      = boxsize

    # Output data with units
    hfdat = hf.create_group('data')
    
    prop = final[['cop_x', 'cop_y', 'cop_z']].to_numpy()
    hfdat.create_dataset('pos',data=prop); prop = []
    hfdat['pos'].dims[0].label = 'x,y,z (Mpc/h)'

    prop = final[['groupnum']].to_numpy()
    hfdat.create_dataset('groupnum',data=prop); prop = []
    hfdat['groupnum'].dims[0].label = 'FoF group number' 

    prop = final[['ms30']].to_numpy()
    hfdat.create_dataset('ms30',data=prop); prop = []
    hfdat['ms30'].dims[0].label = 'log10(M/Msun/h)' 

    prop = final[['HMRdm']].to_numpy()
    hfdat.create_dataset('HMRdm',data=prop); prop = []
    hfdat['HMRdm'].dims[0].label = 'cMpc/h'
    
    prop = final[['partmass']].to_numpy()
    hfdat.create_dataset(nompartmass,data=prop); prop = []
    hfdat[nompartmass].dims[0].label = 'log10(M/Msun/h)' 

    hf.close()

    # Retrurn name of file with output
    return outfile


def old_get_subBH_file(outdir,sim,snap,part=False,addp=False,nhmr=2.,cop=True):
    '''
    TO BE REMOVED once test plots are updated
    Get the name and existance check of the files generated with either
    get_subBH (particle info.) or map_subBH (mapped particles)

    Parameters
    ----------
    outdir: string
       Directory to write or find the file
    sim: string
       Simulation name or path
    snap: string
       Snapshot of the simulation
    part : boolean
        True for the particle file name; False for the mapped particles into haloes file 
    addp : boolean
        True for added particle information.
    nhrm: float
       Times the HalfMassRadius is considered
    cop: boolean
       If True, using CentreOfPotential, otherwise CentreOfMass
    
    Returns
    -------
    outfile: string
       Name of the map_HMR file
    file_exists: boolean
       True if file exists
    '''
    
    if (sim==None):
        outdir2 = outdir
    else:
        outdir2 = outdir+sim+'/'
    dir_exists = io.create_dir(outdir2)

    if part:
        outfile = outdir2+'part_subBH_snap'+str(snap)+'.hdf5'
        if addp:
            outfile = outdir2+'part_add_subBH_snap'+str(snap)+'.hdf5'
    else:
        snhmr = ('%f' % nhmr).rstrip('0').rstrip('.').replace('.','_')
        if cop:
            outfile = outdir2+'subBH_'+snhmr+'HMRmap_cop_snap'+str(snap)+'.hdf5'
        else:
            outfile = outdir2+'subBH_'+snhmr+'HMRmap_com_snap'+str(snap)+'.hdf5'

    file_exists = io.check_file(outfile)

    return outfile, file_exists



def old_get_subBH(snap,sim,env,addp=False,dirz=None,outdir=None,Testing=True,verbose=False):
    '''
    TO BE REMOVED once testing plots are updated
    Produce a file joining subgrid BH properties with their positions 
    and information on halo identifier. Joining done on PartID.
    There is an option to add the masses and accretion for particles in the same position.

    Parameters
    -----------
    snap : int
        Snapshot number
    sim : string
        Name of the simulation
    env : string
        ari, arilega or cosma to use the adecuate paths
    addp : boolean
        True to add properties of particles in the same position.
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    Testing : boolean
        Calculations on part or all the simulation
    verbose : boolean
        To output extra information

    Returns
    -----
    outfile : string
        Name of the output file

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.get_subBH(31,sim,'arilega')
    '''

    # Stop for environments different to arilega
    if (env != 'arilega'):
        print('\n WARNING (b.get_subBH): Function developed for env=arilega.')
        return None

    # Black hole particles to be read, 5:BH
    itype = 5 
    inptype = 'PartType'+str(itype)

    # Output file
    outfile, file_exists = old_get_subBH_file(outdir,sim,snap,part=True,addp=addp)

    # Get subgrid particle information from snapshots------------------------
    files, allfiles = get_particle_files(snap,sim,env,subfind=False)
    if (not allfiles):
        print('\n WARNING (b.get_subBH): no adequate particle files found, {}, {}'.
              format(snap,env))
        return None
    if Testing: files = [files[0],files[1]]
    if verbose: print('Particles: {} \n'.format(files[0]))

    # Loop over the particle files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff,inptype)
        p0 = f[inptype]

        # Read particle information
        if (iff == 0):
            # Check that there is data to be read
            try:
                partID  = p0['ParticleIDs'][:]
            except:
                print('\n WARNING (b.get_subBH): empty data {}'.format(ff+'/'+inptype+'/ParticleIDs'))
                return None

            # Read the data
            BH_Mass = p0['BH_Mass'][:]  # 1e10 Msun/h
            BH_Mdot = p0['BH_Mdot'][:]*10**7*ast.s_in_year/ast.m_in_pc  # Msun/year
        else:
            partID  = np.append(partID,p0['ParticleIDs'][:]) 
            BH_Mass = np.append(BH_Mass,p0['BH_Mass'][:])
            BH_Mdot = np.append(BH_Mdot,p0['BH_Mdot'][:]*10**7*ast.s_in_year/ast.m_in_pc)
            
    if verbose:
        print('BH: seed={:.2e}; min={:.2e}, max={:.2e}'.format(BH_seed_mass,
                                                               min(BH_Mass)*10**10,
                                                               max(BH_Mass)*10**10))

    # Get subgrid information into a pandas dataset to facilitate merging options
    data = np.vstack([partID,BH_Mass,BH_Mdot]).T
    df_psnap = pd.DataFrame(data=data,columns=['partID','BH_Mass','BH_Mdot'])
    partID,BH_Mass,BH_Mdot=[[] for i in range(3)] #Empty individual arrays    
    
    # Get Subfind particle files----------------------------------------------
    files, allfiles = get_particle_files(snap,sim,env)
    if (not allfiles):
        print('\n WARNING (b.get_subBH): no adequate particle files found, {}, {}'.
              format(snap,env))
        return None
    if Testing: files = [files[0]]
    if verbose: print('\n Subfind particles: {} \n'.format(files[0]))

    # Loop over the particle files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff,inptype)
        p0 = f[inptype]

        # Read particle information
        if (iff == 0):
            partID = p0['ParticleIDs'][:] 
            groupnum = p0['GroupNumber'][:] # FoF group number particle is in
            # Negative values: particles within r200 but not part of the halo
            subgroupnum = p0['SubGroupNumber'][:]
            partx = p0['Coordinates'][:,0]      # Mpc/h
            party = p0['Coordinates'][:,1]
            partz = p0['Coordinates'][:,2] 
            pvx = p0['Velocity'][:,0]           # km/s
            pvy = p0['Velocity'][:,1]
            pvz = p0['Velocity'][:,2] 
        else:
            partID      = np.append(partID,p0['ParticleIDs'][:]) 
            groupnum    = np.append(groupnum,p0['GroupNumber'][:])
            subgroupnum = np.append(subgroupnum,p0['SubGroupNumber'][:])
            partx       = np.append(partx,p0['Coordinates'][:,0])
            party       = np.append(party,p0['Coordinates'][:,1])
            partz       = np.append(partz,p0['Coordinates'][:,2])
            pvx         = np.append(pvx,p0['Velocity'][:,0])
            pvy         = np.append(pvy,p0['Velocity'][:,1])
            pvz         = np.append(pvz,p0['Velocity'][:,2])

    if verbose:
        print('x (Mpc/h): min={:.2f}, max={:.2f}'.format(min(partx),max(partx)))
        print('vx (km/s): min={:.2f}, max={:.2f}'.format(min(pvy),max(pvy)))
        print('GroupNum: min={:d}, max={:d}'.format(min(groupnum),max(groupnum)))
        print('SubGroupNum: min={:d}, max={:d} (diff={:d})'.format(min(subgroupnum),max(subgroupnum),
                                                                   min(subgroupnum)-max(subgroupnum)))

    # If all groupnum are less than 0, take abs()
    allgneg = False
    ind = np.where(groupnum<0)
    if(np.shape(ind)[1] == len(groupnum)):
        allgneg = True
        groupnum = abs(groupnum)-1

    # Get particle information into a pandas dataset to facilitate merging options
    #here: This operation changes groupnum and subgroupnum into floats, but doesn't seem to matter
    data = np.vstack([partID,groupnum,subgroupnum,partx,party,partz,pvx,pvy,pvz]).T 
    df_psub = pd.DataFrame(data=data,columns=['partID','groupnum','subgroupnum',
                                              'partx','party','partz',
                                              'pvx','pvy','pvz'])
    partID,groupnum,subgroupnum,partx,party,partz,pvx,pvy,pvz=[[] for i in range(9)]

    # Join the particle information---------------------------------------------
    df_part = pd.merge(df_psub, df_psnap, on=['partID'])
    df_part.sort_values(by=['groupnum', 'subgroupnum'], inplace=True)
    #df_part.sort_values(by=['groupnum', 'partx'], inplace=True)
    df_part.reset_index(inplace=True, drop=True)  

    # Add properties of particles in the same position--------------------------
    if (addp):
        # Group particles in the same position
        groups = df_part.groupby(['groupnum','partx','party','partz'], as_index=False)
        if(groups.ngroups > len(df_part.index)):
            print('\n WARNING (b.get_subBH): problem grouping particle information, {}, {}'.
              format(snap,sim))
            return None

        # Number of particles with the same position
        df_nboson = groups.size().reset_index(name='nboson')
        if verbose: print('Max. BH bosons = {}'.format(df_nboson['nboson'].max()))

        # Number of particles with the same position and velocity
        groupv = df_part.groupby(['groupnum','partx','party','partz','pvx','pvy','pvz'], as_index=False)
        if(groupv.ngroups > len(df_part.index)):
            print('\n WARNING (b.get_subBH): problem grouping particle information, {}, {}'.
              format(snap,sim))
            return None
        df_nbosonv = groupv.size().reset_index(name='nbosonv')
        if verbose: print('Max. BH w equal pos. and vel. = {}'.format(df_nbosonv['nbosonv'].max()))

        # Add the properties for particles in the same position
        df_addM_BH = groups.BH_Mass.sum() # 1e10 Msun/h
        df_addMdot = groups.BH_Mdot.sum() # Msun/year
        df_addM1 = pd.merge(df_addM_BH, df_addMdot, on=['groupnum','partx','party','partz'])
        del df_addM_BH, df_addMdot

        df_addM = pd.merge(df_addM1, df_nboson, on=['groupnum','partx','party','partz'])
        del df_addM1, df_nboson

        # Remove duplicated columns and row from initial particle information
        df1 = df_part[['groupnum','partx','party','partz','partID','pvx','pvy','pvz','subgroupnum']]
        df2 = df1.drop_duplicates(subset=['groupnum','partx','party','partz'],
                                  keep='last',ignore_index=True)
        df3 = df2.sort_values(by=['groupnum','partx','party','partz'], ignore_index=True)
        del df_part, df1, df2

        # Generate the final data set with the merge
        final = pd.merge(df3, df_addM, on=['groupnum','partx','party','partz'])
        del df3, df_addM
    else:
        final = df_part
        del df_part
    if verbose: print(final)

    # Write output file---------------------------------------------------------    
    hf = h5py.File(outfile, 'w') # Generate the file

    # Get the cosmological parameters and boxsize
    omega0, omegab, lambda0, h0, boxsize = get_cosmology(sim,env)
    
    # Output header
    headnom = 'header'
    head = hf.create_dataset(headnom,(100,))
    head.attrs[u'sim']          = sim
    head.attrs[u'snapshot']     = snap
    head.attrs[u'redshift']     = get_z(snap,sim,env,dirz=dirz)
    head.attrs[u'omega0']       = omega0
    head.attrs[u'omegab']       = omegab
    head.attrs[u'lambda0']      = lambda0        
    head.attrs[u'h0']           = h0
    head.attrs[u'boxsize']      = boxsize

    # Output data with units
    hfdat = hf.create_group('data')

    prop = final[['partID']].to_numpy()
    hfdat.create_dataset('partID',data=prop); prop = []
    hfdat['partID'].dims[0].label = 'Particle ID' 
    
    prop = final[['partx', 'party', 'partz']].to_numpy()
    hfdat.create_dataset('pos',data=prop); prop = []
    hfdat['pos'].dims[0].label = 'x,y,z (Mpc/h)'

    prop = final[['groupnum']].to_numpy()
    hfdat.create_dataset('groupnum',data=prop); prop = []
    hfdat['groupnum'].dims[0].label = 'FoF group number' 

    prop = final[['subgroupnum']].to_numpy()
    hfdat.create_dataset('subgroupnum',data=prop); prop = []
    hfdat['subgroupnum'].dims[0].label = 'Subgroup number particle is in' 
    
    prop = final[['BH_Mass']].to_numpy()
    hfdat.create_dataset('BH_Mass',data=prop); prop = []
    hfdat['BH_Mass'].dims[0].label = '10**10 Msun/h' 

    prop = final[['BH_Mdot']].to_numpy()
    hfdat.create_dataset('BH_Mdot',data=prop); prop = []
    hfdat['BH_Mdot'].dims[0].label = 'Msun/year' 

    if (addp):
        prop = final[['nboson']].to_numpy()
        hfdat.create_dataset('nboson',data=prop); prop = []
        hfdat['nboson'].dims[0].label = 'Number of BH particles at the same position.'     
    
    hf.close()
    
    # Retrurn name of file with output
    return outfile


def get_subBH_file(outdir,sim,snap,part=True,nhmr=2.,cop=True):
    '''
    Get the name and existance check of the files generated with either
    get_subBH (particle info.) or map_subBH (mapped particles)

    Parameters
    ----------
    outdir: string
       Directory to write or find the file
    sim: string
       Simulation name or path
    snap: string
       Snapshot of the simulation
    part : boolean
        True for the particle file name; False for the mapped particles into haloes file 
    nhrm: float
       Times the HalfMassRadius is considered
    cop: boolean
       If True, using CentreOfPotential, otherwise CentreOfMass
    
    Returns
    -------
    outfile: string
       Name of the map_HMR file
    file_exists: boolean
       True if file exists
    '''
    
    if (sim==None):
        outdir2 = outdir
    else:
        outdir2 = outdir+sim+'/'
    dir_exists = io.create_dir(outdir2)

    if part:
        outfile = outdir2+'subBH_part_snap'+str(snap)+'.hdf5'
    else:
        snhmr = ('%f' % nhmr).rstrip('0').rstrip('.').replace('.','_')
        if cop:
            outfile = outdir2+'subBH_'+snhmr+'HMRmap_cop_snap'+str(snap)+'.hdf5'
        else:
            outfile = outdir2+'subBH_'+snhmr+'HMRmap_com_snap'+str(snap)+'.hdf5'

    file_exists = io.check_file(outfile)

    return outfile, file_exists


def get_subnum(sim,snap,env,Testing=False,nfiles=2,verbose=False):
    """
    Create an array with indexes 0 (central) to N satellites, 
    within each halo, given by the "groupnum"

    Parameters
    -----------
    sim : string
        Simulation name
    snap : integer
        Snapshot 
    env : string
        ari or cosma, to use the adecuate paths
    Testing: boolean
        True or False
    nfiles : integer
        Number of files to be considered for testing
    verbose : boolean
        True to write first Subfind file out

    Returns
    -----
    subnum : array of (long) integers
        Indexes 0 (central) to N galaxies within each subhalo

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_subnum('HIRES/AGN_TUNED_nu0_L050N256_WMAP9',31,'arilega',Testing=True)
    """

    # Read the FOF Group Number subhalo belongs to and its mass (1e10Msun/h)
    groupnum = get_subfind_prop(snap,sim,env,'Subhalo/GroupNumber',
                                Testing=Testing,verbose=verbose)
    msubh    = get_subfind_prop(snap,sim,env,'Subhalo/Mass',
                                Testing=Testing,verbose=verbose)
    #msubh    = get_subfind_prop(snap,sim,env,'Subhalo/Mass_030kpc',proptype='DM',Testing=Testing)
    #mh = get_subfind_prop(snap,sim,env,'FOF/Group_M_Crit200',Testing=Testing)

    # Index of first subhalo in SubHalo list (starts at 0) 
    cind = get_cenids(snap,sim,env,Testing=Testing)

    if (max(cind) > (len(groupnum)-1) and Testing):
        ind = np.where(cind < len(groupnum)-1) 
        if(np.shape(ind)[1]<1):
            print('STOP (b.get_subnum): no centrals in Subfind file.')
            return None
        cind = cind[ind]
    elif (max(cind) > (len(groupnum)-1) and not Testing):
        print('STOP (b.get_subnum): problem with centrals indexes.')
        return None

    # Initialize the index for subhaloes and store a 0 for centrals
    subgroupnum = np.zeros(len(groupnum),dtype=int); subgroupnum.fill(np.nan)
    subgroupnum[cind] = 0
    print(subgroupnum[0:3],subgroupnum[-10:])
    print(groupnum[0:3],groupnum[-10:])
    print(print(len(cind),min(cind),max(cind))); exit() ###here

    #    # Haloes with no substructure are not well defined
    #    nsubh = f['FOF/NumOfSubhalos'][:]
    #    ind = np.where(nsubh>0)
    #    if (np.shape(ind)[1]>0):
    
    inotord = 0; notord = []     
    for i, subg in enumerate(subgroupnum):
        if subg == -1:
            continue
        elif subg == 0:
            isub = 0;
            maxmass = msubh[i];
        else:
            isub += 1;
            subgroupnum[i] = isub;
            if (msubh[i]>maxmass):
                notord.append(groupnum[i])
                inotord += 1

    if (inotord>0 and verbose):
        print('\n WARNING b.subgroupnum: mass not ordered for {} haloes: {}'.format(inotord,notord))

    return subgroupnum


def get_subhalo4BH(snap,sim,env,dirz=None,outdir=None,rewrite=False,Testing=False,nfiles=2,verbose=False):
    """
    Check if a file with the subhalo properties relevant for the
    analysis of subgrid BH particles is there, and otherwise create it.

    Parameters
    -----------
    snap : integer
        Snapshot 
    sim : string
        Simulation name
    env : string
        ari, arilega or cosma to use the adecuate paths
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    rewrite: boolean
        True or False for rewriting the output file if it already exists
    Testing: boolean
        True or False
    nfiles : integer
        Number of files to be considered for testing
    verbose : boolean
        True to write first Subfind file out

    Returns
    -----
    outfile : string
        Name of the file with subhalo information
    file_exists: boolean
       True if file exists

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_subhalo4BH(27,'L400N1024/WMAP9/Sims/BAHAMAS')
    """

    # Output file
    outfile, file_exists = get_subBH_file(outdir,sim,snap,part=True)
    if (file_exists and not rewrite):
        return outfile, file_exists

    # Get halo information from FOF&Subfind files-------------------------------
    files, allfiles = get_subfind_files(snap,sim,env)
    if allfiles is False:
        print('\n WARNING (b.get_subhalo4BH): no adequate Subfind files found, {}, {}, {}'.
              format(files[0],snap,env))
        return outfile, False
    if verbose: print('b.get_subhalo4BH: First Subfind file is {}'.format(files[0]))

    # Prop index
    stype = ptypes.index('star')
    dmtype = ptypes.index('DM')
    
    # Loop over the FOF&Subfind files
    new_subnum = False
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff)
        fof = f['FOF']
        sh = f['Subhalo']

        # Read subhalo information
        if (iff == 0):
            mhalo = fof['Group_M_Crit200'][:]           # 1e10Msun/h 
            rhalo = fof['Group_R_Crit200'][:]           # cMpc/h 
            fof_x = fof['GroupCentreOfPotential'][:,0]  # cMpc/h 
            fof_y = fof['GroupCentreOfPotential'][:,1]  # cMpc/h
            fof_z = fof['GroupCentreOfPotential'][:,2]  # cMpc/h
            
            groupnum  = sh['GroupNumber'][:]      # FOF GroupNumber
            ms30  = sh['Mass_030kpc'][:,stype]    # 1e10Msun/h
            SFR   = sh['StarFormationRate'][:]    # Msun/h/yr
            cop_x = sh['CentreOfPotential'][:,0]  # cMpc/h
            cop_y = sh['CentreOfPotential'][:,1]  # cMpc/h
            cop_z = sh['CentreOfPotential'][:,2]  # cMpc/h
            shv_x = sh['Velocity'][:,0]           # km/s
            shv_y = sh['Velocity'][:,1]           # km/s
            shv_z = sh['Velocity'][:,2]           # km/s

            try:
                # SubGroup Number of subhalo, begins at 0 for most massive subhalo within a group 
                subgroupnum  = sh['SubGroupNumber'][:]
            except:
                new_subnum = True
        else:
            mhalo = np.append(mhalo,fof['Group_M_Crit200'][:])
            rhalo = np.append(rhalo,fof['Group_R_Crit200'][:])
            fof_x = np.append(fof_x,fof['GroupCentreOfPotential'][:,0])
            fof_y = np.append(fof_y,fof['GroupCentreOfPotential'][:,1])
            fof_z = np.append(fof_z,fof['GroupCentreOfPotential'][:,2])
            
            groupnum  = np.append(groupnum,sh['GroupNumber'][:])
            ms30  = np.append(ms30,sh['Mass_030kpc'][:,stype])
            SFR   = np.append(SFR,sh['StarFormationRate'][:])
            cop_x = np.append(cop_x,sh['CentreOfPotential'][:,0])
            cop_y = np.append(cop_y,sh['CentreOfPotential'][:,1])
            cop_z = np.append(cop_z,sh['CentreOfPotential'][:,2])
            shv_x = np.append(shv_x,sh['Velocity'][:,0])
            shv_y = np.append(shv_y,sh['Velocity'][:,1])
            shv_z = np.append(shv_z,sh['Velocity'][:,2])   

            if (not new_subnum):
                subgroupnum = np.append(subgroupnum,sh['SubGroupNumber'][:])

    if (not new_subnum):
        # Check the IDs of subhaloes make sense
        if (abs(max(subgroupnum)-min(subgroupnum))<1):
            new_subnum = True
    if (new_subnum):
        subgroupnum = get_subnum(sim,snap,env,verbose=verbose)
        if Testing and (len(subgroupnum)>len(groupnum)):
            subgroupnum = subgroupnum[0:len(groupnum)]
        if (len(groupnum) != len(subgroupnum)):
            print('STOP b.get_subhalo4BH: different lenghts for subgroupnum and groupnum, {}'.format(ff))

    # Only consider haloes with some mass enclosed in 30kpc
    ind = np.where(ms30 > 0.)
    if (np.shape(ind)[1] < 1):
        print('STOP (b.map_subBH): no centrals with stellar mass.')
        return None
    gn = groupnum[ind]

    cop_x = cop_x[ind]; cop_y = cop_y[ind]; cop_z = cop_z[ind]
    shv_x = shv_x[ind]; shv_y = shv_y[ind]; shv_z = shv_z[ind]    
    sgn = subgroupnum[ind]

    # Halo properties
    mh, rh = [np.zeros(shape=len(gn)) for i in range(2)]
    mh = mhalo[gn]
    rh = rhalo[gn]

    # Initialize relative distances and velocities
    dr, dvr, dvphi, dvlos = [np.zeros(shape=len(gn)) for i in range(4)]
    dvr.fill(np.nan)
    dvphi.fill(np.nan)
    dvlos.fill(np.nan)

    # Get the cosmological parameters and boxsize
    omega0, omegab, lambda0, h0, boxside = get_cosmology(sim,env)

    # Distance to halo center
    dr_fof = get_r(cop_x,cop_y,cop_z,
                   fof_x[gn],fof_y[gn],fof_z[gn],box=boxside)

    # Get indexes for centrals
    cind = np.where(sgn == 0)
    dr[cind] = dr_fof[cind]

    # Distance and velocities from satellite to central galaxies
    inotcen = 0; notcen = []     
    for ig in np.unique(gn):
        ii = np.where(gn == ig)
        if (np.shape(ii)[1] < 1): continue

        isats = np.where((gn == ig) & (sgn > 0))
        nsats = np.shape(isats)[1]
        if (nsats < 1): continue

        sx = cop_x[isats]; sy = cop_y[isats]; sz = cop_z[isats]
        svx = shv_x[isats]; svy = shv_y[isats]; svz = shv_z[isats]    

        icens = np.where((gn == ig) & (sgn <1))
        ncens = np.shape(icens)[1]
        if (ncens > 1):
            if verbose:
                print('\n WARNING (get_subhalo4BH): {} centrals in halo {}'.format(ncens,ig))
            continue
        else:
            cx,cy,cz,cvx,cvy,cvz = [np.zeros(shape=nsats,dtype=float) for i in range(6)]
            if (ncens < 1): # Distance to halo COP
                cx.fill(fof_x[ig]); cy.fill(fof_y[ig]); cz.fill(fof_z[ig])
                notcen.append(ig)
                inotcen += 1

            else: # Distance and velocities relative to central galaxy
                cx.fill(cop_x[icens][0]); cy.fill(cop_y[icens][0]); cz.fill(cop_z[icens][0])
                cvx.fill(shv_x[icens][0]); cvy.fill(shv_y[icens][0]); cvz.fill(shv_z[icens][0])

                # Radial velocity
                satd = get_vr(cx,cy,cz,sx,sy,sz,
                              cvx,cvy,cvz,svx,svy,svz,box=boxside)
                dvr[isats] = satd

                # Tangential velocity
                satd = get_vphi(cx,cy,cz,sx,sy,sz,
                                cvx,cvy,cvz,svx,svy,svz,box=boxside)
                dvphi[isats] = satd

                # Line of site (z-axis) velocity
                satd = get_vlos(cx,cy,cz,sx,sy,sz,
                                cvx,cvy,cvz,svx,svy,svz,box=boxside)
                dvlos[isats] = satd

        # Radial distance
        satd = get_r(cx,cy,cz,sx,sy,sz,box=boxside)
        dr[isats] = satd

    if (inotcen>0 and verbose):
        print('\n WARNING (get_subhalo4BH): no central selected in {} haloes, {}'.format(inotcen,notcen))

    # Save data in a dataframe
    data = np.vstack([gn,sgn,mh,rh,
                      dr,dvr,dvphi,dvlos,
                      cop_x,cop_y,cop_z,
                      shv_x,shv_y,shv_z,
                      ms30[ind],SFR[ind]]).T
    df_sh = pd.DataFrame(data=data,columns=['groupnum','subgroupnum',
                                            'M200C','R200C',
                                            'dr','dvr','dvphi','dvlos',
                                            'cop_x','cop_y','cop_z',
                                            'shv_x','shv_y','shv_z',
                                            'ms30','SFR'])
    data, gn,sgn,mh,rh,dr,dvr,dvphi,dvlos = [[] for i in range(9)]
    cop_x,cop_y,cop_z,shv_x,shv_y,shv_z,ms30,SFR = [[] for i in range(8)]
            
    # Write output file---------------------------------------------------------    
    hf = h5py.File(outfile, 'w') # Generate the file
    
    # Output header
    headnom = 'header'
    head = hf.create_dataset(headnom,(100,))
    head.attrs[u'sim']          = sim
    head.attrs[u'snapshot']     = snap
    head.attrs[u'redshift']     = get_z(snap,sim,env,dirz=dirz)
    head.attrs[u'omega0']       = omega0
    head.attrs[u'omegab']       = omegab
    head.attrs[u'lambda0']      = lambda0        
    head.attrs[u'h0']           = h0
    head.attrs[u'boxside']      = boxside
    head.attrs[u'new_subnum']   = new_subnum

    # Output data with units
    hfdat = hf.create_group('data')
    shdat = hfdat.create_group('Subhalo')

    # Subhalo/ groupnum, subnum, sat, pos, vel, M*30kpc, SFR
    #here include also dr, dv, dvr (from FOF)   
    noms = ['groupnum','subgroupnum',
            'M200C','R200C',
            'dr','dvr','dvphi','dvlos',
            'cop_x','cop_y','cop_z',
            'shv_x','shv_y','shv_z','ms30','SFR']
    desc = ['FoF group number','Subhalo index, cen:0',
            '1e10Msun/h, Mass within Rcrit200',
            'cMpc/h, Co-moving radius within which density is 200 times critical density',
            'Relative distance to FOF COP for centrals, and to central COP for sat. gal. (cMpc/h)',
            'Relative radial velocity to central gal. (km/s)',
            'Relative tangential (phi) velocity to central gal. (km/s)',
            'Relative velocity to central gal. on the line of sight, assumed z-axis (km/s)',
            'cMpc/h','cMpc/h','cMpc/h','km/s','km/s','km/s',
            '1e10 Msun/h','Msun/h/yr'] 
    for ip, iprop in enumerate(noms):
        nom = noms[ip]
        prop = df_sh[[nom]].to_numpy()
        shdat.create_dataset(nom,data=np.squeeze(prop)); prop = []
        shdat[nom].dims[0].label = desc[ip]
    
    hf.close()
    file_exists = io.check_file(outfile)

    return outfile, file_exists


def get_subBH(snap,sim,env,dirz=None,outdir=None,rewrite=True,Testing=True,verbose=False):
    '''
    Produce a file joining subgrid BH properties with their positions 
    and information on halo identifier. Joining done on PartID.

    Parameters
    -----------
    snap : int
        Snapshot number
    sim : string
        Name of the simulation
    env : string
        ari, arilega or cosma to use the adecuate paths
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    rewrite: boolean
        True or False for rewriting the output file if it already exists
    Testing : boolean
        Calculations on part or all the simulation
    verbose : boolean
        To output extra information

    Returns
    -----
    outfile : string
        Name of the output file

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.get_subBH(31,sim,'arilega')
    '''
    
    # Stop for environments different to arilega
    if (env != 'arilega'):
        print('\n WARNING (b.get_subBH): Function developed for env=arilega.')
        return None

    # Black hole particles to be read, 5:BH
    itype = 5 
    inptype = 'PartType'+str(itype)
    
    # Get subgrid particle information from snapshots------------------------
    files, allfiles = get_particle_files(snap,sim,env,subfind=False)
    if (not allfiles):
        print('\n WARNING (b.get_subBH): no adequate particle files found, {}, {}'.
              format(snap,env))
        return None
    if Testing: files = [files[0],files[1]]
    if verbose: print('Particles: {} \n'.format(files[0]))
    
    # File with information on subhaloes and to output BH information---------
    outfile, file_exists = get_subhalo4BH(snap,sim,env,dirz=dirz,outdir=outdir,
                                          rewrite=rewrite,Testing=Testing,verbose=verbose)
    if (not rewrite): return outfile
    if verbose: print('Outfile (Testing={}): {} \n'.format(Testing,outfile))

    f = h5py.File(outfile, 'r')
    new_subnum = f['header'].attrs['new_subnum']
    boxside = f['header'].attrs['boxside']
    sh = f['data/Subhalo/']
    gn    = sh['groupnum'][:]
    sgn   = sh['subgroupnum'][:]
    cop_x = sh['cop_x'][:]
    cop_y = sh['cop_y'][:]
    cop_z = sh['cop_z'][:]
    shv_x = sh['shv_x'][:]
    shv_y = sh['shv_y'][:]
    shv_z = sh['shv_z'][:]
    f.close()

    # Loop over the particle files-------------------------------------------
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff,inptype)
        p0 = f[inptype]
    
        # Read particle information
        if (iff == 0):
            # Check that there is data to be read
            try:
                partID  = p0['ParticleIDs'][:]
            except:
                print('\n WARNING (b.get_subBH): empty data {}'.format(ff+'/'+inptype+'/ParticleIDs'))
                return None
    
            # Read the data
            BH_Mass = p0['BH_Mass'][:]  # 1e10 Msun/h
            BH_Mdot = p0['BH_Mdot'][:]*10**7*ast.s_in_year/ast.m_in_pc  # Msun/year
            partx = p0['Coordinates'][:,0]      # Mpc/h
            party = p0['Coordinates'][:,1]
            partz = p0['Coordinates'][:,2]
            pvx = p0['Velocity'][:,0]           # km/s
            pvy = p0['Velocity'][:,1]
            pvz = p0['Velocity'][:,2] 
        else:
            partID  = np.append(partID,p0['ParticleIDs'][:]) 
            BH_Mass = np.append(BH_Mass,p0['BH_Mass'][:])
            BH_Mdot = np.append(BH_Mdot,p0['BH_Mdot'][:]*10**7*ast.s_in_year/ast.m_in_pc)
            partx       = np.append(partx,p0['Coordinates'][:,0])
            party       = np.append(party,p0['Coordinates'][:,1])
            partz       = np.append(partz,p0['Coordinates'][:,2])
            pvx         = np.append(pvx,p0['Velocity'][:,0])
            pvy         = np.append(pvy,p0['Velocity'][:,1])
            pvz         = np.append(pvz,p0['Velocity'][:,2])
            
    if verbose:
        print('BH: seed={:.2e}; min={:.2e}, max={:.2e}'.format(BH_seed_mass,
                                                               min(BH_Mass)*10**10,
                                                               max(BH_Mass)*10**10))

    # Get subgrid information into a pandas dataset to facilitate merging options
    #here: This operation changes groupnum and subgroupnum into floats, but doesn't seem to matter
    #      tried dtype=[np.uint32,np.int32,np.int32,np.float64,np.float64,np.float64,np.float32,np.float32,np.float32]
    data = np.c_[partID,BH_Mass,BH_Mdot,partx,party,partz,pvx,pvy,pvz]
    partID,BH_Mass,BH_Mdot,partx,party,partz,pvx,pvy,pvz=[[] for i in range(9)] 
    df_pbh = pd.DataFrame(data=data,columns=['partID','BH_Mass','BH_Mdot',
                                             'partx','party','partz',
                                             'pvx','pvy','pvz'])
    data= []

    # Get the halo the BH particles belong to ----------------------------------
    files, allfiles = get_particle_files(snap,sim,env,subfind=True)
    if (not allfiles):
        print('\n WARNING (b.get_subBH): no adequate particle files found, {}, {}'.
              format(snap,env))
        return None
    if Testing: files = [files[0]]
    if verbose: print('\n Subfind particles: {} \n'.format(files[0]))

    # Loop over the particle files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff,inptype)
        p0 = f[inptype]
    
        # Read particle information
        if (iff == 0):
            partID = p0['ParticleIDs'][:] 
            groupnum = p0['GroupNumber'][:] # FoF group number particle is in
            # Negative values: particles within r200 but not part of the halo
            subnum = p0['SubGroupNumber'][:]
        else:
            partID    = np.append(partID,p0['ParticleIDs'][:]) 
            groupnum  = np.append(groupnum,p0['GroupNumber'][:])
            subnum = np.append(groupnum,p0['SubGroupNumber'][:])
    
    if verbose:
        print('GroupNum: min={:d}, max={:d}'.format(min(groupnum),max(groupnum)))

    # If all groupnum are less than 0, take abs(), returning the Subhalo/GroupNumber values
    allgneg = False
    ind = np.where(groupnum<0)
    if(np.shape(ind)[1] == len(groupnum)):
        allgneg = True
        groupnum = abs(groupnum)-1

    # Get particle information into a pandas dataset to facilitate merging options
    data = np.c_[partID,groupnum,subnum]
    partID,groupnum,subnum = [[] for i in range(3)]
    df_psub = pd.DataFrame(data=data,columns=['partID','groupnum','subgroupnum'])

    # Join the particle information
    df_part = pd.merge(df_psub, df_pbh, on=['partID'])
    df_part.sort_values(by=['groupnum'], inplace=True)
    df_part.reset_index(inplace=True, drop=True)  # Reset index from 0
    del df_psub, df_pbh
    
    if (new_subnum):  # Assign subgrid BH particles to subhaloes
        inilen = len(df_part)
        inr = 0
        
        # Loop over groupnum and the substructure within
        for index, row in df_part.iterrows():
            ind = np.where(gn == row['groupnum'])
            if (np.shape(ind)[1]>=1):
                subnum = sgn[ind]
                subx = cop_x[ind]
                suby = cop_y[ind]
                subz = cop_z[ind]

                mindr = 2*boxside
                minsub  = -99
                for ii,sub in enumerate(subnum):
                    # Measure dr and find the min
                    dr1 = get_r(subx[ii],suby[ii],subz[ii],
                               row['partx'],row['party'],row['partz'],box=boxside)
                    if (dr1 <= mindr):
                        mindr = dr1
                        minsub = sub
                        
                # Assign the subhalo number to the BH particle
                df_part['subgroupnum'].at[index] = minsub

                # Fill vector with relative distances
                if inr == 0:
                    dr = np.array([mindr])
                    inr += 1
                else:
                    dr = np.append(dr,mindr)
            else:
                # Remove particles associated to haloes without stellar mass
                df_part = df_part.drop(labels=index, axis=0)

        ndrop = inilen-len(df_part)
        if ndrop>0:
            df_part.reset_index(inplace=True, drop=True)  # Reset index from 0
            if verbose:print(' {} rows have been deleted from df_part'.format(ndrop))

    # Data frame from subhaloes for BHs, to be able to merge with particle info.
    data = np.c_[gn,sgn,cop_x,cop_y,cop_z,shv_x,shv_y,shv_z]
    gn,sgn,cop_x,cop_y,cop_z,shv_x,shv_y,shv_z=[[] for i in range(8)]
    df_s4bh = pd.DataFrame(data=data,columns=['groupnum','subgroupnum','cop_x','cop_y','cop_z',
                                             'shv_x','shv_y','shv_z'])
    data=[]
    
    # Join information from subhaloes and particles
    df_all = pd.merge(df_s4bh, df_part, on=['groupnum','subgroupnum'])
    if verbose: print(df_all[['partID','groupnum','subgroupnum','cop_x','partx']],
                      df_all.columns.tolist())
    del df_s4bh,df_part
    
    if (not new_subnum): 
        # Calculate relative distances and velocities
        dr = get_r(df_all['cop_x'].to_numpy(),df_all['cop_y'].to_numpy(),df_all['cop_z'].to_numpy(),
                   df_all['partx'].to_numpy(),df_all['party'].to_numpy(),df_all['partz'].to_numpy(),
                   box=boxside)

    # Radial velocity
    dvr = get_vr(df_all['cop_x'].to_numpy(),df_all['cop_y'].to_numpy(),df_all['cop_z'].to_numpy(),
                 df_all['partx'].to_numpy(),df_all['party'].to_numpy(),df_all['partz'].to_numpy(),
                 df_all['shv_x'].to_numpy(),df_all['shv_y'].to_numpy(),df_all['shv_z'].to_numpy(),
                 df_all['pvx'].to_numpy(),df_all['pvy'].to_numpy(),df_all['pvz'].to_numpy(),
                 box=boxside)

    # Tangential velocity
    dvphi = get_vphi(df_all['cop_x'].to_numpy(),df_all['cop_y'].to_numpy(),df_all['cop_z'].to_numpy(),
                 df_all['partx'].to_numpy(),df_all['party'].to_numpy(),df_all['partz'].to_numpy(),
                 df_all['shv_x'].to_numpy(),df_all['shv_y'].to_numpy(),df_all['shv_z'].to_numpy(),
                 df_all['pvx'].to_numpy(),df_all['pvy'].to_numpy(),df_all['pvz'].to_numpy(),
                 box=boxside)
    
    # Line of site (z-axis) velocity
    dvlos = get_vlos(df_all['cop_x'].to_numpy(),df_all['cop_y'].to_numpy(),df_all['cop_z'].to_numpy(),
                 df_all['partx'].to_numpy(),df_all['party'].to_numpy(),df_all['partz'].to_numpy(),
                 df_all['shv_x'].to_numpy(),df_all['shv_y'].to_numpy(),df_all['shv_z'].to_numpy(),
                 df_all['pvx'].to_numpy(),df_all['pvy'].to_numpy(),df_all['pvz'].to_numpy(),
                 box=boxside)


    # Add properties of particles in the same position--------------------------
    groups = df_all.groupby(['groupnum','subgroupnum','partx','party','partz'], as_index=False)
    if(groups.ngroups > len(df_all.index)):
        print('\n WARNING (b.get_subBH): problem grouping particle information, {}, {}'.
          format(snap,sim))
        return None
    
    # Number of particles with the same position
    df_nboson = groups.size().reset_index(name='nboson')
    if verbose: print('Max. BH bosons = {}'.format(df_nboson['nboson'].max()))

    # Add the properties for particles in the same position
    df_addM_BH = groups.BH_Mass.sum() # 1e10 Msun/h
    df_addMdot = groups.BH_Mdot.sum() # Msun/year
    df_addM1 = pd.merge(df_addM_BH, df_addMdot, on=['groupnum','subgroupnum','partx','party','partz'])
    del df_addM_BH, df_addMdot

    df_addM = pd.merge(df_addM1, df_nboson, on=['groupnum','subgroupnum','partx','party','partz'])
    del df_addM1, df_nboson

    # Generate a bhnum column from 1 to the last group of BH particles in the same position
    df_addM.index += 1
    df_addM.index.names = ['bhnum'] 
    df_addM.reset_index(inplace=True)

    # Add bhnum to df_all with BH info
    df_bh = pd.merge(df_addM[['bhnum','groupnum','subgroupnum','partx','party','partz']], df_all,
                     on=['groupnum','subgroupnum','partx','party','partz'])
    if verbose: print(df_bh[['bhnum','groupnum','subgroupnum','partx','party','partz']],
                      df_bh.columns.tolist())
    del df_all
    
    # Add information to output file------------------------------------------------------    
    hf = h5py.File(outfile, 'a')
    hfdat = hf['data']

    # BH/  partID, groupnum, subgroupnum, MBH, Mdot, bhnum, pos, vel, dr, dvr, dvphi, dvlos
    bhdat = hfdat.create_group('BH')

    noms = ['partID','groupnum','subgroupnum',
            'BH_Mass','BH_Mdot','bhnum']
    desc = ['Particle ID','FoF group number','Subhalo index, cen:0',
            '1e10Msun/h', 'Msun/year','Group of BH particles in the same position']
    for ip, iprop in enumerate(noms):
        nom = noms[ip]
        prop = df_bh[[nom]].to_numpy()
        bhdat.create_dataset(nom,data=np.squeeze(prop)); prop = []
        bhdat[nom].dims[0].label = desc[ip]
    
    prop = df_bh[['partx', 'party', 'partz']].to_numpy()
    bhdat.create_dataset('pos',data=prop); prop = []
    bhdat['pos'].dims[0].label = 'x,y,z (Mpc/h)'

    prop = df_bh[['pvx', 'pvy', 'pvz']].to_numpy()
    bhdat.create_dataset('vel',data=prop); prop = []
    bhdat['vel'].dims[0].label = 'vx,vy,vz (km/s)'

    prop = dr; nom = 'dr'
    desc = 'Relative distance to FOF COP for centrals, and to central COP for sat. gal. (cMpc/h)'
    bhdat.create_dataset(nom,data=prop); prop = []
    bhdat[nom].dims[0].label = desc

    prop = dvr; nom = 'dvr'
    desc = 'Relative radial velocity to central gal. (km/s)'
    bhdat.create_dataset(nom,data=prop); prop = []
    bhdat[nom].dims[0].label = desc
    
    prop = dvphi; nom = 'dvphi'
    desc = 'Relative tangential (phi) velocity to central gal. (km/s)'
    bhdat.create_dataset(nom,data=prop); prop = []
    bhdat[nom].dims[0].label = desc
    
    prop = dvlos; nom = 'dvlos'
    desc = 'Relative velocity to central gal. on the line of sight, assumed z-axis (km/s)'
    bhdat.create_dataset(nom,data=prop); prop = []
    bhdat[nom].dims[0].label = desc

    
    #BHadd/  bhnum nboson, added_MBH, added_Mdot
    addat = hfdat.create_group('BHadd')

    noms = ['bhnum','nboson', 'BH_Mass', 'BH_Mdot']
    outn = ['bhnum','nboson', 'added_MBH', 'added_MdotBH']
    desc = ['Group of BH particles in the same position',
            'Number of BH particles in the same position (equal groupnum,subgroupnum,x,y,z)',
            'Added masses of BH particles in the same position (10**10 Msun/h)',
            'Added acretion rates for BH particles in the same position (Msun/year)']
    for ip, iprop in enumerate(noms):
        nom = noms[ip]
        prop = df_addM[[nom]].to_numpy()
        addat.create_dataset(outn[ip],data=np.squeeze(prop)); prop = []
        addat[outn[ip]].dims[0].label = desc[ip]
    
    hf.close()
    
    # Retrurn name of file with output
    return outfile


def map_subBH(snap,sim,env,nhmr=2.,cop=True,addp=False,
             dirz=None,outdir=None,Testing=True,verbose=False):
    '''
    Map subgrid BH properties into the half mass radius (HMR) of central subhaloes.
    Within 'arilega' there is not enough information to map on satellite subhaloes.

    Parameters
    -----------
    snap : int
        Snapshot number
    sim : string
        Name of the simulation
    env : string
        ari, arilega or cosma to use the adecuate paths
    nhmr : float
        Enclosure radius = nhmr*HalfMassRadius(DM)
    cop  : boolean
        True to use the CentreOfPotential, False for CentreOfMass 
    addp : boolean
        True to use the file with added properties for particles in the same position
    dirz : string
        Alternative path to table with z and snapshot.
    outdir : string
        Path to output file
    Testing : boolean
        Calculations on part or all the simulation
    verbose : boolean
        To output extra information

    Returns
    -------
    outfile : string
        Name of the output file

    Examples
    ---------
    >>> import bahamas as b
    >>> sim = 'HIRES/AGN_TUNED_nu0_L050N256_WMAP9'
    >>> b.map_subBH(31,sim,'arilega')
    '''

    # Stop for environments different to arilega
    if (env != 'arilega'):
        print('STOP: Function bahamas.map_subBH developed for env=arilega.')
        return None

    # Black hole particles to be read, 5:BH
    itype = 5 
    inptype = 'PartType'+str(itype)

    # Output file
    outfile, file_exists = get_subBH_file(outdir,sim,snap,nhmr=nhmr,cop=cop)
    
    # Get subgrid BH particle information
    partfile, file_exists = get_subBH_file(outdir,sim,snap,part=True,addp=addp)
    if not file_exists:
        partfile = get_subBH(snap,sim,env,addp=addp,
                             dirz=dirz,outdir=outdir,Testing=Testing,verbose=verbose)
        if (not partfile):
            print('STOP (b.map_subBH): no adequate particle files found, {}, {}'.
                  format(snap,env))
            return None

    f = h5py.File(partfile, 'r') #; print(ff,inptype)

    # Read the data
    groupnum= f['data/groupnum'][:,0] # FoF group number
    partx   = f['data/pos'][:,0]      # x (Mpc/h)
    party   = f['data/pos'][:,1]      # y (Mpc/h)
    partz   = f['data/pos'][:,2]      # z (Mpc/h)    
    BH_Mass = f['data/BH_Mass'][:,0]  # 1e10 Msun/h
    BH_Mdot = f['data/BH_Mdot'][:,0]  # Msun/year

    data = np.vstack([groupnum,partx,party,partz,BH_Mass,BH_Mdot]).T
    df_part = pd.DataFrame(data=data,columns=['groupnum','partx','party','partz',
                                              'BH_Mass','BH_Mdot'])
    groupnum,partx,party,partz,BH_Mass,BH_Mdot=[[] for i in range(6)]

    # Get halo information from FOF&Subfind files-------------------------------
    files, allfiles = get_subfind_files(snap,sim,env)
    if (not allfiles):
        print('STOP (b.map_subBH): no adequate Subfind files found, {}, {}'.
              format(snap,env))
        return None
    if Testing: files = [files[0],files[1]]
    if verbose: print('\n Haloes: {} \n'.format(files[0]))

    # Prop index
    stype = ptypes.index('star')
    dmtype = ptypes.index('DM')

    # Loop over the FOF&Subfind files
    for iff, ff in enumerate(files):
        f = h5py.File(ff, 'r') #; print(ff)
        sh = f['Subhalo']

        # Read halo information
        if (iff == 0):
            groupnum  = sh['GroupNumber'][:]      #FOF GroupNumber
            ms30  = sh['Mass_030kpc'][:,stype]    #1e10Msun/h
            HMRdm = sh['HalfMassRad'][:,dmtype]   #cMpc/h
            if cop:
                cop_x = sh['CentreOfPotential'][:,0]  #cMpc/h
                cop_y = sh['CentreOfPotential'][:,1]  #cMpc/h
                cop_z = sh['CentreOfPotential'][:,2]  #cMpc/h
            else:
                cop_x = sh['CentreOfMass'][:,0]  #cMpc/h
                cop_y = sh['CentreOfMass'][:,1]  #cMpc/h
                cop_z = sh['CentreOfMass'][:,2]  #cMpc/h
        else:
            groupnum  = np.append(groupnum,sh['GroupNumber'][:])
            ms30  = np.append(ms30,sh['Mass_030kpc'][:,stype])
            HMRdm = np.append(HMRdm,sh['HalfMassRad'][:,dmtype])
            if cop:
                cop_x = np.append(cop_x,sh['CentreOfPotential'][:,0])
                cop_y = np.append(cop_y,sh['CentreOfPotential'][:,1])
                cop_z = np.append(cop_z,sh['CentreOfPotential'][:,2])
            else:
                cop_x = np.append(cop_x,sh['CentreOfMass'][:,0])
                cop_y = np.append(cop_y,sh['CentreOfMass'][:,1])
                cop_z = np.append(cop_z,sh['CentreOfMass'][:,2])

    if verbose:
        print('{:d} galaxies read; GroupNum: min={:d}, max={:d}'.format(len(groupnum),
                                                                        min(groupnum),max(groupnum)))

    # Get indexes for centrals---------------------------------------
    cind = get_cenids(snap,sim,env)
    if (max(cind) > len(groupnum) and Testing):
        ind = np.where(cind < len(groupnum)-1) 
        if(np.shape(ind)[1]<1):
            print('STOP (b.map_subBH): no centrals in Subfind file.')
            return None
        cind = cind[ind]
    elif (max(cind) > len(groupnum) and not Testing):
        print('STOP (b.map_subBH): problem with centrals indexes.')
        return None
    
    if verbose:
        print('{:d} central galaxies'.format(len(cind)))
        print('Min. HMR = {:.3f} Mpc/h; Max. = {:.3f} Mpc/h'.format(
            min(HMRdm[cind]),max(HMRdm[cind])))
        print('Min. n*HMR = {:.3f} Mpc/h; Max. = {:.3f} Mpc/h'.format(
            min(nhmr*HMRdm[cind]),max(nhmr*HMRdm[cind])))

    # Mapping for central galaxies with stellar mass-----------------------------------------------
    data = np.vstack([groupnum[cind],ms30[cind],HMRdm[cind],cop_x[cind],cop_y[cind],cop_z[cind]]).T
    df_sh = pd.DataFrame(data=data,columns=['groupnum','ms30','HMRdm','cop_x','cop_y','cop_z'])
    data,groupnum,ms30,HMRdm,cop_x,cop_y,cop_z=[[] for i in range(7)]

    df_sh = df_sh.loc[df_sh.ms30 > 0.] # With stellar mass
    if df_sh.empty:
        print('STOP (b.map_subBH): no centrals with stellar mass.')
        return None

    # Join the particle and FoF information----------
    merge = pd.merge(df_part, df_sh, on=['groupnum'])
    del df_part

    # Get the boxsize
    omega0, omegab, lambda0, h0, boxsize = get_cosmology(sim,env)
    lbox2 = boxsize/2.

    # Position of particles relative to the center of the group
    merge['partx'] = merge.partx - merge.cop_x
    merge['party'] = merge.party - merge.cop_y
    merge['partz'] = merge.partz - merge.cop_z

    # Correct for periodic boundary conditions (for gal. in groups)
    merge.partx.loc[merge.partx < -lbox2] = merge.partx.loc[merge.partx < -lbox2] + boxsize
    merge.party.loc[merge.party < -lbox2] = merge.party.loc[merge.party < -lbox2] + boxsize
    merge.partz.loc[merge.partz < -lbox2] = merge.partz.loc[merge.partz < -lbox2] + boxsize

    merge.partx.loc[merge.partx >= lbox2] = merge.partx.loc[merge.partx >= lbox2] - boxsize
    merge.party.loc[merge.party >= lbox2] = merge.party.loc[merge.party >= lbox2] - boxsize
    merge.partz.loc[merge.partz >= lbox2] = merge.partz.loc[merge.partz >= lbox2] - boxsize

    # Distances to selected particles
    merge['distance'] = (merge.partx**2 +     
                         merge.party**2 +
                         merge.partz**2) ** 0.5
    if verbose: print('Min. distance to centre = {:.3f} Mpc/h; Max. = {:.3f} Mpc/h'.format(
            merge['distance'].min(),merge['distance'].max()))

    # Particles enclosed in radius n*HMR(DM)
    radius = nhmr*merge.HMRdm
    merge['inside_HMRdm'] = merge.distance <= radius;
    merge = merge.loc[merge.inside_HMRdm == True]

    if merge.empty:
        print('STOP (b.map_subBH): no particles within DM HMR.')
        return None
    groups = merge.groupby(['groupnum'], as_index=False)

    # BH mass and mdot of particles within that radius
    #            Output: pandas df w groupnum and prop
    massinHMRdm = groups.BH_Mass.sum() # 1e10 Msun/h
    mdotinHMRdm = groups.BH_Mdot.sum() # Msun/year

    minHMRdm = pd.merge(massinHMRdm, mdotinHMRdm, on=['groupnum'])
    del massinHMRdm, mdotinHMRdm

    final = pd.merge(minHMRdm, df_sh, on=['groupnum'])
    del minHMRdm, df_sh
    if verbose: print(final)
    
    # Write properties to output file
    hf = h5py.File(outfile, 'w') # Generate the file
    if verbose: print('\n Output: {} \n'.format(outfile))

    # Output header
    headnom = 'header'
    head = hf.create_dataset(headnom,(100,))
    head.attrs[u'sim']           = sim
    head.attrs[u'snapshot']      = snap
    head.attrs[u'redshift']      = get_z(snap,sim,env,dirz=dirz)
    head.attrs[u'omega0']        = omega0
    head.attrs[u'omegab']        = omegab
    head.attrs[u'lambda0']       = lambda0        
    head.attrs[u'h0']            = h0
    head.attrs[u'boxsize']       = boxsize
    head.attrs[u'units_boxsize'] = u'Mpc/h'

    # Output data with units
    hfdat = hf.create_group('data')
    
    prop = final[['cop_x', 'cop_y', 'cop_z']].to_numpy()
    hfdat.create_dataset('pos',data=prop); prop = []
    hfdat['pos'].dims[0].label = 'x,y,z (Mpc/h)'

    prop = final[['groupnum']].to_numpy()
    hfdat.create_dataset('groupnum',data=prop); prop = []
    hfdat['groupnum'].dims[0].label = 'FoF group number' 

    prop = final[['ms30']].to_numpy()
    hfdat.create_dataset('ms30',data=prop); prop = []
    hfdat['ms30'].dims[0].label = '1e10 Msun/h' 

    prop = final[['HMRdm']].to_numpy()
    hfdat.create_dataset('HMRdm',data=prop); prop = []
    hfdat['HMRdm'].dims[0].label = 'cMpc/h'
    
    prop = final[['BH_Mass']].to_numpy()
    hfdat.create_dataset('BH_Mass',data=prop); prop = []
    hfdat['BH_Mass'].dims[0].label = '1e10 Msun/h' 

    prop = final[['BH_Mdot']].to_numpy()

    hfdat.create_dataset('BH_Mdot',data=prop); prop = []
    hfdat['BH_Mdot'].dims[0].label = 'Msun/year'  
   
    hf.close()

    # Retrurn name of file with output
    return outfile


if __name__== "__main__":
    dirz = None ; outdir = None
    #snap = 31
    #zz = 0.
    snap = 27 #31
    zz = 0.75 #0.

    env = 'arilega'
    #env = 'cosmalega'
    #env = 'lap'

    
    if (env == 'cosmalega'):
        sim = None #'L400N1024/WMAP9/Sims/BAHAMAS'
        dirz = '/cosma7/data/dp004/dc-gonz3/BAHAMAS/'
        outdir = '/cosma7/data/dp004/dc-gonz3/Junk/'
    if (env == 'arilega'):
        #sim = 'HIRES/AGN_RECAL_nu0_L100N512_WMAP9'
        sim = 'AGN_TUNED_nu0_L400N1024_WMAP9'
        dirz = '/users/arivgonz/output/BAHAMAS/'
        outdir = '/users/arivgonz/output/Junk/'
    if (env == 'ari'):
        sim = 'L050N256/WMAP9/Sims/ws_324_23_mu_7_05_dT_8_35_n_75_BH_beta_1_68_msfof_1_93e11'
    if (env == 'lap'):
        sim = 'HIRES/AGN_RECAL_nu0_L100N512_WMAP9'
        #sim = 'AGN_TUNED_nu0_L400N1024_WMAP9'
        dirz = '/home/violeta/soil/BAHAMAS/'
        outdir = '/home/violeta/Downloads/'
        
    #print(get_particle_files(snap,sim,env,subfind=False))
    print(get_subnum(sim,snap,env,Testing=False))
    #print(old_get_subBH(snap,sim,env,addp=True,dirz=dirz,outdir=outdir,Testing=True,verbose=True))
    #print(get_subBH_file(outdir,sim,snap)) #,part=True))
    #print(get_subBH(snap,sim,env,dirz=dirz,outdir=outdir,Testing=True,verbose=True))
    #print(map_subBH(snap,sim,env,dirz=dirz,outdir=outdir,Testing=True,verbose=True))
    #print(get_mHMRmap_file(outdir,sim,snap))
    #print(map_mHMR(snap,sim,env,ptype='BH',nhmr=2.,cop=True,dirz=dirz,outdir=outdir,verbose=True))
    #print(get_m500_file(outdir,sim,snap))
    #print(map_m500(snap,sim,env,overwrite=True,dirz=dirz,outdir=outdir,Testing=True,verbose=True))
    #print(get_zminmaxs([0.,1.],dz=0.5))
    #print(get_simlabels(['AGN_TUNED_nu0_L100N256_WMAP9',
    #               'HIRES/AGN_RECAL_nu0_L100N512_WMAP9',
    #               'L400N1024/WMAP9/Sims/BAHAMAS']))
    #print(get_outdirs(env,dirz=dirz,outdir=outdir))
    #print(get_dirb(env),get_dirobs(env))
    #print(table_z_sn(sim,env,dirz=dirz))
    #print(get_z(27,sim,env,dirz=dirz))
    #print(get_z(-1,sim,env,dirz=dirz))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz)
    #print(boundary_correction(np.array([110.,2,-0.5,100.]),100.))
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz,zmax=[3.8])
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz,zmin=[2.8])
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #snap, zsnap = get_snap(3.2,sim,env,dirz=dirz,zmin=[2.8],zmax=[3.8])
    #print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    #print(get_allparticle_files(snap,sim,env))
    #print(get_cenids(snap,sim,env))
    #print(get_subfind_prop(snap,sim,env,'Subhalo/Mass_030kpc',proptype='star',Testing=True))
    #print('-------'); print(get_subfind_prop(snap,sim,env,'FOF/Group_M_Crit200',Testing=True))
    #print(resolution(sim,env,dirz=dirz))
    #print('log10(SFR (Msun/Gyr)) = {:2f}'.format(np.log10(get_min_sfr(sim,env,dirz=dirz))+9))
    #print(get_nh(zz,'FOF/Group_M_Mean200',sim,env,dirz=dirz,outdir=outdir))
    #print(get_propfunc(zz,['FOF/Group_M_Mean200','FOF/m2'],
    #                   'mass',sim,env,ptype='DM',dirz=dirz,outdir=outdir))

    #infile = '/hpcdata0/simulations/BAHAMAS/AGN_TUNED_nu0_L100N256_WMAP9/Data/Snapshots/snapshot_026/snap_026.27.hdf5'
    #infile = '/hpcdata0/simulations/BAHAMAS/AGN_TUNED_nu0_L100N256_WMAP9/Data/EagleSubGroups_5r200/groups_026/eagle_subfind_tab_026.0.hdf5'
    #print(print_h5attributes(infile,'Constants'))
