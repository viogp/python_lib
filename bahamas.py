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

def table_z_sn(sim,outdir):
    output = outdir+sim+'/'+tblz
    path = dirbahamas+sim+'/Data/EagleSubGroups_5r200/'
    dirs = glob.glob(path+'groups_*')

    # Initialize arrays for z and sn
    zzs = np.zeros(shape=len(dirs)) ; zzs.fill(-999.)
    sns = np.zeros(shape=len(dirs), dtype=int)

    for ii, dir in enumerate(dirs):
        snap = dir.split('groups_0')[1]
        sns[ii] = float(snap)

        infile = dir+'/eagle_subfind_tab_0'+snap+'.0.hdf5'
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
    with open(output, 'w') as f:
        f.write('# Redshift Snapshot \n')
        for ii in range(len(zz)):
            tofile = '{:.2f} {:d}'.format(zz[ii],sn[ii])
            f.write("%s\n" % tofile)

    return output


def get_z(snap,sim,outdir):
    """
    Get the redshift for a snapshot number and
    a simulation name

    Parameters
    -----------
    snap : int
        Snapshot number
    sim : string
        Name of the Bahamas simulation
    outdir : string
        Path to table with z and snapshot

    Returns
    -----
    snap : string
        Snapshot value

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_z(26,'AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/')
    >>>
    """

    # Table with increasing redshifts and corresponding snapshots
    tablez = outdir+sim+'/'+tblz
    if (not os.path.isfile(tablez)):
        # Generate the table if it doesn't exist
        tablez = table_z_sn(sim,outdir)

    # Read the table:
    zzs, x = np.loadtxt(tablez, unpack=True)
    sns = np.asarray(x,dtype=int)
    
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


def get_snap(zz,zmin,zmax,sim,outdir):
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
    outdir : string
        Path to table with z and snapshot

    Returns
    -----
    snap : int
        Snapshot value
    z_snap : float
        Redshift for that snapshot

    Examples
    ---------
    >>> import bahamas as b
    >>> b.get_snap(3.2,2.8,3.6,'AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/')
    >>>
    """

    # Table with increasing redshifts and corresponding snapshots
    tablez = outdir+sim+'/'+tblz
    if (not os.path.isfile(tablez)):
        # Generate the table if it doesn't exist
        tablez = table_z_sn(sim,outdir)

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


def cenids(snap,sim):
    """
    Get the list of indexes for central galaxies

    Parameters
    -----------
    snap : integer
        Snapshot 
    sim : string
        Simulation name

    Returns
    -----
    cenids : numpy array int
        Global indexes for central galaxies

    Examples
    ---------
    >>> import bahamas as b
    >>> b.cenids(31,AGN_TUNED_nu0_L050N256_WMAP9)
    """

    path = dirbahamas+sim+'/Data/EagleSubGroups_5r200/groups_0'+str(snap)
    root = path+'/eagle_subfind_tab_0'+str(snap)+'.'
    files = glob.glob(root+'*.hdf5')
    lenf = len(files) ; files = []
    if (lenf<1):
        print('STOP(~/python_lib/bahamas): Make sure you can see the path {}'.format(path))
        sys.exit()

    for ii in range(lenf):
        ff = root+str(ii)+'.hdf5'
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
    print(get_z(-1,'AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/'))
    print(get_z(26,'AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/'))

    snap,zsnap = get_snap(3.2,2.8,3.8,'AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/')
    print('target z={} -> snap={}, z_snap={}'.format(3.2,snap,zsnap))
    print(get_snap(-100.,-200.,-5.,'AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/'))
    print(get_snap(0.28,0.26,0.3,'AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/'))

    print(get_zminmaxs([0.]))
    print(get_zminmaxs([0.,1.],dz=0.5))
    print(get_zminmaxs([0.,1.]))
    print(get_zminmaxs([0.,1.,0.5]))

