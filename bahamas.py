import os.path
import numpy as np
import h5py
import glob
#print('\n \n')

dirbahamas = '/hpcdata0/simulations/BAHAMAS/'

ptypes = ['gas','DM','bp1','bp2','star','BH']

dirobs = '/hpcdata0/Obs_Data/'

tblz = 'snap_z.txt'

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
    >>> b.get_z('26','AGN_TUNED_nu0_L100N256_WMAP9','/hpcdata3/arivgonz/bahamas/')
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
    a simulation name

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
    >>> b.get_snap(3.2,3.,3.4,'AGN_TUNED_nu0_L100N256_WMAP9')
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

if __name__== "__main__":
    print(get_z('-1','AGN_TUNED_nu0_L100N256_WMAP9'))
    print(get_z('32','AGN_TUNED_nu0_L100N256_WMAP9'))

    print(get_snap(100.,'AGN_TUNED_nu0_L100N256_WMAP9'))
    print(get_snap(-100.,'AGN_TUNED_nu0_L100N256_WMAP9'))
    print(get_snap(1.2,'AGN_TUNED_nu0_L100N256_WMAP9'))
