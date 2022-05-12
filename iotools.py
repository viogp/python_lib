import os.path
import sys
import h5py

class nf(float):
    '''
    Define a class that forces representation of float to look a certain way
    This removes trailing zero so '1.0' becomes '1'
    '''
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__() 
        else:
            return '%.1f' % self.__float__()
    

def stop_if_no_file(infile):
    '''
    Stop if the file does not exist
    '''
    if (not os.path.isfile(infile)):
        print('STOP: no input file {}'.format(infile)) 
        sys.exit()
    return


def check_file(infile,verbose=False):
    '''
    Return True if the file exists
    '''
    file_fine = True  
    if (not os.path.isfile(infile)):
        file_fine = False
        if verbose:
            print('WARNING (io.check_file): file not found {}'.format(infile))

    return file_fine


def create_dir(outdir):
    '''
    Return True if directory already exists or it has been created
    '''
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            print('WARNING (iotools.create_dir): problem creating directory ',outdir)
            return False
    return True


def is_sorted(a):
    '''
    Return True if the array is sorted
    '''
    for i in range(len(a)-1): 
        if a[i+1] < a[i] : 
            return False
    return True


def print_h5attr(infile,inhead='Header'):
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
    >>> import iotools as io
    >>> infile = '/hpcdata0/simulations/BAHAMAS/AGN_TUNED_nu0_L100N256_WMAP9/Data/Snapshots/snapshot_026/snap_026.27.hdf5'
    >>> io.print_header5(infile)
    """

    filefine = check_file(infile) #print(filefine)
    if (not filefine):
        print('WARNING (iotools.printh5attr): Check that the file provided is correct')
        return ' '
    
    f = h5py.File(infile, 'r')
    header = f[inhead]
    for hitem in list(header.attrs.items()): 
        print(hitem)
    f.close()

    return ' '



if __name__== "__main__":
    infile = 'blu'
    outdir = 'remove_dir/blu'

    print(nf('4.0'))
    print('Check file {}: {}'.format(infile,check_file(infile)))
    print('Create dir {}: {}'.format(outdir,create_dir(outdir)))
    print(is_sorted([1,3,6]))
    print(is_sorted([1,9,6]))

    print(stop_if_no_file(infile))
