import os.path
import sys

g = 6.673*10**(-11)  #Nm**2kg-2
k = 1.381*10**(-23)  #J/K
c = 2.998*10**8      #m/s
h = 6.626*10**(-34)  #Js
u = 1.661*10**(-27)  #kg
mp = 1.673*10**(-27) #kg
sb = 5.67*10**(-8)   #Wm^-2K^-4
msun = 1.99*10**(30) #kg
rsun = 6.96*10**8   #m
lsun = 3.83*10**(26) #W
teff = 5780          #K
tsun = 15.6*10**6    #K
rhosun = 1.48*10**5  #km m**-3
ev =1.602*10**(-19)  #J
au = 1.496*10**(11)  #m

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
    Return True if directory has been created
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


if __name__== "__main__":
    infile = 'blu'
    outdir = 'remove_dir/blu'

    print(nf('4.0'))
    print('Check file {}: {}'.format(infile,check_file(infile)))
    print('Create dir {}: {}'.format(outdir,create_dir(outdir)))
    print(is_sorted([1,3,6]))
    print(is_sorted([1,9,6]))

    print(stop_if_no_file(infile))
