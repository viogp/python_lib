import sys
import os.path

class nf(float):
        # Define a class that forces representation of float to look a certain way
        # This remove trailing zero so '1.0' becomes '1'
        def __repr__(self):
                str = '%.1f' % (self.__float__(),)
                if str[-1] == '0':
                        return '%.0f' % self.__float__()
                else:
                        return '%.1f' % self.__float__()

def stop_if_no_file(infile):
	if (not os.path.isfile(infile)):
		print('STOP: no input file {}'.format(infile))
		sys.exit()
	return 

def check_file(infile):
	file_fine = True
	if (not os.path.isfile(infile)):
		file_fine = False
	return file_fine

def is_sorted(a):
    for i in range(a.size-1):
         if a[i+1] < a[i] :
               return False
    return True
