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

                
def check_file(infile):
	if (not os.path.isfile(infile)):
		print('STOP: no input file {}'.format(infile))
		sys.exit()
	return

