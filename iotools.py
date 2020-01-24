import sys
import os.path

def check_file(infile):
	if (not os.path.isfile(infile)):
		print('STOP: no input file {}'.format(infile))
		sys.exit()
	return
