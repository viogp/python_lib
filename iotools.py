import sys
import os.path

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
