#! /usr/bin/env python
#########################UNFINISHED####################################
def params_file_cute(cosmo,infile=input_file,output=corr_data,random=None)

h0 = cosmo[0]
omegaM = cosmo[1]
omegaL = cosmo[2]

if (random is None):
    random = '/gpfs/data/Lightcone/CUTE/CUTE-1.2/CUTE/test/random.dat'

params = '''
# input-output files and parameters
data_filename= '''+input_file+'''
random_filename= '''+random+'''
input_format= 1
mask_filename= none 
z_dist_filename= none
#IMPORTANT: CUTE assumes that the values of z are monotonously
#increasing. The normalization of N(z) is irrelevant.
output_filename= '''+corr_data+''' 
num_lines= all

# estimation parameters
corr_type= angular
corr_estimator= LS
np_rand_fact= 8

# cosmological parameters
omega_M= '''+str(omegaM)+'''
omega_L= '''+str(omegaL)+'''
w= -1          

# binning
log_bin= 0
n_logint= 10
dim1_max= 10.
dim1_nbin= 64
dim2_max= 0.1
dim2_nbin= 32
dim3_min= 0.4
dim3_max= 0.7
dim3_nbin= 1

# pixels for radial correlation
radial_aperture= 1

# pm parameters
use_pm= 1
n_pix_sph= 2048
'''
f.write(params)
f.close()
print 'Written parameters to file', param_file
