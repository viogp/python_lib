#!/usr/bin/env python

import os
from datetime import datetime
now = datetime.now()

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    
# Job description
jobnames = ['test']
time = '01:00:00'  #Format: d-hh:mm:ss
nodes = 1
ntasks = 1
cputask = 1
partition = 'compute'  # test/compute 

# Get current working directory
path2program = os.getcwd()+'/'+'bahamas.py'  

# Logs output directory
log_dir = '/users/arivgonz/output/logs/'
mkdir_p(log_dir)

# For checking if the module is loaded
script_content1 = 'module_to_check="apps/anaconda3/2023.03/bin"'
script_content2 = 'if [[ $module_list_output != *"$module_to_check"* ]]; then'

# Submission script
for jobname in jobnames:
    #script = log_dir+jobname+now.strftime('_%Y%m%d_%H%M')+'.sh'
    script = log_dir+jobname+'.sh'
    
    with open(script,'w') as fh:
        fh.write("#!/bin/bash \n")
        fh.write("\n")
        fh.write("#SBATCH --job-name=job.{}.%J \n".format(jobname))
        fh.write("#SBATCH --output={}out.{}.%J \n".format(log_dir,jobname))
        fh.write("#SBATCH --error={}err.{}.%J \n".format(log_dir,jobname))
        fh.write("#SBATCH --time={} \n".format(time))
        fh.write("#SBATCH --nodes={} \n".format(str(nodes)))
        fh.write("#SBATCH --ntasks={} \n".format(str(ntasks)))
        fh.write("#SBATCH --cpus-per-task={} \n".format(str(cputask)))
        fh.write("#SBATCH --partition={} \n".format(partition))
        #fh.write("#SBATCH --mail-type=END,FAIL \n")         #Not working
        #fh.write("#SBATCH --mail-user=violetagp@protonmail.com \n")
        fh.write("\n")
        fh.write("flight env activate gridware \n")
        fh.write("\n")
        fh.write("{} \n".format(script_content1))
        fh.write("module_list_output=$(module list 2>&1) \n")
        fh.write("{} \n".format(script_content2))
        fh.write("  module load $module_to_check \n")
        fh.write("fi \n")
        fh.write("\n")    
        fh.write("python3 {} \n".format(path2program))
    
    print("Run {}".format(script))
    os.system("sbatch {}".format(script))
