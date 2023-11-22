import os

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    
# Job description
nom = 'nd_bahamas'
time = '00:15:00'  #Format: d-hh:mm:ss
nodes = 1
ntasks = 1
cputask = 1
partition = 'test'  # test/compute 

# Get current working directory
job = os.getcwd()+'/'+'bahamas.py'  

# Logs output directory
log_dir = '/users/arivgonz/output/logs/'
mkdir_p(log_dir)

# For checking if the module is loaded
script_content1 = 'module_to_check="apps/anaconda3/2023.03/bin" \n'

# Submission script
script = os.path.join(log_dir,"%s.sh" % nom)
with open(script,'w') as fh:
    fh.write("#!/bin/bash \n")
    fh.write("\n")
    fh.write("#SBATCH --job-name=%s.job \n" % nom)
    fh.write("#SBATCH --output=%sout.%s \n" % (log_dir,nom))
    fh.write("#SBATCH --error=%serr.%s \n" % (log_dir,nom))
    fh.write("#SBATCH --time=%s \n" % time)
    fh.write("#SBATCH --nodes=%s \n" % str(nodes))
    fh.write("#SBATCH --ntasks=%s \n" % str(ntasks))
    fh.write("#SBATCH --cpus-per-task=%s \n" % str(cputask))
    fh.write("#SBATCH --partition=%s \n" % partition)
    fh.write("#SBATCH --mail-type=END,FAIL \n")
    fh.write("#SBATCH --mail-user=violetagp@protonmail.com \n")
    fh.write("\n")
    fh.write("flight env activate gridware \n")
    fh.write("\n")
    fh.write("%s \n" % script_content1)
    fh.write("module_list_output=$(module list 2>&1) \n")
    fh.write("if [[ $module_list_output !~ $module_to_check ]]; then \n")
    fh.write("  module load $module_to_check \n")
    fh.write("fi \n")
    fh.write("\n")    
    fh.write("python3 %s \n" % job)

print("Log dir.: %s" % log_dir)
os.system("sbatch %s" % script)
