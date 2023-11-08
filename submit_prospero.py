import os

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    

# Get current working directory
job_directory = "%s/.job" %os.getcwd()  

# Logs output directory
log_dir = '/tmp/users/arivgonz/logs/'
mkdir_p(log_dir)

# Job file
nom = 'bahamas'
job_file = os.path.join(job_directory,"%s.job" % nom)
print(job_file)
#    nom_data = os.path.join(data_dir, nom)
#
#    # Create nom directories
#    mkdir_p(nom_data)
#
#    with open(job_file) as fh:
#        fh.writelines("#!/bin/bash\n")
#        fh.writelines("#SBATCH --job-name=%s.job\n" % nom)
#        fh.writelines("#SBATCH --output=.out/%s.out\n" % nom)
#        fh.writelines("#SBATCH --error=.out/%s.err\n" % nom)
#        fh.writelines("#SBATCH --time=2-00:00\n")
#        fh.writelines("#SBATCH --mem=12000\n")
#        fh.writelines("#SBATCH --qos=normal\n")
#        fh.writelines("#SBATCH --mail-type=ALL\n")
#        fh.writelines("#SBATCH --mail-user=$USER@stanford.edu\n")
#        fh.writelines("Rscript $HOME/project/NomLips/run.R %s potato shiabato\n" %nom_data)
#
#    os.system("sbatch %s" %job_file)
