import sys,os.path
import numpy as np
import h5py
import glob
import pandas as pd
import time
from datetime import datetime
import Cosmology as cosmo
import bahamas as b
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)
#print('\n \n')

def wctime(sims,labels,env,dirplot=None):
    """
    Plot the Wall clock time versus redshift and age
    a simulation name.

    Parameters
    -----------
    sim : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
    dirplot : string
        Path to plots

    Returns
    -----
    plotf : string
        Path to plot

    Examples
    ---------
    >>> import plotbahamas as pb
    >>> pb.wctime(['AGN_TUNED_nu0_L100N256_WMAP9'],['REF'],'ari')
    >>> pb.wctime(['L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_125_ntask128'],['tmax0_125_ntask128'],'cosma')
    """ 

    # Check that the size of the arrays for the simulations and labels is the same
    lensims = len(sims)
    if (len(labels) != lensims):
        print('WARNING! Labels array does not match Sims array')
        # Generate labels
        labels = [x.split('/')[-1] for x in sims]

    # Set up plot variables
    fig = plt.figure()
    ax = plt.subplot()
    cols = get_distinct(lensims)
    xtit = 'Age (Gyr)'
    ytit = 'Wallclock time (s)'
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    #ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 

    # Loop over all the simulations to be compared
    for ii, sim in enumerate(sims):
        # Simulation input
        if (env == 'ari'):
            dirb = b.dirbahamasari
            path = dirb+sim+'/Data/EagleSubGroups_5r200/'
        elif (env == 'cosma'):
            dirb = b.dirbahamascosma
        path = dirb+sim+'/data/'

        # Read file with wall clock information
        infof = path+'info.txt'
        start_time = time.time()
        redshift = [] ; wctime = []
        with open(infof, "r") as ff:
            for line in ff:
                if line.strip():
                    char1 = line.strip()[0]
                    if (char1 == 'S'):
                        s1 = line.split('Redshift:')[1]
                        redshift.append(float(s1.split(',')[0]))
                    elif (char1 == 'D'):
                        s1 = line.split('[')[1]
                        s2 = s1.split(']')[0]
                        date = datetime.strptime(s2, '%Y-%m-%d %H:%M:%S.%f')
                        wctime.append(datetime.timestamp(date))
                        
        wc0 = wctime[0]
        wctime = [x - wc0 for x in wctime]

        # Set the cosmology for this simulation and calculate the ages
        omega0, omegab, lambda0, h0 = b.get_cosmology(sim,env)
        cosmo.set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0,h0=h0,
                            universe="Flat",include_radiation=False)
        age = [cosmo.age_of_universe(x) for x in redshift]

        # Plot wall clock time vs z
        ax.plot(age,wctime,color=cols[ii],label=labels[ii])

    ## Top axis with redshift   NOT WORKING
    #axz = ax.twiny()
    #axz.set_ylabel('z')
    #zs = [100.,20,10.,5.,4.,3.,2.,1.,0.5,0.1,0.]
    #zticks = [cosmo.age_of_universe(x) for x in zs]
    #axz.set_xticks(zticks)
    #axz.set_xticklabels(['{:g}'.format(x) for x in zs])

    # Add a legend
    leg = ax.legend(loc=0)
    leg.draw_frame(False)
    for ii,text in enumerate(leg.get_texts()):
        text.set_color(cols[ii])

    # Path to plot
    if (dirplot == None):
        dirp = dirb+'plots/'+sim+'/'
    else:
        dirp = dirplot+sim+'/'

    if (not os.path.exists(dirp)):
        os.makedirs(dirp)
    plotf = dirp+'wctime.pdf'
    fig.savefig(plotf)

    return plotf
    

if __name__== "__main__":
    env = 'cosma'

    if (env == 'cosma'):
        sim1 = 'L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_01_ntask128'
        sim2 = 'L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_01_ntask64'
        sim3 = 'L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_01_ntask32'
        sims = [sim1, sim2, sim3]
        labels = ['ntask128','ntask64','ntask32']
        print(wctime(sims,labels,'cosma'))

