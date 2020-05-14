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

zs0 = [100.,5.,4.,3.,2.,1.,0.5,0.1,0.]
epsz = 0.0001

def get_zticks(lowestz,xmin,xmax):
    """
    Get the labels and positions of a secundary x-axis for redshifts, 
    given the age of the Universe as main x-axis

    Parameters
    -----------
    lowestz : float
        Lowest z read
    xmin : float
        Minimum age value of the main x-axis
    xmax : float
        Maximum age value of the main x-axis

    Returns
    -----
    zs : list of floats
        Values of z to be used as axis labels
    zticks : list of floats
        Position of the ticks, following the main x-axis

    Examples
    ---------
    >>> import plotbahamas as pb
    >>> import Cosmology as cosmo
    >>> cosmo.set_cosmology()
    >>> zs, zticks = pb.get_zticks(5.,0.12,12.)
    >>> print(zs)
    [100.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 5.0]
    """ 

    zticks0 = [cosmo.age_of_universe(x) for x in zs0]
    agelz = cosmo.age_of_universe(lowestz)

    if (lowestz < epsz):
        zs = zs0
        zticks = zticks0
    else:
        ind = next(x[0] for x in enumerate(zticks0) if x[1] > xmax) - 1
        if (ind == 0):
            if (lowestz == zs0[0]):
                zs = [zs0[0]]
                zticks = [zticks0[0]]
            else:
                zs = [zs0[0],lowestz]
                zticks = [zticks0[0],min(agelz,xmax)]
        elif (ind > 0):
            zs = zs0[0:ind]
            zticks = zticks0[0:ind]
            if (agelz < xmax):
                zs.append(lowestz)
                zticks.append(agelz)

    return zs,zticks

def wctime(sims,labels,env,dirplot=None,zrange=None):
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
    zrange : list of length 2 with floats
        zrange = [zmin, zmax], limits in redshift for the plot

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

    # Loop over all the simulations to be compared
    lowestz = 999
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

        # Find lowest redshift
        if (lowestz > float(redshift[-1])): lowestz = float(redshift[-1])

        # Set the cosmology for this simulation and calculate the ages
        omega0, omegab, lambda0, h0 = b.get_cosmology(sim,env)
        cosmo.set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0,h0=h0,
                            universe="Flat",include_radiation=False)
        age = [cosmo.age_of_universe(x) for x in redshift]

        # Plot wall clock time vs age
        ax.plot(age,wctime,color=cols[ii],label=labels[ii])

    # If specified, set z range using the last loaded cosmology
    if zrange is not None:
        if (len(zrange) == 2):            
            lowestz = zrange[0]
            xmax = cosmo.age_of_universe(lowestz)
            xmin = cosmo.age_of_universe(zrange[1])
            ax.set_xlim(xmin,xmax)
        else:
            print('WARNING (plotbahamas): zrange = [zmin,zmax]')

    # Top axis with redshift
    xmin, xmax = ax.get_xlim()
    zs, zticks = get_zticks(lowestz,xmin,xmax)
    axz = ax.twiny()
    axz.minorticks_off()
    axz.set_xticks(zticks)
    axz.set_xticklabels(['{:g}'.format(x) for x in zs])
    axz.set_xlim(xmin, xmax)
    axz.set_xlabel('z')

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

def cputime(sims,labels,env,dirplot=None,zrange=None):
    """
    Plot the CPU percentages time versus redshift and age
    a simulation name.

    Parameters
    -----------
    sim : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
    dirplot : string
        Path to plots
    zrange : list of length 2 with floats
        zrange = [zmin, zmax], limits in redshift for the plot

    Returns
    -----
    plotf : string
        Path to plot

    Examples
    ---------
    >>> import plotbahamas as pb
    >>> pb.cputime(['AGN_TUNED_nu0_L100N256_WMAP9'],['REF'],'ari')
    >>> pb.cputime(['L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_125_ntask128'],['tmax0_125_ntask128'],'cosma')
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
    ytit = '%CPU time'
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    #ax.set_xlim(xmin,xmax)
    ax.set_ylim(0.,100.) 
    plot_lines = [] ; plot_colors = []

    # Properties to be read each step
    props = ['treegrav','pmgrav','sph','eagle_total']
    ls = ['-',':','--','-.']

    # Loop over all the simulations to be compared
    lowestz = 999
    for ii, sim in enumerate(sims):
        # Simulation input
        if (env == 'ari'):
            dirb = b.dirbahamasari
            path = dirb+sim+'/Data/EagleSubGroups_5r200/'
        elif (env == 'cosma'):
            dirb = b.dirbahamascosma
        path = dirb+sim+'/data/'

        # Initialize lists
        redshift = [] 
        #times = [[] for i in range(len(props))]
        percentages = [[] for i in range(len(props))]

        # Read file with cpu information
        infof = path+'cpu.txt'
        start_time = time.time()
        icount = 0
        with open(infof, "r") as ff:
            for line in ff:
                if line.strip():
                    word = line.split()[0]
                    if (word == 'Step'):
                        s1 = line.split('z:')[1]
                        redshift.append(float(s1.split('CPUs')[0]))
                    elif (word in props):
                        ind = props.index(word)
                        #s1 = float(line.split()[1])
                        #times[ind].append(s1)

                        s1 = line.split()[2]
                        s2 = float(s1.split('%')[0])
                        percentages[ind].append(s2)

        # Find lowest redshift
        if (lowestz > float(redshift[-1])): lowestz = float(redshift[-1])
        if (lowestz > zs0[0]): lowestz = zs0[0]

        # Set the cosmology for this simulation and calculate the ages
        omega0, omegab, lambda0, h0 = b.get_cosmology(sim,env)
        cosmo.set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0,h0=h0,
                            universe="Flat",include_radiation=False)
        age = [cosmo.age_of_universe(x) for x in redshift]

        # Plot properties
        for jj, jprop in enumerate(props):
            y = percentages[jj]
            l1, = ax.plot(age,y,color=cols[ii],linestyle=ls[jj])

            if (ii==0):
                plot_lines.append(l1)
        plot_colors.append(l1)

    # If specified, set z range using the last loaded cosmology
    if zrange is not None:
        if (len(zrange) == 2):            
            lowestz = zrange[0]
            xmax = cosmo.age_of_universe(lowestz)
            xmin = cosmo.age_of_universe(zrange[1])
            ax.set_xlim(xmin,xmax)
        else:
            print('WARNING (plotbahamas): zrange = [zmin,zmax]')

    # Legend
    legend1 = ax.legend(plot_lines, props, loc=2)
    legend1.draw_frame(False)
    for h in legend1.legendHandles:
        h.set_color('k')
    plt.gca().add_artist(legend1)

    leg = ax.legend(plot_colors, labels, loc=1,
                    handlelength=0, handletextpad=0)
    leg.draw_frame(False)
    for ii,text in enumerate(leg.get_texts()):
        text.set_color(cols[ii])
    for item in leg.legendHandles:
        item.set_visible(False)

    # Top axis with redshift
    xmin, xmax = ax.get_xlim()
    zs, zticks = get_zticks(lowestz,xmin,xmax)
    axz = ax.twiny()
    axz.minorticks_off()
    axz.set_xticks(zticks)
    axz.set_xticklabels(['{:g}'.format(x) for x in zs])
    axz.set_xlim(xmin, xmax)
    axz.set_xlabel('z')

    # Path to plot
    if (dirplot == None):
        dirp = dirb+'plots/'+sim+'/'
    else:
        dirp = dirplot+sim+'/'

    if (not os.path.exists(dirp)):
        os.makedirs(dirp)
    plotf = dirp+'cputime.pdf'
    fig.savefig(plotf)

    return plotf
    

if __name__== "__main__":
    env = 'cosma'

    if (env == 'cosma'):
        sim1 = 'L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_01_ntask128'
        sim2 = 'L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_047619_ntask128'
        sim3 = 'L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_tmax0_125_ntask128' 

        sims = [sim1]#, sim2, sim3]
        labels = ['t$_{max}$=0.01']#,'t$_{max}$=0.05','t$_{max}$=0.125']
        print(wctime(sims,labels,'cosma'))
        #print(cputime(sims,labels,'cosma'))

