import sys,os.path
import numpy as np
import h5py
import glob
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

def logformat(y,pos):
    '''
    To be used as follows:
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))
    (https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting/33213196)
    '''
    # Find the number of decimal places required
    decimalplaces = int(np.maximum(-np.log10(y),0))     # =0 for numbers >=1
    # Insert that number into a format string
    formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
    # Return the formatted tick label
    return formatstring.format(y)


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
    >>> import bahamasplot as bp
    >>> import Cosmology as cosmo
    >>> cosmo.set_cosmology()
    >>> zs, zticks = bp.get_zticks(5.,0.12,12.)
    >>> print(zs)
    [100.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 5.0]
    """ 

    if (lowestz > xmax):
        print('STOP (bahamasplot/get_zticks): lowest z is above {}'.format(xmax))
        sys.exit()

    zticks0 = [cosmo.age_of_universe(x) for x in zs0]
    agelz = cosmo.age_of_universe(lowestz)

    if (lowestz < epsz):
        return zs0,zticks0
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


def wctime(sims,env,labels=None,dirplot=None,zrange=None):
    """
    Plot the Wall clock time versus redshift and age
    a simulation name.

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
    labels : list of strings
        Array with the labels to be used
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
    >>> import bahamasplot as bp
    >>> bp.wctime(['AGN_TUNED_nu0_L100N256_WMAP9'],'ari')
    >>> bp.wctime(['L050N256/WMAP9/Sims/ex'],'cosma')
    """ 

    # Check labels
    labels = b.get_simlabels(sims,labels=labels)

    # Set up plot variables
    fig = plt.figure()
    ax = plt.subplot()
    cols = get_distinct(len(sims))
    xtit = 'Age (Gyr)'
    ytit = 'Wallclock time (h)'
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

    # Loop over all the simulations to be compared
    lowestz = 999
    for ii, sim in enumerate(sims):
        # Simulation input
        path = b.get_path2data(sim,env)

        # Read file with wall clock information
        infof = path+'info.txt'
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
        wctime = [x/3600. for x in wctime]  # From s to hours

        # Find lowest redshift
        if (lowestz > float(redshift[-1])): lowestz = float(redshift[-1])

        # Set the cosmology for this simulation and calculate the ages
        omega0, omegab, lambda0, h0 = b.get_cosmology(sim,env)
        cosmo.set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0,h0=h0,
                            universe="Flat",include_radiation=False)
        age = [cosmo.age_of_universe(x) for x in redshift]
        if (len(wctime) != len(age)):
            print('WARNING (wctime): problem with ages for sim:{}'.format(sim))
        else:
            # Plot wall clock time vs age
            ax.plot(age,wctime,c=cols[ii],label=labels[ii])

    # If specified, set z range using the last loaded cosmology
    if zrange is not None:
        if (len(zrange) == 2):            
            lowestz = zrange[0]
            xmax = cosmo.age_of_universe(lowestz)
            xmin = cosmo.age_of_universe(zrange[1])
            ax.set_xlim(xmin,xmax)
        else:
            print('WARNING (bahamasplot): zrange = [zmin,zmax]')

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
        dirp = b.get_dirb(env)+'plots/'+sim+'/'
    else:
        dirp = dirplot+sim+'/'

    if (not os.path.exists(dirp)):
        os.makedirs(dirp)
    plotf = dirp+'wctime.pdf'
    fig.savefig(plotf)

    return plotf

def cputime(sims,env,labels=None,dirplot=None,zrange=None):
    """
    Plot the CPU percentages time versus redshift and age
    a simulation name.

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
    labels : list of strings
        Array with the labels to be used
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
    >>> import bahamasplot as bp
    >>> bp.cputime(['AGN_TUNED_nu0_L100N256_WMAP9'],'ari')
    >>> bp.cputime(['L050N256/WMAP9/Sims/ex'],'cosma')
    """ 

    # Check labels
    labels = b.get_simlabels(sims,labels=labels)

    # Set up plot variables
    fig = plt.figure()
    ax = plt.subplot()
    cols = get_distinct(len(sims))
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
        path = b.get_path2data(sim,env)

        # Initialize lists
        redshift = [] 
        percentages = [[] for i in range(len(props))]

        # Read file with cpu information
        infof = path+'cpu.txt'
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
            if (len(y) != len(age)):
                print('WARNING (cputime): problem with the percentages for sim:{}'.format(sim))
            else:
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
            print('WARNING (bahamasplot): zrange = [zmin,zmax]')

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
        dirp = b.get_dirb(env)+'plots/'+sim+'/'
    else:
        dirp = dirplot+sim+'/'

    if (not os.path.exists(dirp)):
        os.makedirs(dirp)
    plotf = dirp+'cputime.pdf'
    fig.savefig(plotf)

    return plotf


def mf_sims(zz,massdef,sims,env,mmin=9.,mmax=16.,dm=0.1,labels=None,
            dirz=None,outdir=None,Testing=False):
    """
    Compare the mass function of different simulations at a given z

    Parameters
    -----------
    zz : float
        Redshift
    massdef : list of string
        Names of the mass definitions to be used
    sims : list of strings
        Array with the names of the simulation
    env : list of string
        ari(lega) or cosma(lega), to use the adecuate paths
    mmin : float
        Mimimum mass to be considered
    mmax : float
        Maximum mass to be considered
    dm : float
        Intervale step in halo mass   
    labels : list of strings
        Array with the labels to be used
    dirz : string
        Alternative path to table with z and snapshot
    outplot : string
        Path to output
    Testing : boolean
        True or False for testing with few subfiles

    Returns
    -----
    plotf : string
        Path to plot

    Examples
    ---------
    >>> import bahamasplot as bp
    >>> bp.mf_sims(27,['FOF/Group_M_Mean200'],['AGN_TUNED_nu0_L100N256_WMAP9'],['ari'])
    """ 

    # Check that the size of the arrays is the same for sims, massdef and env
    if (len(sims) != len(env) or len(sims) != len(massdef)):
        print('WARNING (bp.mf_sims): Input arrays have different lengths {}, {}, {}'.format(
            sims,env,massdef))
        return None

    outdir, dirz, dirplots = b.get_outdirs(env[0],dirz=dirz,outdir=outdir,sim_label=None)

    # Get labels or check the arrays
    labels = b.get_simlabels(sims,labels=labels)

    # Set up plot variables
    fig = plt.figure()
    ax = plt.subplot()
    cols = get_distinct(len(sims))
    xtit = '${\\rm log}_{10}(\\rm{M/M_{\odot}}h^{-1})$' 
    ytit = '${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dlog}_{10}M)$'
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

    xmin = 10. ; xmax = 16.
    ymin = -6.5 ; ymax = 0.
    ax.set_xlim(xmin,xmax) ;  ax.set_ylim(ymin,ymax) 

    # Loop over all the simulations to be compared
    lowestz = 999 ; files2plot = 0
    for ii, sim in enumerate(sims):
        snap = b.get_snap(snap,sim,env[ii]) #here 
        print(snap); exit()
        mnom = massdef[ii]
        prop = mnom.split('/')[1]
        print(prop) ; exit()

        if ('FOF' in mnom):
            mass = b.get_fofprop(snap,sim,env[ii],prop,Testing=Testing)
        print(mass)
        print(mnom) ; exit()        
        
        if outfil is None: continue
        files2plot += 1 ; print(outfil)

        # Read data
        mhist, nh = np.loadtxt(outfil,usecols=(0,3),unpack=True)

        # Read volume and dm
        ff = open(outfil, 'r')
        line1 = ff.readline()
        line2 = ff.readline()
        ff.close()
        volume = float(line2.split(',')[0].split('=')[-1])
        dm = float(line2.split(',')[1].split('=')[-1])

        # Calculate HMF
        print('Side of sim box = {:.2f} Mpc^3/h^3'.format(np.power(volume,1./3.)))
        if (volume>0.):
            hmf = nh/volume/dm  # In Mpc^3/h^3

            ind = np.where(hmf>0.)
            ax.plot(mhist[ind],np.log10(hmf[ind]),c=cols[ii],label=labels[ii])

    if (files2plot<1):
        print('WARNING (bp.mf_sims): No mf_sims plot made for snapshot={}'.format(snap))
        return None

    # Legend
    ax.annotate('z='+str(zz),xy=(xmax-0.17*(xmax-xmin),ymax-0.07*(ymax-ymin)))
    leg = ax.legend(loc=3, handlelength=0, handletextpad=0)
    leg.draw_frame(False)
    for ii,text in enumerate(leg.get_texts()):
        text.set_color(cols[ii])
    for item in leg.legendHandles:
        item.set_visible(False)

    plotf = dirplots+'mf_z'+str(zz)+'_sims'+str(len(sims))+'_mdef'+\
            str(len(np.unique(massdef)))+'.pdf'
    fig.savefig(plotf)

    return plotf
    

if __name__== "__main__":
    dirz = None ; outdir = None
    snap = 18
    zz = 3.

    env = 'arilega'
    #env = 'cosmalega'

    if (env == 'cosmalega'):
        sim = 'L400N1024/WMAP9/Sims/BAHAMAS'
        dirz = '/cosma6/data/dp004/dc-gonz3/BAHAMAS/'
        outdir = '/cosma6/data/dp004/dc-gonz3/Junk/'
    if (env == 'arilega'):
        sim = 'HIRES/AGN_RECAL_nu0_L100N512_WMAP9'
        dirz = '/hpcdata4/arivgonz/BAHAMAS/'
        outdir = '/hpcdata4/arivgonz/Junk/'         
    if (env == 'ari'):
        sim = 'L050N256/WMAP9/Sims/ws_324_23_mu_7_05_dT_8_35_n_75_BH_beta_1_68_msfof_1_93e11'

    print(mf_sims(zz,['FOF/Group_M_Mean200','Subhalo/ApertureMeasurements/Mass/030kpc'],
                  [sim,sim],[env,env],
                  dirz=dirz,outdir=outdir, Testing=True))

