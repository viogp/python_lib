import sys,os.path
import numpy as np
import h5py
import glob
import Cosmology as cosmo
import bahamas as b
import bahamasplot as pb
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)
#print('\n \n')

def cal_plots(sims,env,zz=0.,massdef='ApertureMeasurements/Mass/030kpc',
              labels=None,dirplot=None,Testing=False):
    """
    Compare the halo mass function of different simulations at a given z

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    env : string
        ari or cosma, to use the adecuate paths
    zz : float
        Redshift to make the calibration plot, set to 0 by default.
    massdef : string
        Name of the stellar mass definition to be used
    labels : list of strings
        Array with the labels to be used
    dirplot : string
        Path to plots
    Testing : boolean
        True or False for testing with few subfiles

    Returns
    -----
    plotf : string
        Path to plot

    Examples
    ---------
    >>> import bahamascal as bc
    >>> sims = ['L050N256/WMAP9_PMGRID512/Sims/ex','L050N256/WMAP9_PMGRID1024/Sims/ex'] 
    >>> labels = ['PMGRID=512','PMGRID=1024']
    >>> bc.cal_plots(sims,'cosma')
    >>> bc.cal_plots(['AGN_TUNED_nu0_L100N256_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9'],'ari')
    """ 

    # Generate labels
    labels = pb.get_simlabels(sims,labels=labels)

    # The subfiles to loop over
    nvols = 'All'
    if Testing: nvols = 2

    # Set up ploting grid 
    fig = plt.figure()
    gs = matplotlib.gridspec.GridSpec(3,2)

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])
    ax4 = plt.subplot(gs[4])
    ax5 = plt.subplot(gs[5])
   
    # Redshift ranges to look for snapshot, zmin<z<zmax
    zmins,zmaxs = b.get_zminmaxs([zz])

    # Initialize the GSMF arrays and plot
    mmin = 8.5 ; mmax = 16. ; dm = 0.1
    edges = np.array(np.arange(mmin,mmax,dm))
    mhist = edges[1:]-0.5*dm

    xtit="${\\rm log}_{10}(M_{*}/h^{-1}{\\rm M}_{\odot})$" 
    ytit="${\\rm log}_{10}(\phi/h^3{\\rm Mpc}^{-3}{\\rm dex}^{-1})$"  
    xmin=10. ; xmax=12.
    ymin=-5. ; ymax=0. 
    ax0.set_xlim(xmin,xmax) ; ax0.set_ylim(ymin,ymax) 
    ax0.set_xlabel(xtit) ; ax0.set_ylabel(ytit)
    ax0.text(xmax-0.15*(xmax-xmin),ymax-0.05*(ymax-ymin), 'z='+str(zz))

    # Initialize the Madau plot
    zmin = 0. ; zmax = 10. ; dz = 1.
    edges = np.array(np.arange(mmin,mmax,dm))
    zhist = edges[1:]-0.5*dz

    xtit="Redshift" 
    ytit="$\\dot{\\rho}_*({\\rm M}_{\\odot}{\\rm yr}^{-1}{\\rm cMpc}^{-3})$"  
    xmin=zmin ; xmax=zmax
    ymin=0.001 ; ymax=0.4 
    ax5.set_xlim(xmin,xmax) ; ax5.set_ylim(ymin,ymax) 
    ax5.set_xlabel(xtit) ; ax5.set_ylabel(ytit)
    ax5.set_yscale('log')

    #
    ## Set up plot variables
    #fig = plt.figure()
    #ax = plt.subplot()
    #cols = get_distinct(lensims)
    #xtit = '${\\rm log}_{10}(\\rm{M/M_{\odot}}h^{-1})$' 
    #ytit = '${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dlog}_{10}M)$'
    #ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    #
    #xmin = 10. ; xmax = 16.
    #ymin = -6.5 ; ymax = 0.
    #ax.set_xlim(xmin,xmax) ;  ax.set_ylim(ymin,ymax) 
    #
    # Loop over all the simulations to be compared
    files2plot = 0
    for ii, sim in enumerate(sims):
        #hmf  = np.zeros(shape=(len(mhist)))
        volume = 0.
    
        # Get the closest snapshot to the input redshift
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env)
    
        # Get subfind files
        files = b.get_subfind_files(snap,sim,env)
        if (len(files)<1):
            print('WARNING (bahamascal): no subfind files at snap={}, {} '.format(snap,sim))
            continue
        files2plot += len(files)
    
        for iff, ff in enumerate(files):
            f = h5py.File(ff, 'r')

            # Read the BoxSize in first iteration
            if (iff == 0):
                header = f['Header']
                boxsize = header.attrs['BoxSize']
    #
    #        haloes = f['FOF']
    #        mh = haloes[massdef][:]  #10^10Msun/h
    #        ind = np.where(mh > 0.)
    #        lmh = np.log10(mh[ind]) + 10.
    #
    #        # HMF
    #        H, bins_edges = np.histogram(lmh,bins=edges)
    #        hmf[:] = hmf[:] + H
    #
            if (nvols != 'All'):
                if (iff>nvols): break
    #
    #    print('Side of sim box = {:.2f} Mpc^3/h^3'.format(boxsize)))
    #    if (boxsize>0.):
    #volume = np.power(boxsize,3.)
    #        hmf = hmf/volume/dm  # In Mpc^3/h^3
    #
    #        ind = np.where(hmf>0.)
    #        ax.plot(mhist[ind],np.log10(hmf[ind]),c=cols[ii],label=labels[ii])

        # Madau plot 
        fil_sfr = b.get_path2data(sim,env)+'sfr.txt'
        aexp, SFR = np.loadtxt(fil_sfr, usecols=(0,2),unpack=True)
        # SFR (total, end) [M0/yr] 
        xz = (1./aexp) - 1.

        #n, bin_edges, binnumber

        #print(xz) ; sys.exit()

    #if (files2plot<1):
    #    print('WARNING (bahamasplot): No mf_sims plot made at z={}'.format(zz))
    #    return ' '
    #
    ## Legend
    #ax.annotate('z='+str(zz),xy=(xmax-0.17*(xmax-xmin),ymax-0.07*(ymax-ymin)))
    #leg = ax.legend(loc=3, handlelength=0, handletextpad=0)
    #leg.draw_frame(False)
    #for ii,text in enumerate(leg.get_texts()):
    #    text.set_color(cols[ii])
    #for item in leg.legendHandles:
    #    item.set_visible(False)
    #
    # Path to plot
    if (dirplot == None):
        dirp = b.get_dirb(env)+'plots/'+sim+'/'
    else:
        dirp = dirplot+sim+'/'

    if (not os.path.exists(dirp)):
        os.makedirs(dirp)

    plotf = dirp+'cal_z'+str(zz)+'.pdf'
    fig.savefig(plotf)

    return plotf
    

if __name__== "__main__":
    env = 'cosma'

    if (env == 'cosma'):
        sims = ['L050N256/WMAP9/Sims/ws_96_84_mu_7_76_dT_7_71_n_24_BH_DensTh_m_2_76_ntask128']
        labels = None

    elif (env == 'ari'):
        sims=['AGN_TUNED_nu0_L100N256_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9']
        labels = None
        
    cal_plots(sims,env,labels=labels,Testing=True)
