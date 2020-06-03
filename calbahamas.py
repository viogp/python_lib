import sys,os.path
import numpy as np
import h5py
import glob
import Cosmology as cosmo
import bahamas as b
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
    >>> import calbahamas as pb
    >>> sims = ['L050N256/WMAP9_PMGRID512/Sims/ex','L050N256/WMAP9_PMGRID1024/Sims/ex'] 
    >>> labels = ['PMGRID=512','PMGRID=1024']
    >>> pb.cal_plots(sims,'cosma')
    >>> pb.cal_plots(['AGN_TUNED_nu0_L100N256_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9'],'ari')
    """ 

    # Generate labels
    labels = [x.split('/')[-1] for x in sims]

    # The subfiles to loop over
    nvols = 'All'
    if Testing: nvols = 2

    ## Bins in halo mass
    #mmin = 9. ; mmax = 16. ; dm = 0.1
    #edges = np.array(np.arange(mmin,mmax,dm))
    #mhist = edges[1:]-0.5*dm  #; print("mhist={}".format(mhist)) 
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
    ## Loop over all the simulations to be compared
    #lowestz = 999 ; files2plot = 0
    #for ii, sim in enumerate(sims):
    #    hmf  = np.zeros(shape=(len(mhist)))
    #    volume = 0.
    #
    #    # Get the closest snapshot to the input redshift
    #    zmins, zmaxs = b.get_zminmaxs([zz])
    #    snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env)
    #
    #    # Get subfind files
    #    files = b.get_subfind_files(snap,sim,env)
    #    if (len(files)<1):
    #        print('WARNING (plotbahamas): no subfind files at snap={}, {} '.format(snap,sim))
    #        continue
    #    files2plot += len(files)
    #
    #    for iff, ff in enumerate(files):
    #        f = h5py.File(ff, 'r')
    #        # Read volume in first iteration
    #        if (iff == 0):
    #            header = f['Header']
    #            volume = np.power(header.attrs['BoxSize'],3.)
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
    #        if (nvols != 'All'):
    #            if (iff>nvols): break
    #
    #    print('Side of sim box = {:.2f} Mpc^3/h^3'.format(np.power(volume,1./3.)))
    #    if (volume>0.):
    #        hmf = hmf/volume/dm  # In Mpc^3/h^3
    #
    #        ind = np.where(hmf>0.)
    #        ax.plot(mhist[ind],np.log10(hmf[ind]),c=cols[ii],label=labels[ii])
    #
    #if (files2plot<1):
    #    print('WARNING (plotbahamas): No mf_sims plot made at z={}'.format(zz))
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
    ## Path to plot
    #if (dirplot == None):
    #    dirp = b.get_dirb(env)+'plots/'+sim+'/'
    #else:
    #    dirp = dirplot+sim+'/'
    #
    #if (not os.path.exists(dirp)):
    #    os.makedirs(dirp)
    #plotf = dirp+'mf_sims_z'+str(zz)+'_'+massdef+'.pdf'
    #fig.savefig(plotf)
    plotf='blublu'
    return plotf
    

if __name__== "__main__":
    env = 'cosma'

    if (env == 'cosma'):
        sims = ['L050N256/WMAP9/Sims/ex','L050N256/WMAP9_PMGRID1024/Sims/ex']
        labels = ['PMGRID=512','PMGRID=1024']

    elif (env == 'ari'):
        sims=['AGN_TUNED_nu0_L100N256_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9']
        labels = None
        
    pb.cal_plots(sims,env,labels=labels)
