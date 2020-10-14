import sys,os.path
import numpy as np
import h5py
import glob
import pandas as pd
import Cosmology as cosmo
import stats
import bahamas as b
import bahamasplot as bp
from scipy.interpolate import interp1d
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
#from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)
#print('\n \n')

npart = 100 # Minimal number of particles to trust an attribute
npartsf = 1. # Minimal number of SF particles to trust an attribute

lw = 4
dcol = 'k'
dsty = ['-','--']
dbig = [3,3]

def compare_dsims(sims,envs,massdefs,highres=True):
    """
    Add default Bahamas sims to compare with new runs

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    envs : list of string
        cosma, ari or arilega, to use the adecuate paths
    massdefs : list of strings
        Array with mass definition names
    highres : boolean
        True or False for including the HighRes Bahamas or not

    Returns
    -----
    sims : list of strings
        Updated array with the names of the simulations
    envs : list of strings
        Updated array with the arilega of the default simularions
    massdefs : list of strings
        Updated array with mass names
    """ 

    if (envs[0] == 'ari'):
        massd = 'Mass_'+massdefs[0].split('/')[-1]
        
        if (highres):
            sims.insert(0,'HIRES/AGN_RECAL_nu0_L100N512_WMAP9')
            envs.insert(0,'arilega')
            massdefs.insert(0,massd)

        sims.insert(0,'AGN_TUNED_nu0_L400N1024_WMAP9')
        envs.insert(0,'arilega')
        massdefs.insert(0,massd)

        return sims,envs,massdefs
    else:
        print('WARNING (bahamascal): Comparisson with default Bahamas to be done within ARI')
        return sims,envs,massdefs

def cal_plots(sims,env,zz=0.,massdef='ApertureMeasurements/Mass/030kpc',
              ndatbin=5,labels=None,dirplot=None,
              compare_default=False,Testing=False):
    """
    Compare the halo mass function of different simulations at a given z

    Parameters
    -----------
    sims : list of strings
        Array with the names of the simulation
    env : string
        cosma, ari or arilega, to use the adecuate paths
    zz : float
        Redshift to make the calibration plot, set to 0 by default.
    massdef : string
        Name of the stellar mass definition to be used
    ndatbin : integer
        Minimum number of data points per bin to calculate medians or means.

    labels : list of strings
        Array with the labels to be used
    dirplot : string
        Path to plots
    compare_default : boolean
        True or False for comparing with published Bahamas runs. Need to be in ari env.
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
    >>> bc.cal_plots(sims,'cosma',labels=labels)
    >>> bc.cal_plots(['AGN_TUNED_nu0_L100N256_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9'],'ari')
    """ 

    # BAHAMAS variables
    nom_sfr = 'StarFormationRate' #Msun/yr
    nom_mass = 'Mass_030kpc' #10^10Msun/h
    mtype = 'star'
    itype = b.ptypes.index(mtype)

    if (labels == None):
        # Generate labels
        labels = bp.get_simlabels(sims,labels=labels)

    # The subfiles to loop over
    nvols = 'All'
    if Testing: nvols = 2

    # Observations
    diro = b.get_dirobs(env)
    dircalo = diro+'calibration/'

    # Comparison with published Bahamas if in ari env.
    nsims = len(sims)
    envs = [env for i in range(nsims)]
    massdefs = [massdef for i in range(nsims)]

    #here deal with massdefs
    if (compare_default):
        sims, envs, massdefs = compare_dsims(sims,envs,massdefs)

    # Set up ploting grid 
    fig = plt.figure(figsize=(14.,21.))
    gs = matplotlib.gridspec.GridSpec(3,2)

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])
    ax4 = plt.subplot(gs[4])
    ax5 = plt.subplot(gs[5])

    ocol = 'grey'
    cm = plt.get_cmap('tab10') # Colour map to draw colours from

    # Redshift ranges to look for snapshot, zmin<z<zmax
    zmins,zmaxs = b.get_zminmaxs([zz])

    # Initialize the fgas plot
    gmin = 13. ; gmax = 16. ; dg = 0.2
    gedges = np.array(np.arange(gmin,gmax+dg,dg))
    ghist = gedges[1:]-0.5*dg

    xtit="${\\rm log}_{10}(M_{500}/{\\rm M}_{\odot})$" 
    ytit="$M_{\\rm gas,500}/M_{500}$"  
    xmin = 13. ; xmax = 15.2
    ymin = 0. ; ymax = 0.2
    ax0.set_xlim(xmin,xmax) ;  ax0.set_ylim(ymin,ymax) 
    ax0.set_xlabel(xtit) ; ax0.set_ylabel(ytit)
    #ax0.text(xmin+0.15*(xmax-xmin),ymax-0.05*(ymax-ymin), 'z='+str(zz))    

    # Initialize the GSMF arrays and plot
    mmin = 8.5 ; mmax = 16. ; dm = 0.2
    medges = np.array(np.arange(mmin,mmax,dm))
    mhist = medges[1:]-0.5*dm

    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$" 
    ytit="${\\rm log}_{10}(\phi/{\\rm Mpc}^{-3}{\\rm dex}^{-1})$"  
    xmin = 9. ; xmax = 12.2
    ymin = -5. ; ymax = -1.
    ax1.set_xlim(xmin,xmax) ;  ax1.set_ylim(ymin,ymax) 
    ax1.set_xlabel(xtit) ; ax1.set_ylabel(ytit)
    #ax1.text(xmax-0.15*(xmax-xmin),ymax-0.05*(ymax-ymin), 'z='+str(zz))

    # Initialize the SMHM relation plot
    xtit="${\\rm log}_{10}(M_{200c}/{\\rm M}_{\odot})$"
    ytit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$"  #ytit="${\\rm log}_{10}(M_{*}/M_{200c})$"
    ax2.set_xlim(10.5,14.) ; ax2.set_ylim(8.,xmax) #; ax2.set_ylim(-3,0.) 
    ax2.set_xlabel(xtit) ; ax2.set_ylabel(ytit)

    # Initialize the sSFR plot (using mass ranges from GSMF)
    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$"
    ytit="${\\rm log}_{10}({\\rm SFR}/M_{*}/{\\rm Gyr}^{-1})$"
    ymin_ssfr = -2 ; ymax_ssfr = 1.
    ax3.set_xlim(xmin,xmax) ; ax3.set_ylim(ymin_ssfr,ymax_ssfr) 
    ax3.set_xlabel(xtit) ; ax3.set_ylabel(ytit)

    # Initialize the passive fraction plot (using mass ranges from GSMF)
    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$"
    ytit="Passive fraction"   
    ax4.set_xlim(xmin,xmax) ; ax4.set_ylim(0.,1.) 
    ax4.set_xlabel(xtit) ; ax4.set_ylabel(ytit)

    # Initialize the Madau plot
    zedges = np.logspace(0, 1, num=15)
    zhist = (zedges[1:]+zedges[:-1])/2

    xtit="1+z" 
    ytit="$\\dot{\\rho}_*({\\rm M}_{\\odot}{\\rm yr}^{-1}{\\rm cMpc}^{-3})$"  
    ax5.set_xlim(1,11) ; ax5.set_ylim(0.001,0.4) 
    ax5.set_xlabel(xtit) ; ax5.set_ylabel(ytit)
    ax5.set_xscale('log') ; ax5.set_yscale('log')
    ax5.xaxis.set_major_formatter(ticker.FuncFormatter(bp.logformat)) 
    ax5.yaxis.set_major_formatter(ticker.FuncFormatter(bp.logformat)) 
    
    # Loop over all the simulations to be compared
    files2plot = 0 ; cols =[]
    for ii, sim in enumerate(sims):
        env = envs[ii] ; massdef = massdefs[ii] 
        print('Starting with sim{}: {} ({}, {})'.format(ii,sim,env,massdef))
        volume = 0.
        col = cm(1.*ii/nsims) ; cols.append(col)

        # Get the closest snapshot to the input redshift
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env)

        # Get particle files
        files = b.get_particle_files(snap,sim,env)
        if (len(files)<1):
            print('WARNING (bahamascal): no subfind files at snap={}, {} '.format(snap,sim))
            continue

        # Get mass resolution
        mdm,mgas = b.resolution(sim,snap,env)
        
        # Look over all the files
        ntot, pftot, pfcen, pfsat, pfres  = [np.full((len(mhist)),0.) for i in range(5)]

        istart = 0
        for iff, ff in enumerate(files):
            f = h5py.File(ff, 'r')
            p0 = f['PartType0'] # Gas particles (0:gas, 1:DM, 4: stars, 5:BH)

            if (iff == 0):
                # Read simulation constants in first iteration
                header = f['Header'] #;print(list(header.attrs.items()))
                boxsize = header.attrs['BoxSize'] #Mpc/h
                h0 = header.attrs['HubbleParam'] 
                omega0 = header.attrs['Omega0']
                omegab = header.attrs['OmegaBaryon']
                lambda0 =header.attrs['OmegaLambda']
                cosmo.set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0,h0=h0,
                              universe="Flat",include_radiation=False)
                slim = 1./cosmo.tHubble(z_snap) #1/Gyr 

                # Read gas properties of the particles
                groupnum = p0['GroupNumber'][:]
                subgroupnum = p0['SubGroupNumber'][:]
                partmass = p0['Mass'][:]*1e10/h0       #Msun
                partx = p0['Coordinates'][:,0]/h0 
                party = p0['Coordinates'][:,1]/h0 
                partz = p0['Coordinates'][:,2]/h0 

            else:
                # Read gas properties of the particles
                groupnum    = np.append(groupnum,p0['GroupNumber'][:])
                # Negative values: particles that don't belong to a halo
                subgroupnum = np.append(subgroupnum,p0['SubGroupNumber'][:])
                partmass    = np.append(partmass,p0['Mass'][:]*1e10/h0)
                partx       = np.append(partx,p0['Coordinates'][:,0]/h0) 
                party       = np.append(party,p0['Coordinates'][:,1]/h0)
                partz       = np.append(partz,p0['Coordinates'][:,2]/h0)
                
            f.close()

        # Get subfind files
        files = b.get_subfind_files(snap,sim,env)
        if (len(files)<1):
            print('WARNING (bahamascal): no subfind files at snap={}, {} '.format(snap,sim))
            continue
        files2plot += len(files)

        minsfr = 1e10
        for iff, ff in enumerate(files):
            f = h5py.File(ff, 'r') 
            fof = f['FOF'] ; subhaloes = f['Subhalo']

            if (env != 'arilega'):
                snum1  = subhaloes['SubGroupNumber'][:]
            
            # Stellar mass
            mass1 = subhaloes[massdef][:,itype]  #10^10Msun/h
            lm1 = np.zeros(shape=len(mass1)) ; lm1.fill(-999.)
            ind = np.where(mass1 > 0.)
            lm1[ind] = np.log10(mass1[ind]) + 10. - np.log10(h0) #Msun
            
            # sSFR
            sfr1  = subhaloes[nom_sfr][:] #Msun/h/yr
            ssfr1 = np.zeros(shape=len(sfr1)) ; ssfr1.fill(-999.)
            ssfr1[ind] = sfr1[ind]/(10.*mass1[ind]) #1/Gyr
            
            ind = np.where(sfr1>0.)
            if (np.shape(ind)[1] > 0.):
                sfGyr = sfr1[ind]*10.**9/h0  #Msun/Gyr
                if (minsfr > min(sfGyr)): minsfr = min(sfGyr)            

            mass1 = [] ; sfr1 = [] ; sfGyr = []
                
            # Read other quantities
            if (iff == 0):
                lm    = lm1
                ssfr  = ssfr1
                if (env != 'arilega'): snum  = snum1
                gnum  = subhaloes['GroupNumber'][:]

                m200  = fof['Group_M_Crit200'][:]*1e10/h0      #Msun
                m500  = fof['Group_M_Crit500'][:]*1e10/h0      #Msun
                r500  = fof['Group_R_Crit500'][:]/h0           #cMpc
                cop_x = fof['GroupCentreOfPotential'][:,0]/h0  #cMpc
                cop_y = fof['GroupCentreOfPotential'][:,1]/h0  #cMpc
                cop_z = fof['GroupCentreOfPotential'][:,2]/h0  #cMpc
            else:
                lm    = np.append(lm,lm1)
                ssfr  = np.append(ssfr,ssfr1)
                if (env != 'arilega'): snum  = np.append(snum,snum1)
                gnum  = np.append(gnum,subhaloes['GroupNumber'][:])
                
                m200  = np.append(m200,fof['Group_M_Crit200'][:]*1e10/h0)
                m500  = np.append(m500,fof['Group_M_Crit500'][:]*1e10/h0)
                r500  = np.append(r500,fof['Group_R_Crit500'][:]/h0)
                cop_x = np.append(cop_x,fof['GroupCentreOfPotential'][:,0]/h0)
                cop_y = np.append(cop_y,fof['GroupCentreOfPotential'][:,1]/h0)
                cop_z = np.append(cop_z,fof['GroupCentreOfPotential'][:,2]/h0)
            
            f.close()
            
            # Total number of galaxies per stellar bin
            H, bins_edges = np.histogram(lm1,bins=medges)
            ntot = ntot + H
            
            # Total number of passive galaxies
            ind = np.where((ssfr1 <= 0.3*slim) & (ssfr1 > -999.))
            if (np.shape(ind)[1]>1):
                H, bins_edges = np.histogram(lm1[ind],bins=medges)
                pftot = pftot + H
            
            if (env != 'arilega'):
                # Passive central galaxies
                ind = np.where((ssfr1 <= 0.3*slim) & (ssfr1 > -999.) &
                               (snum1 == 0))
                if (np.shape(ind)[1] > 0): 
                    H, bins_edges = np.histogram(lm1[ind],bins=medges) 
                    pfcen = pfcen + H

                # Passive satellite galaxies
                ind = np.where((ssfr1 <= 0.3*slim) & (ssfr1 > -999.) &
                               (snum1 > 0))
                if (np.shape(ind)[1] > 0): 
                    H, bins_edges = np.histogram(lm1[ind],bins=medges) 
                    pfsat = pfsat + H
                
            # Sample a set of subvolumes
            if (nvols != 'All'):
                if (iff>nvols): break

        print('Side of sim box = {:.2f} Mpc^3/h^3'.format(boxsize))
        if (boxsize<=0.):
            continue
        volume = np.power(boxsize/h0,3.) # In Mpc^3

        if (nsims == 1 and env != 'arilega'):
            # Number density for ngal galaxies within the given box
            ngal = 100. ;  ndlim = ngal/volume
            print('* Number density for {} gal. = {:.2e}'.format(ngal,ndlim))

        #--------------------------------------------------
        # fgas observations (1E13 Msun)
        file = diro+'calibration/all_fgas.txt'
        om500, omin, omax, ofgas, ofmax, ofmin = np.loadtxt(file, usecols=[0,1,2,3,4,5], unpack=True)
        ox = np.log10(om500) + 13.
        oxmin = np.log(omin) + 13.
        oxmax = np.log(omax) + 13.
        #ax0.errorbar(ox, ofgas, yerr=[ofmin,ofmax], xerr=[oxmin,oxmax],
        #             fmt='o',ecolor=ocol,color=ocol,mec=ocol,alpha=0.75)
        ax0.errorbar(ox, ofgas, yerr=[ofmin,ofmax],
                     fmt='o',ecolor=ocol,color=ocol,mec=ocol,alpha=0.75)

        
        # fgas model
        df_part = pd.DataFrame(data=np.vstack([groupnum,subgroupnum,partmass,
                                               partx,party,partz]).T,
                               columns=['groupnum','subgroupnum','partmass',
                                        'partx','party','partz'])
        groupnum,subgroupnum,partmass,partx,party,partz=[[] for i in range(6)]
        df_part.sort_values(by=['groupnum', 'subgroupnum'], inplace=True)
        df_part.reset_index(inplace=True, drop=True)

        df_fof = pd.DataFrame(data=np.vstack([m500,r500,cop_x,cop_y,cop_z]).T,
                              columns=['m500','r500','cop_x','cop_y','cop_z'])
        m500,r500,cop_x,cop_y,cop_z=[[] for i in range(5)]
        df_fof.index += 1
        df_fof.index.names = ['groupnum']
        df_fof.reset_index(inplace=True)

        merge = pd.merge(df_part, df_fof, on=['groupnum'])

        # Positions of gas particles relative to the center of the group
        lbox2 = boxsize/2.
        merge['partx'] = merge.partx - merge.cop_x
        merge['party'] = merge.party - merge.cop_y
        merge['partz'] = merge.partz - merge.cop_z

        # Correct for periodic boundary conditions (for gal. in groups)
        merge.partx.loc[merge.partx < -lbox2] = merge.partx.loc[merge.partx < -lbox2] + boxsize
        merge.party.loc[merge.party < -lbox2] = merge.party.loc[merge.party < -lbox2] + boxsize
        merge.partz.loc[merge.partz < -lbox2] = merge.partz.loc[merge.partz < -lbox2] + boxsize

        merge.partx.loc[merge.partx >= lbox2] = merge.partx.loc[merge.partx >= lbox2] - boxsize
        merge.party.loc[merge.party >= lbox2] = merge.party.loc[merge.party >= lbox2] - boxsize
        merge.partz.loc[merge.partz >= lbox2] = merge.partz.loc[merge.partz >= lbox2] - boxsize

        # Distances within groups and clusters
        merge = merge.loc[merge.m500 > 1e13] 
        merge['distance'] = (merge.partx**2 + 
                             merge.party**2 + 
                             merge.partz**2) ** 0.5

        # Gas mass enclosed in r500
        merge['inside_r500'] = merge.distance <= merge.r500
        merge = merge.loc[merge.inside_r500 == True]
        groups = merge.groupby(['groupnum'], as_index=False)
        massinr500 = groups.partmass.sum() # partmass now = gas mass
        final = pd.merge(massinr500, df_fof, on=['groupnum'])

        # Plot median gas_mass(within r500)/m500 vs m500
        mass500 = np.log10(final.m500.values)
        gas_mh = 10**(np.log10(final.partmass.values)-mass500)

        if (nsims == 1 and env != 'arilega'): #Quartiles
            per1 = stats.perc_2arrays(gedges,mass500,gas_mh,0.1,nmin=ndatbin)
            per9 = stats.perc_2arrays(gedges,mass500,gas_mh,0.9,nmin=ndatbin)
            ind = np.where((per1 != -999.) & (per9 != -999.))
            if (np.shape(ind)[1] > 0):
                ax0.fill_between(ghist[ind],per1[ind],per9[ind],alpha=0.2,color=col)
        
        medians = stats.perc_2arrays(gedges,mass500,gas_mh,0.5,nmin=ndatbin)
        ind = np.where(medians != -999.)
        if (np.shape(ind)[1] > 0):
            ax0.plot(ghist[ind],medians[ind],color=col)

            if (np.shape(ind)[1] < len(medians)):
                # Plot individual cases where there's not enough data
                val = gedges[ind[0][-1]+1]
                ind = np.where(mass500 >= val)
                ax0.scatter(mass500[ind], gas_mh[ind], s=40, zorder=10,color=col)

        #--------------------------------------------------
        # GSMF observations
        file = diro+'gsmf/baldry_2012_z0_cha.txt'
        xobs, p3, dp3 = np.loadtxt(file, usecols=[0,1,2], unpack=True)
        yobs = np.zeros(shape=(len(xobs))) ; yobs.fill(-999.)
        ind = np.where(p3>0)
        yobs[ind] = np.log10(p3[ind]) -3.
        lerr = np.zeros(shape=(len(xobs))) ; lerr.fill(-999.)
        ind = np.where((p3-dp3)>0)
        lerr[ind] = np.log10(p3[ind] - dp3[ind])-3.
        herr = np.zeros(shape=(len(xobs))) ; herr.fill(-999.)
        ind = np.where((p3+dp3)>0)
        herr[ind] = np.log10(p3[ind] + dp3[ind])-3.
        ax1.errorbar(xobs,yobs, yerr=[yobs-lerr,herr-yobs], fmt='o',
                    ecolor=ocol,color=ocol,mec=ocol,
                    label='Baldry+2012, z>0.06')

        if (ii==0): # Obs legend
            leg = ax1.legend(loc=3) ; leg.draw_frame(False)

        # GSMF model
        gsmf = ntot/volume/dm  # In Msun/Mpc^3 
        sels = [((gsmf>0) & (mhist >= np.log10(npart*mgas)) & (ntot>=ndatbin)),
                ((gsmf>0) & (mhist >= np.log10(npart*mgas)) & (ntot<=ndatbin)),
                ((gsmf>0) & (mhist <= np.log10(npart*mgas)))]
        lsty = ['-','-.',':']
        for isel,sel in enumerate(sels):
            ind = np.where(sel)
            x = mhist[ind] ; y = np.log10(gsmf[ind])
            if (env == 'arilega'):
                if (isel==0):
                    styl = dsty[isel]
                else:
                    styl = lsty[isel]
                ax1.plot(x,y,color=dcol,linestyle=styl,linewidth=dbig[ii])
            else:
                ax1.plot(x,y,color=col,linestyle=lsty[isel])

        if (nsims == 1 and env != 'arilega'):
            # Find the stellar mass corresponding to ndlim
            x = np.cumsum(gsmf[::-1])[::-1]
            y=medges[:-1]
            if (ndlim >= min(x) and ndlim <= max(x)):
                f = interp1d(x,y)
                m_ndlim = f(ndlim)
            elif (ndlim > max(x)):
                m_ndlim = -999.
            elif (ndlim < min(x)):
                m_ndlim = 999.
            print('* Stellar mass corresponding to ndlim = {:.4f}'.format(m_ndlim)) 

            # Region where there are two few massive objects
            if (m_ndlim < xmax):
                ax1.axvspan(m_ndlim, xmax+1, facecolor='0.2', alpha=0.3)
            
        #--------------------------------------------------
        # SMHM-relation model
        if(env != 'arilega'):
           ind = np.where((snum==0) & (lm>-999.))
           if (np.shape(ind)[1]>0):
               lmc = lm[ind]
               mhc = m200[gnum[ind]-1]
               
               ind = np.where(mhc > 0.)
               if (np.shape(ind)[1]>0):
                   x = np.log10(mhc[ind])
                   #y = 10**(lmc[ind] - x)
                   y = 10**(lmc[ind])
           
                   if (nsims == 1 and env != 'arilega'): #Quartiles
                       per1 = stats.perc_2arrays(medges,x,y,0.1,nmin=ndatbin)
                       per9 = stats.perc_2arrays(medges,x,y,0.9,nmin=ndatbin)
           
                       ind = np.where((per1 != -999.) & (per9 != -999.))
                       if (np.shape(ind)[1] > 0):
                           ax2.fill_between(mhist[ind],
                                            np.log10(per1[ind]),
                                            np.log10(per9[ind]),alpha=0.2,color=col)
           
                   # Median values
                   medians = stats.perc_2arrays(medges,x,y,0.5,nmin=ndatbin)
                   ind = np.where((medians != -999.) & (mhist >= np.log10(npart*mdm)))
                   if (np.shape(ind)[1] > 0):
                       ax2.plot(mhist[ind],np.log10(medians[ind]),color=col) #here
                       
                       if (np.shape(ind)[1] < len(medians)):
                           # Plot individual cases where there's not enough data
                           val = medges[ind[0][-1]+1]
                           ind = np.where(x >= val)
                           ax2.scatter(x[ind], np.log10(y[ind]), s=40, zorder=10,color=col)

                   ind1 = np.where((medians != -999.) & (mhist <= np.log10(npart*mdm)))
                   if (np.shape(ind1)[1] > 0):
                       ind = ind1 #ind1.append(ind[0][0]) #here
                       ax2.plot(mhist[ind],np.log10(medians[ind]),color=col,linestyle=':')
                           
        #--------------------------------------------------
        # sSFR Obs
        fobs = diro+'calibration/G2010.dat'
        lgmstar, vmax, compl, sfr_haBD, sn = np.loadtxt(fobs, unpack=True,
                                                                 usecols=[8,9,10,14,30])
        hobs = 0.7
        corrsf = 1.49/1.57 # Kroupa to Chabrier IMF
        corrm = 0.74/0.81

        weight = np.zeros(shape=(len(vmax)))
        ind = np.where((vmax>0) & (compl>0))
        if (np.shape(ind)[1]>0):
            weight[ind] = 1.0 / (vmax[ind]/hobs**3) / compl[ind]

        lg_ms =  lgmstar[weight>0] + np.log10(corrm)
        lgssfr = np.log10(sfr_haBD) + np.log10(corrsf) - lgmstar + 9.0
        lgssfr[sn<2.5] = -999.        
        lg_ssfr = lgssfr[weight>0]

        w = weight[weight > 0]
        ok =  (10**lg_ssfr > 0.3*slim)
        x = lg_ms[ok] ; y = 10**lg_ssfr[ok]
        
        per1 = stats.perc_2arrays(medges,x,y,0.1,nmin=ndatbin)
        per9 = stats.perc_2arrays(medges,x,y,0.9,nmin=ndatbin)
        medians = stats.perc_2arrays(medges,x,y,0.5,nmin=ndatbin)

        ind = np.where(medians != -999.)
        if (np.shape(ind)[1] > 0): 
            ax3.errorbar(medges[ind], np.log10(medians[ind]),
                         yerr=[np.log10(medians[ind])-np.log10(medians[ind]-per1[ind]),
                               np.log10(medians[ind]+per9[ind])-np.log10(medians[ind])],
                         fmt='*', color=ocol)

        if (ii==0): # Obs legend
            ax3.text(xmin+(xmax-xmin)*0.05,
                     ymin_ssfr+(ymax_ssfr-ymin_ssfr)*0.05,
                     r'* Gilbank et al. 2010, H$_{\alpha}$',{'color': ocol})

        # sSFR limit used
        ax3.plot([xmin-1,xmax+1],[np.log10(0.3*slim),np.log10(0.3*slim)],color=ocol)

        if (nsims == 1 and env != 'arilega'):
            # Region where resolution might affect the results
            if (m_ndlim < xmax):
                ax4.axhspan(m_ndlim, xmax+1, facecolor='0.2', alpha=0.3)
        
        # sSFR model
        ind = np.where(ssfr>0.3*slim)
        x = lm[ind] ; y = ssfr[ind]
        
        if (nsims == 1 and env != 'arilega'): #Quartiles
            per1 = stats.perc_2arrays(medges,x,y,0.1,nmin=ndatbin)
            per9 = stats.perc_2arrays(medges,x,y,0.9,nmin=ndatbin)
            ind = np.where((per1 != -999.) & (per9 != -999.))
            if (np.shape(ind)[1] > 0):
                if (env == 'arilega'):
                    ax3.fill_between(mhist[ind],np.log10(per1[ind]),
                                     np.log10(per9[ind]),alpha=0.2,
                                     color=dcol,linestyle=dsty[ii],linewidth=dbig[ii])
                else:
                    ax3.fill_between(mhist[ind],np.log10(per1[ind]),
                                     np.log10(per9[ind]),alpha=0.2,color=col)
        
        medians = stats.perc_2arrays(medges,x,y,0.5,nmin=ndatbin)
        ind = np.where((medians != -999.) & (mhist >= np.log10(npart*mgas)))
        if (np.shape(ind)[1] > 0):
            if (env == 'arilega'):
                ax3.plot(mhist[ind],np.log10(medians[ind]),label=labels[ii],
                         color=dcol,linestyle=dsty[ii],linewidth=dbig[ii])
            else:
                ax3.plot(mhist[ind],np.log10(medians[ind]),
                         label=labels[ii],color=col)

            if (np.shape(ind)[1] < len(medians)):
                # Plot individual cases where there's not enough data
                val = medges[ind[0][-1]+1]
                ind = np.where(x >= val)
                if (env == 'arilega'):
                    ax3.scatter(x[ind],np.log10(y[ind]), s=40, zorder=10,color=dcol)
                else:
                    ax3.scatter(x[ind],np.log10(y[ind]), s=40, zorder=10,color=col)


        ind = np.where((medians != -999.) & (mhist <= np.log10(npart*mgas)))
        if (np.shape(ind)[1] > 0):
            if (env == 'arilega'):
                ax3.plot(mhist[ind],np.log10(medians[ind]),label=labels[ii],
                         color=dcol,linestyle=':',linewidth=dbig[ii])
            else:
                ax3.plot(mhist[ind],np.log10(medians[ind]),
                         label=labels[ii],color=col,linestyle=':')

        #--------------------------------------------------           
        # Obs Passive fraction
        fobs = diro+'passivef/z0_gilbank10.txt'
        lmo, p, erro = np.loadtxt(fobs, unpack=True)
        xo = lmo + np.log10(0.7)
        ax4.errorbar(xo,p, yerr=erro, color=ocol, ecolor=ocol,
                     label ='Gilbank+10', fmt = 'o')

        fobs = diro+'passivef/z0_bauer13.txt'
        lmo, p = np.loadtxt(fobs, unpack=True)
        xo = lmo + np.log10(0.7)
        erro = lmo*0.
        ax4.errorbar(xo,p, yerr=erro, color=ocol, ecolor=ocol,
                    label ='Bauer+11', fmt = '^')
        
        if (ii==0): # Obs legend
            leg = ax4.legend(loc=0) ; leg.draw_frame(False)

        if (nsims == 1 and env != 'arilega'):
            # Region where there are two few massive objects
            if (m_ndlim < xmax):
                ax4.axvspan(m_ndlim, xmax+1, facecolor='0.2', alpha=0.3)

        # Model Passive fraction

        # Accounting for resolution
        ind = np.where((ssfr <= 0.3*slim) & (ssfr > -999.) &
                       (lm >= np.log10(npartsf*minsfr/(0.3*slim))))
        if (np.shape(ind)[1]>1):
            H, bins_edges = np.histogram(lm[ind],bins=medges)
            pfres = pfres + H

        for i in range(len(mhist)):
            if (ntot[i]>ndatbin):
                pfres[i] = pfres[i]/ntot[i]
                pftot[i] = pftot[i]/ntot[i]
                pfcen[i] = pfcen[i]/ntot[i]
                pfsat[i] = pfsat[i]/ntot[i]
            else:
                pftot[i] = -999. ; pfres[i] = -999.
                pfcen[i] = -999. ; pfsat[i] = -999.

        ind = np.where(pftot>-999.)
        if (np.shape(ind)[1]>1):
            if (env == 'arilega'):
                ax4.plot(mhist[ind],pftot[ind],
                         color=dcol,linestyle=dsty[ii],linewidth=dbig[ii])
            else:
                ax4.plot(mhist[ind],pftot[ind],linestyle=':',linewidth=lw,color=col)
                iind = np.where(pfres > 0.)
                ax4.plot(mhist[iind[0][1:]],pfres[iind[0][1:]],linestyle='-',linewidth=lw,color=col)

            if(nsims == 1):
                ax4.plot(mhist[ind],pfcen[ind],linestyle='-',color=col)
                ax4.plot(mhist[ind],pfsat[ind],linestyle='--',color=col)

        #--------------------------------------------------
        # Madau plot: obs
        filo = dircalo + 'sfr_burgarella2013.txt' # SFRD_Tot [10^-2 M_sun yr^-1 Mpc^-3]
        hobs = 0.7
        corrcosmo = np.log10(h0) - np.log10(hobs)
        z,rhostar,err = np.loadtxt(filo, usecols=(0,10,11), unpack=True)
        corr = 0.76/1.20 # From Salpeter to Chabrier (Table B1, Lacey 2016)
        lrhostar = np.log10(rhostar * 1e-2 * corr) + corrcosmo 
        err = err*1e-2*corr
        ax5.errorbar(1+z, 10**lrhostar, yerr=err, fmt='o', color=ocol, label='Burgarella+2013')

        if (ii==0): # Obs legend
            leg = ax5.legend(loc=0) ; leg.draw_frame(False)
        
        # Madau plot: model
        if (env != 'arilega'):
            fil_sfr = b.get_path2data(sim,env)+'sfr.txt' # SFR(total, end)[M0/yr]
            aexp, sfr = np.loadtxt(fil_sfr, usecols=(0,2),unpack=True)
            zval = (1./aexp) #1+z

            per1 = stats.perc_2arrays(zedges,zval,sfr/volume,0.1,nmin=ndatbin)
            per9 = stats.perc_2arrays(zedges,zval,sfr/volume,0.9,nmin=ndatbin)
            ind = np.where((per1 != -999.) & (per9 != -999.))
            if (np.shape(ind)[1] > 0):
                ax5.fill_between(zhist[ind],per1[ind],per9[ind],alpha=0.2,color=col)
            
            medians = stats.perc_2arrays(zedges,zval,sfr/volume,0.5,nmin=ndatbin)
            ax5.plot(zhist,medians,label=labels[ii],linewidth=lw,color=col)


    #if (files2plot<1):
    #    print('WARNING (bahamasplot): No mf_sims plot made at z={}'.format(zz))
    #    return ' '
    #

    #leg = ax3.legend(loc=0,fontsize='small')
    leg = ax3.legend(loc=0,fontsize='small', handlelength=0, handletextpad=0)
    leg.draw_frame(False)
    for ii,text in enumerate(leg.get_texts()):
        text.set_color(cols[ii])
    for item in leg.legendHandles:
        item.set_visible(False)
    
    # Path to plot
    if (dirplot == None):
        dirp = b.get_dirb(env)+'plots/'+sim+'/'
    else:
        dirp = dirplot+sim+'/'

    if (not os.path.exists(dirp)):
        os.makedirs(dirp)

    plt.tight_layout()
    plotf = dirp+'cal_z'+str(zz)+'.pdf'
    fig.savefig(plotf)

    return plotf
    

if __name__== "__main__":
    env = 'cosma'

    if (env == 'cosma'):
        sims = ['L050N256/WMAP9/Sims/ws_108_35_mu_6_20_dT_8_13_n_74_BH_beta_1_20_msfof_1_93e11']
        labels = None

    elif (env == 'ari'):
        sims=['L050N256/WMAP9/Sims/ws_269_70_mu_3_13_dT_8_15_n_152_beta_2_67_mfof_0_5']
        labels = None
        
    print(cal_plots(sims,env,labels=labels,Testing=True))
