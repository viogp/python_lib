import sys,os.path
import numpy as np
import h5py
import glob
import pandas as pd
import Cosmology as cosmo
import stats
import bahamas as b
import bahamasplot as pb
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
#from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)
#print('\n \n')

def simname(innom,msfof=False):
    inleg = innom.split('/')[-1]
    if not msfof:
        val = '_msfof'
        if val in inleg:
            inleg = inleg.split(val)[0]
        
    return inleg
    
    
def cal_plots(sims,env,zz=0.,massdef='ApertureMeasurements/Mass/030kpc',
              labels=None,dirplot=None,Testing=False):
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
    >>> bc.cal_plots(sims,'cosma',labels=labels)
    >>> bc.cal_plots(['AGN_TUNED_nu0_L100N256_WMAP9','HIRES/AGN_RECAL_nu0_L100N512_WMAP9'],'ari')
    """ 

    # BAHAMAS variables
    nom_sfr = 'StarFormationRate' #Msun/yr
    nom_mass = 'Mass_030kpc' #10^10Msun/h
    mtype = 'star'
    itype = b.ptypes.index(mtype)
    
    # Generate labels
    labels = pb.get_simlabels(sims,labels=labels)

    # The subfiles to loop over
    nvols = 'All'
    if Testing: nvols = 2

    # Observations
    diro = b.get_dirobs(env)
    dircalo = diro+'calibration/'

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
    mmin = 8.5 ; mmax = 16. ; dm = 0.1
    medges = np.array(np.arange(mmin,mmax,dm))
    mhist = medges[1:]-0.5*dm

    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$" 
    ytit="${\\rm log}_{10}(\phi/{\\rm Mpc}^{-3}{\\rm dex}^{-1})$"  
    xmin = 9. ; xmax = 12.5
    ymin = -5. ; ymax = -1.
    ax1.set_xlim(xmin,xmax) ;  ax1.set_ylim(ymin,ymax) 
    ax1.set_xlabel(xtit) ; ax1.set_ylabel(ytit)
    #ax1.text(xmax-0.15*(xmax-xmin),ymax-0.05*(ymax-ymin), 'z='+str(zz))

    # Initialize the SMHM relation plot
    xtit="${\\rm log}_{10}(M_{200c}/{\\rm M}_{\odot})$"
    ytit="${\\rm log}_{10}(M_{*}/M_{200c})$"
    ax2.set_xlim(10.,14.) ; ax2.set_ylim(-3,0.) 
    ax2.set_xlabel(xtit) ; ax2.set_ylabel(ytit)

    # Initialize the sSFR plot (using mass ranges from GSMF)
    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$"
    ytit="${\\rm log}_{10}({\\rm sSFR}/{\\rm Gyr}^{-1})$"
    ax3.set_xlim(xmin,xmax) ; ax3.set_ylim(-2,5.) 
    ax3.set_xlabel(xtit) ; ax3.set_ylabel(ytit)

    # Initialize the passive fraction plot (using mass ranges from GSMF)
    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$"
    ytit="Passive fraction"   
    ax4.set_xlim(xmin,xmax) ; ax4.set_ylim(0.,1.) 
    ax4.set_xlabel(xtit) ; ax4.set_ylabel(ytit)

    # Initialize the Madau plot
    zmin = 0. ; zmax = 10. ; dz = 1.
    zedges = np.array(np.arange(zmin,zmax+dz,dz))
    zhist = zedges[1:]-0.5*dz

    xtit="Age(Gyr)" 
    ytit="$\\dot{\\rho}_*({\\rm M}_{\\odot}{\\rm yr}^{-1}{\\rm cMpc}^{-3})$"  
    ax5.set_xlim(zmin,zmax) ; ax5.set_ylim(0.001,0.4) 
    ax5.set_xlabel(xtit) ; ax5.set_ylabel(ytit)
    ax5.set_yscale('log')

    # Loop over all the simulations to be compared
    files2plot = 0
    for ii, sim in enumerate(sims):
        print('Starting with sim{}: {}'.format(ii,sim))
        volume = 0.

        # Get the closest snapshot to the input redshift
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env)

        # Get the indexes for centrals
        cenids = b.cenids(snap,sim,env)
        
        # Get particle files
        files = b.get_particle_files(snap,sim,env)
        if (len(files)<1):
            print('WARNING (bahamascal): no subfind files at snap={}, {} '.format(snap,sim))
            continue

        # Look over all the files
        ntot, pftot, pfcen, pfsat  = [np.zeros(shape=(len(mhist))) for i in range(4)]
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
    
        for iff, ff in enumerate(files):
            f = h5py.File(ff, 'r')
            fof = f['FOF'] ; subhaloes = f['Subhalo']
            
            # Stellar mass
            mass1 = subhaloes[massdef][:,itype]  #10^10Msun/h
            lm1 = np.zeros(shape=len(mass1)) ; lm1.fill(-999.)
            ind = np.where(mass1 > 0.)
            lm1[ind] = np.log10(mass1[ind]) + 10. - np.log10(h0) #Msun
            
            # sSFR
            sfr1  = subhaloes[nom_sfr][:] #Msun/h/yr
            ssfr1 = np.zeros(shape=len(sfr1)) ; ssfr1.fill(-999.)
            ssfr1[ind] = sfr1[ind]/(10.*mass1[ind]) #1/Gyr

            totlen = len(mass1)
            mass1 = [] ;sfr1 = []
                
            # Read other quantities
            if (iff == 0):
                lm    = lm1
                ssfr  = ssfr1
                m200  = fof['Group_M_Crit200'][:]*1e10/h0      #Msun
                m500  = fof['Group_M_Crit500'][:]*1e10/h0      #Msun
                r500  = fof['Group_R_Crit500'][:]/h0           #cMpc
                cop_x = fof['GroupCentreOfPotential'][:,0]/h0  #cMpc
                cop_y = fof['GroupCentreOfPotential'][:,1]/h0  #cMpc
                cop_z = fof['GroupCentreOfPotential'][:,2]/h0  #cMpc
            else:
                lm    = np.append(lm,lm1)
                ssfr  = np.append(ssfr,ssfr1)
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
            ind = np.where((ssfr1<=0.3*slim) & (ssfr1 != -999.))
            if (np.shape(ind)[1]>1):
                H, bins_edges = np.histogram(lm1[ind],bins=medges)
                pftot = pftot + H

            # Get central indices
            iglob = np.arange(totlen) + istart
            icen = np.intersect1d(iglob,cenids) - istart
            istart = istart + totlen
            if (len(icen)>=1):
                if(icen[-1]>totlen-1):
                    print('STOP: something wrong with central index')

                # Passive central galaxies
                lmp = lm1[icen]
                ssfrp = ssfr1[icen]
                ind = np.where((ssfrp<=0.3*slim) & (ssfrp != -999.))
                if (np.shape(ind)[1]>1):
                    H, bins_edges = np.histogram(lmp[ind],bins=medges) 
                    pfcen = pfcen + H

            # Get satellite indices
            if (len(icen)>1):
                isat = np.array([i for i in range(totlen) if i not in icen],dtype=int)
            else:
                isat = np.arange(totlen)

            if (len(isat)>=1):
                # Passive satellite galaxies
                lmp = lm1[isat]
                ssfrp = ssfr1[isat]
                ind = np.where((ssfrp<=0.3*slim) & (ssfrp != -999.))
                if (np.shape(ind)[1]>1):
                    H, bins_edges = np.histogram(lmp[ind],bins=medges) 
                    pfsat = pfsat + H
                
            # Sample a set of subvolumes
            if (nvols != 'All'):
                if (iff>nvols): break
    
        print('Side of sim box = {:.2f} Mpc^3/h^3'.format(boxsize))
        if (boxsize<=0.):
            continue
        volume = np.power(boxsize/h0,3.) # In Mpc^3

        #--------------------------------------------------
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
        medians = stats.perc_2arrays(gedges,mass500,gas_mh,0.5,nmin=5)
        ind = np.where(medians > 0.)
        if (np.shape(ind)[1] > 0):
            ax0.plot(ghist[ind],medians[ind])
        
        # Plot individual cases where there's not enough data
        gind = np.where(medians == -999.)
        if (np.shape(gind)[1] > 0):
            val = gedges[gind[0][0]]
            ind = np.where(mass500 > val)
            ax0.scatter(mass500[ind],gas_mh[ind], s=40, zorder=10)

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

        #--------------------------------------------------
        # GSMF model
        gsmf = ntot/volume/dm  # In Msun/Mpc^3 
        ind = np.where(gsmf>0.)
        ax1.plot(mhist[ind],np.log10(gsmf[ind]))

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
                    label='Baldry*2012, z>0.06')

        if (ii==0): # Obs legend
            leg = ax1.legend(loc=0) ; leg.draw_frame(False)


        #--------------------------------------------------
        # SMHM-relation model
        if (len(cenids)>=1):
            if(cenids[-1]>len(lm)-1):
                print('STOP: something wrong with central index')

            lmc = lm[cenids]
            mhc = m200[cenids]
            ind = np.where((lmc > 0) & (mhc > 0))
            x = np.log10(mhc[ind])
            y = 10**(lmc[ind] - x)            
            medians = stats.perc_2arrays(medges,x,y,0.5,nmin=5)
            ind = np.where(medians != -999.)
            if (np.shape(ind)[1] > 0):
                ax2.plot(mhist[ind],np.log10(medians[ind]))
        
            ## Plot individual cases where there's not enough data
            #if (np.shape(ind)[1] < len(medians)):
            #    val = medges[ind[0][-1]+1]
            #    iind = np.where(x >= val)
            #    ax2.scatter(x[iind], np.log10(y[iind]), s=40, zorder=10)

        
        #--------------------------------------------------
        # sSFR model
        ind = np.where(ssfr>0.3*slim)
        medians = stats.perc_2arrays(medges,lm[ind],ssfr[ind],0.5,nmin=5)
        ind = np.where(medians != -999.)
        if (np.shape(ind)[1] > 0):
            ax3.plot(mhist[ind],np.log10(medians[ind]))

        # Plot individual cases where there's not enough data
        if (np.shape(ind)[1] < len(medians)):
            val = medges[ind[0][-1]+1]
            iind = np.where((lm >= val) & (ssfr>0.3*slim))
            ax3.scatter(lm[iind],np.log10(ssfr[iind]), s=40, zorder=10)

        # Obs
        fobs = diro+'passivef/z0_gilbank10.txt'
        #ax4.errorbar(xo,p, yerr=erro, color=ocol, ecolor=ocol,
        #             label ='Gilbank+10', fmt = 'o')
        #if (ii==0): # Obs legend
        #    leg = ax4.legend(loc=0) ; leg.draw_frame(False)

                    
        #--------------------------------------------------
        # Passive fraction 
        for i in range(len(mhist)):
            if (ntot[i]>5.):
                pftot[i] = pftot[i]/ntot[i]
                pfcen[i] = pfcen[i]/ntot[i]
                pfsat[i] = pfsat[i]/ntot[i]
            else:
                pftot[i] = -999.
                pfcen[i] = -999. ; pfsat[i] = -999.
        ax4.plot(mhist,pftot)
        if(len(sims)==1):
            ax4.plot(mhist,pfcen,linestyle='--')
            ax4.plot(mhist,pfsat,linestyle=':')
            
        # Obs
        fobs = diro+'passivef/z0_gilbank10.txt'
        lm, p, erro = np.loadtxt(fobs, unpack=True)
        xo = lm + np.log10(0.7)
        ax4.errorbar(xo,p, yerr=erro, color=ocol, ecolor=ocol,
                     label ='Gilbank+10', fmt = 'o')

        fobs = diro+'passivef/z0_bauer13.txt'
        lm, p = np.loadtxt(fobs, unpack=True)
        xo = lm + np.log10(0.7)
        erro = lm*0.
        ax4.errorbar(xo,p, yerr=erro, color=ocol, ecolor=ocol,
                    label ='Bauer+10', fmt = '^')
        

        if (ii==0): # Obs legend
            leg = ax4.legend(loc=0) ; leg.draw_frame(False)

        #--------------------------------------------------
        # Madau plot 
        fil_sfr = b.get_path2data(sim,env)+'sfr.txt'
        aexp, sfr = np.loadtxt(fil_sfr, usecols=(0,2),unpack=True)
        # SFR (total, end) [M0/yr] 
        zval = (1./aexp) - 1.

        medians = stats.perc_2arrays(zedges,zval,sfr/volume,0.5)
        ax5.plot(zhist,medians,label=simname(sims[ii]))
        #print(medians) ###HERE

    #if (files2plot<1):
    #    print('WARNING (bahamasplot): No mf_sims plot made at z={}'.format(zz))
    #    return ' '
    #
    # Legend
    if (len(sims) == 1):
        # Title when showing only one sim
        fig.suptitle(simname(sims[0])+', z='+str(zz))
    else:
        leg = ax5.legend(loc=0)
        #leg = ax5.legend(loc=0, handlelength=0, handletextpad=0)
        leg.draw_frame(False)
        #for ii,text in enumerate(leg.get_texts()):
        #    text.set_color(leg.get_color([ii]))
        #for item in leg.legendHandles:
        #    item.set_visible(False)
    
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
