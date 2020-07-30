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
   
    # Redshift ranges to look for snapshot, zmin<z<zmax
    zmins,zmaxs = b.get_zminmaxs([zz])

    # Initialize the fgas plot
    gmin = 12.5 ; gmax = 16. ; dg = 0.1
    xtit="${\\rm log}_{10}(M_{500}/{\\rm M}_{\odot})$" 
    ytit="$M_{\\rm gas,500}/M_{500}$"  
    xmin = 13. ; xmax = 15.2
    ymin = 0. ; ymax = 0.2
    ax0.set_xlim(xmin,xmax) ;  ax0.set_ylim(ymin,ymax) 
    ax0.set_xlabel(xtit) ; ax0.set_ylabel(ytit)
    ax0.text(xmin+0.15*(xmax-xmin),ymax-0.05*(ymax-ymin), 'z='+str(zz))    

    # Initialize the GSMF arrays and plot
    mtype = 'star' 
    itype = b.ptypes.index(mtype) 

    mmin = 8.5 ; mmax = 16. ; dm = 0.1
    medges = np.array(np.arange(mmin,mmax,dm))
    mhist = medges[1:]-0.5*dm

    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$" 
    ytit="${\\rm log}_{10}(\phi/{\\rm Mpc}^{-3}{\\rm dex}^{-1})$"  
    xmin = 9. ; xmax = 12.5
    ymin = -5. ; ymax = -1.
    ax1.set_xlim(xmin,xmax) ;  ax1.set_ylim(ymin,ymax) 
    ax1.set_xlabel(xtit) ; ax1.set_ylabel(ytit)
    ax1.text(xmax-0.15*(xmax-xmin),ymax-0.05*(ymax-ymin), 'z='+str(zz))

    # Initialize the SMHM relation plot
    xtit="${\\rm log}_{10}(M_{200c}/{\\rm M}_{\odot})$"
    ytit="${\\rm log}_{10}(M_{*}/M_{200c})$"
    ax2.set_xlim(xmin,xmax) ; ax2.set_ylim(-2,1.) 
    ax2.set_xlabel(xtit) ; ax2.set_ylabel(ytit)

    # Initialize the sSFR plot (using mass ranges from GSMF)
    xtit="${\\rm log}_{10}(M_{*}/{\\rm M}_{\odot})$"
    ytit="${\\rm log}_{10}({\\rm sSFR}/{\\rm Gyr}^{-1})$"
    ax3.set_xlim(xmin,xmax) ; ax3.set_ylim(-2,10.) 
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
        gsmf  = np.zeros(shape=(len(mhist)))
    
        # Get the closest snapshot to the input redshift
        snap, z_snap = b.get_snap(zz,zmins[0],zmaxs[0],sim,env)

        # Get particle files
        files = b.get_particle_files(snap,sim,env)
        if (len(files)<1):
            print('WARNING (bahamascal): no subfind files at snap={}, {} '.format(snap,sim))
            continue

        for iff, ff in enumerate(files):
            f = h5py.File(ff, 'r')
            p0 = f['PartType0']

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
                #slim = 1./cosmo.tHubble(redshift) #1/Gyr ##here

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

            # Read M500 quantities
            fof = f['FOF']
            if (iff == 0):
                m500  = fof['Group_M_Crit500'][:]*1e10/h0      #Msun
                r500  = fof['Group_R_Crit500'][:]/h0           #cMpc
                cop_x = fof['GroupCentreOfPotential'][:,0]/h0  #cMpc
                cop_y = fof['GroupCentreOfPotential'][:,1]/h0  #cMpc
                cop_z = fof['GroupCentreOfPotential'][:,2]/h0  #cMpc
            else:
                m500  = np.append(m500,fof['Group_M_Crit500'][:]*1e10/h0)
                r500  = np.append(r500,fof['Group_R_Crit500'][:]/h0)
                cop_x = np.append(cop_x,fof['GroupCentreOfPotential'][:,0]/h0)
                cop_y = np.append(cop_y,fof['GroupCentreOfPotential'][:,1]/h0)
                cop_z = np.append(cop_z,fof['GroupCentreOfPotential'][:,2]/h0)
                
            # Read the stellar mass
            subhaloes = f['Subhalo']
            mass = subhaloes[massdef][:,itype]  #10^10Msun/h
            ind = np.where(mass > 0.)
            lm = np.log10(mass[ind]) + 10. - np.log10(h0) #Msun

            # Read the SFR 
            #sfr = subhaloes[nom_sfr][:] #Msun/h/yr  ###here
            
            f.close()
    
            # GSMF
            H, bins_edges = np.histogram(lm,bins=medges)
            gsmf[:] = gsmf[:] + H
    
            if (nvols != 'All'):
                if (iff>nvols): break
    
        print('Side of sim box = {:.2f} Mpc^3/h^3'.format(boxsize))
        if (boxsize<=0.):
            continue
        volume = np.power(boxsize/h0,3.) # In Mpc^3

        # fgas model
        df_part = pd.DataFrame(data=np.vstack([groupnum,subgroupnum,partmass,
                                               partx,party,partz]).T,
                               columns=["groupnum","subgroupnum","partmass",
                                        "partx","party","partz"])
        groupnum,subgroupnum,partmass,partx,party,partz=[[] for i in range(6)]

        df_fof = pd.DataFrame(data=np.vstack([m500,r500,cop_x,cop_y,cop_z]).T,
                              columns=["m500","r500","cop_x","cop_y","cop_z"])
        #m500,r500,cop_x,cop_y,cop_z=[[] for i in range(5)]
        #df_fof.index += 1
        #df_fof.index.names = ['groupnum']
        #df_fof.reset_index(inplace=True)
        #df_part.sort_values(by=['groupnum', 'subgroupnum'], inplace=True)
        #df_part.reset_index(inplace=True, drop=True)
        #merge = pd.merge(df_part, df_fof, on=['groupnum'])
        #merge['partx'] = merge.partx - merge.cop_x + boxsize / 2
        #merge['party'] = merge.party - merge.cop_y + boxsize / 2
        #merge['partz'] = merge.partz - merge.cop_z + boxsize / 2
        #merge.x.loc[merge.partx < 0] = merge.partx.loc[merge.x < 0] + boxsize
        #merge.y.loc[merge.party < 0] = merge.party.loc[merge.y < 0] + boxsize
        #merge.z.loc[merge.partz < 0] = merge.partz.loc[merge.z < 0] + boxsize
        #
        #merge = merge.loc[merge.m500 > 1e13]
        #merge['distance'] = (((boxsize / 2) - merge.partx) ** 2 + 
        #                     ((boxsize / 2) - merge.party) ** 2 + 
        #                     ((boxsize / 2) - merge.partz) ** 2) ** 0.5
        #merge['inside_r500'] = merge.distance <= merge.r500
        #merge = merge.loc[merge.inside_r500 == True]
        #groups = merge.groupby(['groupnum'], as_index=False)
        #gas_mass = groups.Mass.sum() ####here

        # fgas observations

        # GSMF model
        gsmf = gsmf/volume/dm  # In Msun/Mpc^3 
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
                    ecolor='grey',color='grey',mec='grey',
                    label='Baldry*2012, z>0.06')

        # Madau plot 
        fil_sfr = b.get_path2data(sim,env)+'sfr.txt'
        aexp, sfr = np.loadtxt(fil_sfr, usecols=(0,2),unpack=True)
        # SFR (total, end) [M0/yr] 
        zval = (1./aexp) - 1.

        medians = stats.perc_2arrays(zedges,zval,sfr/volume,0.5)
        ax5.plot(zhist,medians)
        print(medians) ###HERE

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
