
# Huifang Wang, INSERM U1106

from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division



import scipy.io
import os.path
import sys
from scipy import signal
#import pygraphviz as pgv

from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib import font_manager


#sys.path.append( os.path.expanduser('/Users/huifangwang/BHBS/NSeeg/Codes/Codes/') )
#import ns_tools as nt

def mln_plt_mat_graph_SEEG(savefig,fMat,mmRule,area_list,Rmat=[],figtitle=[]):
    
    nchan,nchan=np.shape(fMat)

    if Rmat==[]:
        flagM=0
        color_wl='DeepPink'
        color_sl='DarkGreen'
    else:
        flagM=1
    
    A=pgv.AGraph(directed=True,strict=True,splines='spline')
    A.layout(prog='dot')
    inode=1
    for iarea, nleed, isposition, ifposition in area_list:
        for ileed in np.arange(nleed):
            ileedpos=[(ifposition[0]-isposition[0])/float(nleed)*ileed+isposition[0],
                      (ifposition[1]-isposition[1])/float(nleed)*ileed+isposition[1]]
            print(ileedpos)
            A.add_node(inode,label=iarea+str(ileed+1)+'-',#color=Color_ch[inode],
                         pos=str(ileedpos[0])+','+str(ileedpos[1])+'!',width='0.5!',fontsize=12,fixedsize=True)
            inode=inode+1

    # add the edges for the nodes
    nodes=np.arange(0,nchan).astype(int)
    for inode in nodes:
        for jnode in nodes:
            if fMat[inode,jnode]==mmRule[0][0]:
                if flagM==1:
                    if Rmat[inode][jnode]<0.1:
                        color_sl='Navy'
                        stype_f='dashed'
                    else:
                        color_sl='DarkGreen'
                        stype_f='solid'
                    A.add_edge(inode+1,jnode+1,penwidth=2.0,label=Rmat[inode][jnode],fontcolor=color_sl,color=color_sl,style=stype_f)
                else:
                    A.add_edge(inode+1,jnode+1,penwidth=2.0,color=color_sl)
        
            elif fMat[inode,jnode]>=mmRule[1][0]:
                if flagM==1:
                    if Rmat[inode][jnode]<0.1:
                        color_wl='Navy'
                        stype_f='dashed'
                    else:
                        color_wl='DeepPink'
                        stype_f='solid'
                    A.add_edge(inode+1,jnode+1,penwidth=1.0,label=Rmat[inode][jnode],fontcolor=color_wl,color=color_wl,style=stype_f)
                else:
                    A.add_edge(inode+1,jnode+1,penwidth=1.0,color=color_wl)
            
            elif flagM==1 and Rmat[inode][jnode]>0.05:
                   A.add_edge(inode+1,jnode+1,penwidth=1.0,label=Rmat[inode][jnode],fontcolor='DarkViolet',color='DarkViolet',style='dashed',arrowType='box')
    
    
    A.draw(savefig)


def mln_plt_mat_graph(savefig,fMat,mmRule,Rmat=[],nodelabels=[],figtitle=[],subnetwork=5):
    
    nchan,nchan=np.shape(fMat)
    if np.size(subnetwork)==1:
        if np.mod(nchan,5)==0:
            subnetwork=(np.ones(nchan/5)*5).astype(int)
        else:
            print("Please specify the subnetwork!!")

    if Rmat==[]:
        flagM=0
        color_wl='DeepPink'
        color_sl='DarkGreen'
    else:
        flagM=1
    Ngroup=len(subnetwork)
    xRange = 500
    yRange = 500
    
    GroupRad=300
    withGroupRad=100
    
    noderadiuscircle=30
    if nodelabels==[]:
        flagnl=0
    else:
        flagnl=1
    
    ThetaG=np.linspace(0.,2.*np.pi,Ngroup+1)
    ThetaG=ThetaG[:-1]
    print(ThetaG)
    phoG=np.ones(Ngroup)*GroupRad
    [X_G,Y_G]=[np.cos(ThetaG)*phoG,np.sin(ThetaG)*phoG]
    print(subnetwork)
    A=pgv.AGraph(directed=True,strict=True,splines='spline')
    A.layout(prog='dot')
    inode=1
    for indG, ncsubnetwork in enumerate(subnetwork):
        iwithGroupRad=withGroupRad*ncsubnetwork/np.max(subnetwork)
        ThetaiG=np.linspace(0.,2.*np.pi,ncsubnetwork+1)
        ThetaiG=ThetaiG[:-1]
        RHO_iG=iwithGroupRad+noderadiuscircle
        for indig,ignode in enumerate(np.arange(ncsubnetwork)+inode):
            ileedpos=[np.cos(ThetaiG[indig])*RHO_iG+X_G[indG],
                      np.sin(ThetaiG[indig])*RHO_iG+Y_G[indG]]
                
            if flagnl==1:
               ilabel=str(nodelabels[inode-1])
            else:
               ilabel=str(inode)
               A.add_node(inode,label=ilabel,
                             pos=str(ileedpos[0])+','+str(ileedpos[1])+'!',width='0.5!',fontsize=18)
            inode=inode+1
            print(ileedpos)

    # add the edges for the nodes
    nodes=np.arange(0,nchan).astype(int)
    for inode in nodes:
        for jnode in nodes:
            if fMat[inode,jnode]==mmRule[0][0]:
                if flagM==1:
                    if Rmat[inode][jnode]<0.1:
                        color_sl='Navy'
                        stype_f='dashed'
                    else:
                        color_sl='DarkGreen'
                        stype_f='solid'
                    A.add_edge(inode+1,jnode+1,penwidth=2.0,label=Rmat[inode][jnode],fontcolor=color_sl,color=color_sl,style=stype_f)
                else:
                    A.add_edge(inode+1,jnode+1,penwidth=2.0,color=color_sl)
        
            elif fMat[inode,jnode]>=mmRule[1][0]:
                if flagM==1:
                    if Rmat[inode][jnode]<0.1:
                        color_wl='Navy'
                        stype_f='dashed'
                    else:
                        color_wl='DeepPink'
                        stype_f='solid'
                    A.add_edge(inode+1,jnode+1,penwidth=1.0,label=Rmat[inode][jnode],fontcolor=color_wl,color=color_wl,style=stype_f)
                else:
                    A.add_edge(inode+1,jnode+1,penwidth=1.0,color=color_wl)

            elif flagM==1 and Rmat[inode][jnode]>0.05:
                A.add_edge(inode+1,jnode+1,penwidth=1.0,label=Rmat[inode][jnode],fontcolor='DarkViolet',color='DarkViolet',style='dashed',arrowType='box')
        
        
    A.draw(savefig)


def mln_plt_str_graph(savefig,nchan,condi,nodelabels=[],figtitle=[],subnetwork=5):
    # this function to plot the structures
    if np.size(subnetwork)==1:
        if np.mod(nchan,5)==0:
            subnetwork=(np.ones(nchan/5)*5).astype(int)
        else:
            print("Please specify the subnetwork!!")
    color_sl='DarkGreen'
    Ngroup=len(subnetwork)
    scale=0.6
    xRange = 500*scale
    yRange = 500*scale
    
    GroupRad=300*scale
    withGroupRad=100*scale
    
    noderadiuscircle=50
    if nodelabels==[]:
        flagnl=0
    else:
        flagnl=1
    
    ThetaG=np.linspace(0.,2.*np.pi,Ngroup+1)
    ThetaG=ThetaG[:-1]
    print(ThetaG)
    phoG=np.ones(Ngroup)*GroupRad
    [X_G,Y_G]=[np.cos(ThetaG)*phoG,np.sin(ThetaG)*phoG]
    print(subnetwork)
    A=pgv.AGraph(directed=True,strict=True,splines='spline')
    A.layout(prog='dot')
    inode=1
    for indG, ncsubnetwork in enumerate(subnetwork):
        iwithGroupRad=withGroupRad*ncsubnetwork/np.max(subnetwork)
        ThetaiG=np.linspace(0.,2.*np.pi,ncsubnetwork+1)
        ThetaiG=ThetaiG[:-1]
        RHO_iG=iwithGroupRad+noderadiuscircle
        for indig,ignode in enumerate(np.arange(ncsubnetwork)+inode):
            ileedpos=[np.cos(ThetaiG[indig])*RHO_iG+X_G[indG],
                      np.sin(ThetaiG[indig])*RHO_iG+Y_G[indG]]
                
            if flagnl==1:
                ilabel=str(nodelabels[inode-1])
            else:
                ilabel=str(inode)
            A.add_node(inode,label=ilabel,
                       pos=str(ileedpos[0])+','+str(ileedpos[1])+'!',width='0.5!',fontsize=18)
            inode=inode+1
            print(ileedpos)

    # add the edges for the nodes
    for icondi in condi:
        A.add_edge(icondi[0],icondi[1],penwidth=2.0,color=color_sl)
    
    A.draw(savefig)

def mln_polar_bar(savefigname,tptnVector,ax=[]):
    if ax==[]:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

    Nbin,Ncolor = np.shape(tptnVector)
    theta = np.arange(0.0, 2*np.pi, 2*np.pi/Nbin)
    width = np.pi*2/Nbin
    Vcs=np.arange(0,1.01,0.1)
    VcsLabel=['CS:\n'+str(iVcs) for iVcs in Vcs]
    thetaTicks = theta + width/2
    if Ncolor==3:
        color3=['DarkGreen','DeepPink','Navy']
    bottom=0.0
    for i, icolor in enumerate(color3):
        radii = tptnVector[:,i]
        print(radii)
        bars = ax.bar(theta, radii, width=width, bottom=bottom)
        for r,bar in zip(radii, bars):
            bar.set_facecolor(icolor)
            bar.set_alpha(0.7)
        bottom=radii+bottom
    ax.xaxis.set_ticks(thetaTicks)
    ax.xaxis.set_ticklabels(VcsLabel,fontsize=20)
    ax._r_label_position._t = (0.0, 0.0)
    ax._r_label_position.invalidate()
    
    yticks_font = font_manager.FontProperties(family='sans-serif', style='italic', size=18, weight='bold',stretch='normal')
    for label in ax.get_yticklabels():
        label.set_fontproperties(yticks_font)
    plt.tick_params(axis='y',labelcolor='Orange')
    plt.savefig(savefigname)

def mln_polar_bar_statics(savefigname,tptnVector3d):
    '''
        for plot the statistical results: median and its percentage
        tptnVector3d is Nbin, Ncolor, Nsamples'''
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    
    tptnVectorm=np.median(tptnVector3d,axis=2)
    Nbin,Ncolor = np.shape(tptnVectorm)
    tptnVectorm75=np.percentile(tptnVector3d,75,axis=2)
    tptnVectorm25=np.percentile(tptnVector3d,25,axis=2)
    theta = np.arange(0.0, 2*np.pi, 2*np.pi/Nbin)
    width = np.pi*2/Nbin
    Vcs=np.arange(0,1.01,0.1)
    bias = width/5
    VcsLabel=['CS:\n'+str(iVcs) for iVcs in Vcs]
    thetaTicks = theta + width/2
    if Ncolor==3:
        color3=['DarkGreen','DeepPink','Navy']
    bottom=0.0
    # draw the radii
    for i, icolor in enumerate(color3):
        radii = tptnVectorm[:,i]
        bars = ax.bar(theta, radii, width=width, bottom=bottom)
        for r,bar in zip(radii, bars):
            bar.set_facecolor(icolor)
            bar.set_alpha(0.7)
        bottom=radii+bottom
    bottom=np.zeros(np.shape(radii))
    # draw the radii
    for i, icolor in enumerate(color3):
        radii = tptnVectorm[:,i]
        radii75=tptnVectorm75[:,i]
        radii25=tptnVectorm25[:,i]
        #opionts=ax.plot(theta+bias*(i+1)+width/5/2,radii+bottom,'o',color=icolor)
        pebars=ax.bar(theta+bias*(i+1),radii75,width=width/5,bottom=radii25+bottom)
            #for per,ipebar in zip(radii75-(radii25+bottom),pebars):
        for idpebar,ipebar in enumerate(pebars):
            ipebar.set_facecolor(icolor)
            ipebar.set_alpha(1)
            ipebar.set_edgecolor('none')
            if radii[idpebar]>0.0000000001:
                ax.plot(theta[idpebar]+bias*(i+1)+width/5/2,radii[idpebar]+bottom[idpebar],'o',color='red')

        bottom=radii+bottom

    ax.xaxis.set_ticks(thetaTicks)
    ax.xaxis.set_ticklabels(VcsLabel,fontsize=20)
    ax._r_label_position._t = (0.0, 0.0)
    ax._r_label_position.invalidate()
    ax.set_xlim([0,1])
    
    yticks_font = font_manager.FontProperties(family='sans-serif', style='italic', size=18, weight='bold',stretch='normal')
    for label in ax.get_yticklabels():
        label.set_fontproperties(yticks_font)
    plt.tick_params(axis='y',labelcolor='Orange')
    ax.set_ylim([0,1])
    plt.savefig(savefigname)


def mln_polar_bar_multi(savefigname,tptnVector,Cmv):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    
    Nbin,Ncolor = np.shape(tptnVector)
    theta = np.arange(0.0, 2*np.pi, 2*np.pi/Nbin)
    width = np.pi*2/Nbin
    Vcs=np.arange(0,1.01,0.1)
    VcsLabel=['CS:\n'+str(iVcs) for iVcs in Vcs]
    thetaTicks = theta + width/2
    if Ncolor==3:
        color3=['DarkGreen','DeepPink','Navy']
    else:
        color3=Cmv
    bottom=0.0
    for i, icolor in enumerate(color3):
        radii = tptnVector[:,-(i+1)]
        print(radii),icolor
        bars = ax.bar(theta, radii, width=width, bottom=bottom)
        for r,bar in zip(radii, bars):
            bar.set_facecolor(icolor)
            bar.set_alpha(0.7)
        bottom=radii+bottom
    ax.xaxis.set_ticks(thetaTicks)
    ax.xaxis.set_ticklabels(VcsLabel,fontsize=20)
    ax._r_label_position._t = (0.0, 0.0)
    ax._r_label_position.invalidate()
    
    yticks_font = font_manager.FontProperties(family='sans-serif', style='italic', size=18, weight='bold',stretch='normal')
    for label in ax.get_yticklabels():
        label.set_fontproperties(yticks_font)
    plt.tick_params(axis='y',labelcolor='Orange')
    plt.savefig(savefigname)


def mln_polar_bar_hf(savefigname,tptnVector,ax=[]):
    # this polar bar is used to highlight the false rate
    if ax==[]:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    
    Nbin,Ncolor = np.shape(tptnVector)
    theta = np.arange(0.0, 2*np.pi, 2*np.pi/Nbin)
    width = np.pi*2/Nbin
    Vcs=np.arange(0,1.01,0.1)
    VcsLabel=['CS:\n'+str(iVcs) for iVcs in Vcs]
    thetaTicks = theta + width/2
    if Ncolor==3:
        color3=['DarkGreen','DeepPink','Navy']
    bottom=0.0
    for i, icolor in enumerate(color3):
        radii = tptnVector[:,i]
        bars = ax.bar(theta[1:], radii[1:], width=width, bottom=bottom)
        
        for r,bar in zip(radii, bars):
            bar.set_facecolor(icolor)
            bar.set_edgecolor(icolor)
            bar.set_alpha(0.7)
        bottom=radii[1:]+bottom
    # for cs=0, Vcs=0
    bottom=0.0
    c3r=[color3[i] for i in np.arange(np.size(color3),0,-1)-1]
    for i,icolor in enumerate(c3r):
        radii=tptnVector[0,2-i]
        bars=ax.bar(theta[0],radii,width=width, bottom=bottom,facecolor=icolor,edgecolor=icolor,linewidth=2,alpha=0.7)
        
        #bars.set_facecolor(icolor)
        #bars.set_alpha(0.7)
        bottom=radii+bottom
    
    
    ax.xaxis.set_ticks(thetaTicks)
    ax.xaxis.set_ticklabels(VcsLabel,fontsize=20)
    ax._r_label_position._t = (0.0, 0.0)
    ax._r_label_position.invalidate()

    yticks_font = font_manager.FontProperties(family='sans-serif', style='italic', size=18, weight='bold',stretch='normal')
    for label in ax.get_yticklabels():
        label.set_fontproperties(yticks_font)
    plt.tick_params(axis='y',labelcolor='Orange')
    plt.savefig(savefigname)

def ns_plot_timeFrequencyS_6r2c_given_wavelet(gs1,titleax,yEEG,fs,yt,freqband,Go=7,stepfreqraw=1,ylabel='Raw',flagsp=True,flagrawlog =1,xtickstep=60):
    """ two colomns
        given seeg
        @huifang:April 20 2015, 3 row """
        
    ts=int(yt[0])
    tf=int(yt[-1])
    
    xticks_given=np.arange((int(yt[0])/xtickstep+1)*xtickstep,int(yt[-1]),xtickstep).astype(int)
    xticks_given=np.append(xticks_given,int(yt[0]))
    xticks_given=np.append(xticks_given,int(np.ceil(yt[-1])))
    
    #stepfreqraw=0.5
    
    
    #ggs1.update(left=0.05, right=0.95, wspace=0.1)
    
    xgrid=np.arange(ts,tf)
    ax=plt.subplot(gs1[0,0])
    ax.plot(yt,yEEG,color='DarkGreen',lw=2)
    ax.grid(True)
    ylim=ax.get_ylim()
    ax.vlines(xgrid,ylim[0],ylim[1],color='DarkGray',linestyles='--')
    ax.axis('tight')
    ax.set_xticks(xticks_given)
    ax.set_xticklabels(xticks_given)
    ax.set_title(titleax, color='DeepPink', fontsize=18)
    ax.set_ylabel(ylabel, fontsize=16)
    
    freqlistRaw=np.arange(freqband[0],freqband[1]+stepfreqraw,stepfreqraw)
    ppro=['raw','Flatting'];
    
    freqbandData=[('Delta',[0.5,3],0.5,'Crimson'),('Theta',[3,6],0.5,'DeepPink'),('gamma',[30,90],1,'Teal')]    
    
    
    for i, ifreqband in enumerate(freqbandData):
        #stepfreq=ifreqband[2]
        
        b, a = signal.butter(3, [2.0*ifreqband[1][0]/fs, ifreqband[1][1]*2.0/fs], 'bandpass')
        #print b,a
        iSEEGfilted= signal.filtfilt(b, a, yEEG)
        ax0=plt.subplot(gs1[i+1,0])
        ax0.plot(yt,iSEEGfilted)
        ax0.axis('tight')
        ax0.set_ylabel(ifreqband[0], fontsize=18)
        ax0.set_xticks([])
        px0=plt.subplot(gs1[i+1,1])
        iscalo=nt.ns_scalogram(iSEEGfilted,fs,freqlistRaw)
        showScalo=iscalo.mean(axis=0)
        px0.plot(freqlistRaw,showScalo,color=ifreqband[3],lw=2)
        sSlim=[showScalo.min()-showScalo.ptp()*0.05,showScalo.max()+showScalo.ptp()*0.05]
        px0.set_ylim([sSlim[0],sSlim[1]])
        px0.vlines(ifreqband[1],sSlim[0],sSlim[1],color='gray',linestyles='dashed')
        #par0.plot(freqlistRaw,iscalo.mean(axis=0),color=ifreqband[3],lw=1.5)
        px0.set_xlim([freqband[0],freqband[1]])
        px0.text(40,showScalo.ptp()/2+sSlim[0],str(ifreqband[1][0])+'-'+str(ifreqband[1][1])+' Hz',fontsize=14)
        px0.set_yticks([])
    
    Nabove=len(freqbandData)+1
    
    
    for indppro, ippro in enumerate(ppro):
        
        
        if ippro=='raw':
            iyEEG=yEEG
            wt=nt.ns_wt(yEEG,fs,freqlistRaw,xi=Go)
            ititlegs='Original'
            wtPowerRaw=np.abs(wt)**2
            if flagrawlog == 1:
                wtPowerRaw[wtPowerRaw<10**-3]=10**-3
                showtm=10*np.log10(wtPowerRaw)
            else:
                showtm=wtPowerRaw
            colorf = 'Orange'
        elif ippro=='diff':
            iyEEG=np.diff(yEEG)
            iyEEG=np.insert(iyEEG,0,0)
            wt=nt.ns_wt(iyEEG,fs,freqlistRaw,xi=Go)
            ititlegs='diff'
            wtPowerRaw=np.abs(wt)**2
            #wtPowerRaw=np.abs(wt)
            showtm=wtPowerRaw
            
        elif ippro=='Flatting':
            iyEEG=yEEG
            wt=nt.ns_wt(yEEG,fs,freqlistRaw,xi=Go)
            wt=wt*freqlistRaw
            ititlegs='Spectral flattening'
            wtPowerRaw=np.abs(wt)**2
            #wtPowerRaw=np.abs(wt)
            showtm=wtPowerRaw
            colorf = 'Navy'
        #print(indppro+Nabove,ippro)
        ax1=plt.subplot(gs1[indppro+Nabove,0])
#        X= yt
#        Y=freqlistRaw
#        #print yt, freqlistRaw
#        X, Y  = np.meshgrid(yt,freqlistRaw)
#        #print(np.shape(X), np.shape(Y), np.shape(showtm.T))
#        #ax1.contourf(X,Y,showtm.T,100,cmap='Spectral_r')
#        ax1.pcolorfast(X,Y,showtm.T[:-1,:-1],cmap='Spectral_r')
#        
#        ax1.set_yscale('log')
#        cxlim=ax1.get_xlim()
#        #ax1.set_ylim([1,30])
#        xvlabels=[3,6,30,90]
#        xvlabelsm=[3,6,30,90]
#        ax1.hlines(xvlabels,cxlim[0],cxlim[1],color='gray',linestyles='dashed')
#        ax1.set_yticks=(np.log(xvlabelsm))
#        ax1.set_yticks=xvlabelsm
#        ax1.set_yticklabels=(xvlabelsm)
#        ax1.set_xticks(xticks_given)
#        ax1.set_xticklabels(xticks_given)
#        ax1.yaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
#        ax1.yaxis.set_minor_locator(plt.FixedLocator(xvlabelsm))
        xvlabels2=[1,3,6,15,30]
        ax1.imshow(showtm.T, cmap='Spectral_r',extent=[yt.min(),yt.max(), freqlistRaw.max(),freqlistRaw.min()],aspect='auto')
        ax1.invert_yaxis()        
        ax1.set_title(ititlegs,color='Black', fontsize=12)

        if flagsp==True:
            par1=plt.subplot(gs1[0,1])
            mshowtm=showtm.mean(axis=0)
            par1.plot(freqlistRaw,mshowtm/np.max(mshowtm),color=colorf,lw=3)
        
            par1.set_yticks([])
            #mshowtm=showtm.mean(axis=0)
            par0=plt.subplot(gs1[indppro+Nabove,1])
            par0.plot(freqlistRaw,mshowtm,color=colorf,lw=3)
            par0.set_xscale('log')
            #plt.savefig(data_path+'TimeF_20'+Pname+data_filename+channelname+'ts_'+str(ts)+'Pattern.png')
            par0.set_xticks([])
            par0.set_yticks([])
            par0.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
            par0.xaxis.set_minor_locator(plt.FixedLocator(xvlabels2))
            cylim=par0.set_ylim()
            par0.vlines(xvlabels2,cylim[0],cylim[1],color='gray',linestyles='dashed')
        
        ax1.set_ylabel('Frequency',fontsize=16)
    ax1.set_xlabel('Time in seconds',fontsize=16)

    
    return gs1
    
    
def ns_plot_timeFrequencyS_3r2c_given_wavelet(gs1,titleax,yEEG,fs,yt,freqband,stepfreqraw=1,ylabel='Raw',flagsp=True,xtickstep=60,Nrow=3):
    """ two colomns
        given seeg
        @huifang:April 20 2015, 3 row """
        
    ts=int(yt[0])
    tf=int(yt[-1])
    
    xticks_given=np.arange((int(yt[0])/xtickstep+1)*xtickstep,int(yt[-1]),xtickstep).astype(int)
    xticks_given=np.append(xticks_given,int(yt[0]))
    xticks_given=np.append(xticks_given,int(np.ceil(yt[-1])))
    
    #stepfreqraw=0.5
    
    
    #ggs1.update(left=0.05, right=0.95, wspace=0.1)
    
    xgrid=np.arange(ts,tf)
    ax=plt.subplot(gs1[0,0])
    ax.plot(yt,yEEG,color='DarkGreen',lw=2)
    ax.grid(True)
    ylim=ax.get_ylim()
    ax.vlines(xgrid,ylim[0],ylim[1],color='DarkGray',linestyles='--')
    ax.axis('tight')
    ax.set_xticks(xticks_given)
    ax.set_xticklabels(xticks_given)
    ax.set_title(titleax, color='DeepPink', fontsize=18)
    ax.set_ylabel(ylabel, fontsize=16)
    
    freqlistRaw=np.arange(freqband[0],freqband[1]+stepfreqraw,stepfreqraw)
    ppro=['raw','Flatting'];
    for indppro, ippro in enumerate(ppro):
        
        
        if ippro=='raw':
            iyEEG=yEEG
            wt=nt.ns_wt(yEEG,fs,freqlistRaw,xi=7)
            ititlegs='Original'
            wtPowerRaw=np.abs(wt)**2
            wtPowerRaw[wtPowerRaw<10**-3]=10**-3
            showtm=10*np.log10(wtPowerRaw)
            colorf = 'Orange'
        elif ippro=='diff':
            iyEEG=np.diff(yEEG)
            iyEEG=np.insert(iyEEG,0,0)
            wt=nt.ns_wt(iyEEG,fs,freqlistRaw,xi=7)
            ititlegs='diff'
            wtPowerRaw=np.abs(wt)**2
            #wtPowerRaw=np.abs(wt)
            showtm=wtPowerRaw
            
        elif ippro=='Flatting':
            iyEEG=yEEG
            wt=nt.ns_wt(yEEG,fs,freqlistRaw,xi=7)
            wt=wt*freqlistRaw
            ititlegs='Spectral flattening'
            wtPowerRaw=np.abs(wt)**2
            #wtPowerRaw=np.abs(wt)
            showtm=wtPowerRaw
            colorf = 'Navy'
        
        ax1=plt.subplot(gs1[indppro+1,0])
        #X= yt
        #Y=freqlistRaw
        #print yt, freqlistRaw
        #X, Y  = np.meshgrid(yt,freqlistRaw)
        #print(np.shape(X), np.shape(Y), np.shape(showtm.T))
        #ax1.contourf(X,Y,showtm.T,100,cmap='Spectral_r')
        #ax1.pcolorfast(X,Y,showtm.T[:-1,:-1],cmap='Spectral_r')
        ax1.imshow(showtm.T, cmap='Spectral_r',extent=[yt.min(),yt.max(), freqlistRaw.max(),freqlistRaw.min()],aspect='auto')
        ax1.invert_yaxis()
        #ax1.set_yscale('log')
        #cxlim=ax1.get_xlim()
        #ax1.set_ylim([2,30])
        #xvlabels=[4,8,15,20,]
        #xvlabelsm=[4,8,15,20]
    #ax1.hlines(xvlabels,cxlim[0],cxlim[1],color='gray',linestyles='dashed')
        #ax1.set_yticks=(np.log(xvlabelsm))
        #ax1.set_yticks=xvlabelsm
        #ax1.set_yticklabels=(xvlabelsm)
        #ax1.set_xticks(xticks_given)
        #ax1.set_xticklabels(xticks_given)
        #ax1.yaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
        #ax1.yaxis.set_minor_locator(plt.FixedLocator(xvlabelsm))
        xvlabels2=[1,4,8,15]
        ax1.set_title(ititlegs,color='Black', fontsize=16)

        if flagsp==True:
            par1=plt.subplot(gs1[0,1])
            mshowtm=showtm.mean(axis=0)
            par1.plot(freqlistRaw,mshowtm/np.max(mshowtm),color=colorf,lw=3)
        
            par1.set_yticks([])
            #mshowtm=showtm.mean(axis=0)
            par0=plt.subplot(gs1[indppro+1,1])
            par0.plot(freqlistRaw,mshowtm,color=colorf,lw=3)
            par0.set_xscale('log')
            #plt.savefig(data_path+'TimeF_20'+Pname+data_filename+channelname+'ts_'+str(ts)+'Pattern.png')
            par0.set_xticks([])
            par0.set_yticks([])
            par0.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
            par0.xaxis.set_minor_locator(plt.FixedLocator(xvlabels2))
            cylim=par0.set_ylim()
            par0.vlines(xvlabels2,cylim[0],cylim[1],color='gray',linestyles='dashed')
        
        ax1.set_ylabel('Frequency',fontsize=16)
    ax1.set_xlabel('Time in seconds',fontsize=16)

    
    return gs1
    
def ns_chosen_channels_fit_LFP(ax,ch_EEG,Chosen_ch,Color_ch,times,Normalized=True,scale=3.,steptime=100):
    ''' plot time series with out top and right axis'''

    Nchanv=np.arange(len(scale))
    for i, y, iscale in zip(Nchanv,ch_EEG[:, :],scale):
        
        if Normalized==True:
            plt.plot(times,(y - y.mean())/y.ptp()*iscale + i,color=Color_ch[i])
        else:
            plt.plot(times,(y - y.mean())/ch_EEG.ptp()*iscale + i,color=Color_ch[i])

#xticks_given=np.arange(0,int(max(times)),100)
    #itimes=[x[0] for x in remarks_times_vedio]
    xticks_given=[i for i in np.arange(min(times),max(times),steptime).astype(int)]
    
    #xticks_given=np.append(xticks_given,array(itimes).astype(int))
    xticks_given=np.append(xticks_given,int(max(times)))
    
    ax.set_yticklabels(Chosen_ch,fontsize=20)
    ax.set_yticks(range(0,len(Chosen_ch)))
    ax.set_xlim([min(times),max(times)])
    ax.set_ylim([-1,len(Chosen_ch)])
    ax.set_xticks(xticks_given)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax.set_xticklabels(xticks_given,fontsize=18)
    #draw_axvlines(ax,remarks_times_vedio)
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.gca().invert_yaxis()

    

