import scipy.io
import os.path
import sys
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm



def mln_plt_dist_score(ax,bins,perhist,medscore,meanscore,spe_cmap):
    sbins=bins[:-1]
    icsV=1
    bottomhist=0
    width=0.5
    
    for indbins,ibins in enumerate(sbins):
        iperhist=perhist[indbins]
        plt.bar(icsV,iperhist,width=width,
            bottom=bottomhist,color=spe_cmap(ibins),edgecolor=spe_cmap(ibins),align='center')
        bottomhist=bottomhist+iperhist

    plt.plot([icsV-width/2,icsV+width/2],[medscore,medscore],
             'red',linewidth=3)
    plt.plot(icsV,meanscore,'b*')


def mln_plt_sta_simulated_dscore(allGMat,fig,Rmat,mmRule,spe_cmap,pltflag=False):
    nchan,nchan,nbt=np.shape(allGMat)
    mdMat=np.zeros([nchan,nchan])
    mnMat=np.zeros([nchan,nchan])
    fMat=np.zeros([nchan,nchan])
    if pltflag==1:
        gs1 = gridspec.GridSpec(nchan, nchan)
        gs1.update(left=0.05, right=0.95, top=0.9,bottom=0.05,wspace=0.1)
    indices=np.arange(1,nchan+1)
    normed_value = 2
    for i, ichan in enumerate(indices):
        for j, jchan in enumerate(indices):
            nbt=20
            if pltflag==1:
                ax0=plt.subplot(gs1[i,j])
            if j==i:
                if pltflag==1:
                    ax0.set_xticks([])
                    ax0.set_yticks([])
                    ax0.patch.set_facecolor('LightGrey')
                    plt.setp(ax0.spines.values(), color='none')
                    #plt.axis('off')
                    if j==0:
                        ax0.set_ylabel(str(ichan),fontsize=12)
                    if j==nchan-1:
                        ax0.set_xlabel(str(ichan),fontsize=12)
                continue
    
            pij=allGMat[i,j,:]
            mdMat[i,j]=medscore=np.median(pij)
            mnMat[i,j]=meanscore=np.mean(pij)
            if mmRule[0][1]<=medscore<=1 and mmRule[0][2]<=meanscore<=1:
                fMat[i,j]=mmRule[0][0]
            elif mmRule[1][1]<=medscore<=1 and mmRule[1][2]<=meanscore<=1:
                fMat[i,j]=mmRule[1][0]
            else:
                fMat[i,j]=0
                            
            
            if pltflag==1:
                if j==0:
                    ax0.set_ylabel(str(ichan),fontsize=14)
                if i==nchan-1:
                    ax0.set_xlabel(str(jchan),fontsize=14)
                
                hist, bins = np.histogram(pij, bins=20, range=(0.,1.), density=False)
                perhist=hist/float(np.sum(hist))

                
                if Rmat[i,j]>0.05:
                    
                    ax0.text(0.25,0.5,'('+str(Rmat[i,j])+')',
                             horizontalalignment='center',
                             verticalalignment='center',color='darkgreen')
                
                
                widths=0.5
                mln_plt_dist_score(ax0,bins,perhist,mdMat[i,j],mnMat[i,j],spe_cmap)
                ax0.set_xlim(0-widths,1.+widths)
                ax0.set_ylim(0-widths,1+widths)
                
                
                ax0.set_xticks([])
                ax0.set_yticks([])
                ax0.patch.set_facecolor('LightGrey')
                plt.setp(ax0.spines.values(), color='g')
                ax0.spines['top'].set_color('none')
                ax0.spines['right'].set_color('none')
    return mdMat,mnMat,fMat

def mln_cal_mat_mltp(basicdir,basicnom,Gs,th1s,th1f,istru,th2s,th2f,mmRule,pltflag=0):



    dirname=basicdir+'ToutResults/'
    savedir=basicdir+'/stamlnvcsfigure/'




    if not os.path.exists(savedir):
        os.makedirs(savedir)


    toutprenom='Tout_'+basicnom+'.mat'
    tout=scipy.io.loadmat(dirname+toutprenom)
    try:
        Rmat=tout['Connectivity']
    except:
        Rmat=[]


    allGMat=[]
    vth1=np.arange(th1s,th1f+1,5)
    vth2=np.arange(th2s,th2f+1,5)
    for th1 in vth1:
        for th2 in vth2:
            
            thnom='th12'+str(th1)+str(th2)
            staprenom='Sta_'+basicnom+thnom+'.mat'
            mln=scipy.io.loadmat(dirname+staprenom)
            iMULAN=mln['iMULAN']
            iMat=iMULAN['Mat']
            iMULANmat=iMat[0][0]
            for igs in Gs:
                
                iG1=igs[0];iG3=igs[1];iG4=igs[2]
                iGMat=iMULANmat[:,:,iG1,iG3,iG4,:]
                if allGMat==[]:
                    allGMat=iGMat
                else:
                    allGMat=np.concatenate((allGMat,iGMat),axis=2)

    if pltflag==1:
        myCcmap='PiYG'
        spe_cmap=plt.get_cmap(myCcmap)
        sfilename='S_'+basicnom+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)
        fig=plt.figure(figsize=(10,10),facecolor='LightGrey')
        fig.suptitle(sfilename,y=0.95,fontsize=14,
                     fontweight='bold',color='DarkGreen')
        mdMat,mnMat,fMat=mln_plt_sta_simulated_dscore(allGMat,fig,Rmat,mmRule,spe_cmap,pltflag=1)
        plt.savefig(savedir+sfilename+'.eps',format='eps',dpi=300)
        plt.close()
    else:
        fig=[]
        spe_cmap=[]
        mdMat,mnMat,fMat=mln_plt_sta_simulated_dscore(allGMat,fig,Rmat,mmRule,spe_cmap,pltflag=0)
        return mdMat, mnMat, fMat,Rmat

def mln_find_subgroup(inode,nodelabels):
    for lind,ilabels in enumerate(nodelabels):
        a=ilabels==inode
        if a[a.nonzero()]:
            return lind

def mln_plt_mat_mltp(basicdir,savedir,basicnom,Gs,th1s,th1f,istru,th2s,th2f,mmRule,nodelabels=[],pltflag=0,forceRMat=1):
    mdMat, mnMat, fMat,Rmat=mln_cal_mat_mltp(basicdir,basicnom,Gs,th1s,th1f,istru,th2s,th2f,mmRule,pltflag)
    if forceRMat==1:
        Rmat=[]
    sfilename='S_'+basicnom+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)
    savefig=savedir+'Graph_'+sfilename+'.eps'
    md.mln_plt_mat_graph(savefig,fMat,mmRule,Rmat,nodelabels,figtitle=sfilename)

# if forceRMat=0, we also output the statistical results
    if forceRMat==0:
        savefig=savedir+'TPTN_'+sfilename+'.pdf'
        tptnVectorR, tptnVector, pnVector=mln_tptn(fMat,Rmat)
        md.mln_polar_bar_hf(savefig,tptnVectorR)

def mln_tptn(fMat,Rmat,thetaw=1,thetas=2):
# reture true posivite true negative for connection strength [0,1]
# Rmat is the ground truth
    csV=np.arange(0,1.01,0.1)
    tptnVector=np.zeros(shape=(len(csV),3))
    pnVector=np.zeros(shape=(len(csV),1))
    nchan,nchan=np.shape(fMat)
    indices=np.arange(nchan).astype(int)
    for i in indices:
        for j in indices:
            if i==j:
                continue
            indexall=np.where(csV>=Rmat[i][j]-0.000001)
            indexcs=indexall[0]
            indexcs=indexcs[0]
            #print Rmat[i][j],indexcs
            
            if fMat[i][j]>=thetas:
                tptnVector[indexcs,0]=tptnVector[indexcs,0]+1
            elif fMat[i][j]>=thetaw:
                tptnVector[indexcs,1]=tptnVector[indexcs,1]+1
            else:
                tptnVector[indexcs,2]=tptnVector[indexcs,2]+1
            pnVector[indexcs]=pnVector[indexcs]+1
    
    tptnVectorR=tptnVector/pnVector
    #print(tptnVector,pnVector,tptnVectorR)
    #raw_input()
    return tptnVectorR, tptnVector, pnVector

def mln_tptn_HC(fMat,Rmat,thetaw=1,thetas=2):
    # reture true posivite true negative for connection strength [0,1] for the human connectome network
    # Rmat is the ground truth
    csV=np.arange(0,1.01,0.1)
    tptnVector=np.zeros(shape=(len(csV),3))
    pnVector=np.zeros(shape=(len(csV),1))
    nchan,nchan=np.shape(fMat)
    indices=np.arange(nchan).astype(int)
    for i in indices:
        for j in indices:
            if i==j:
                continue
            indexall=np.where(csV<=Rmat[i][j])
            indexcs=indexall[0]
            indexcs=indexcs[-1]
            if csV[indexcs]<0.001 and Rmat[i][j]>0:
                continue
            
            if fMat[i][j]>=thetas:
                tptnVector[indexcs,0]=tptnVector[indexcs,0]+1
            elif fMat[i][j]>=thetaw:
                tptnVector[indexcs,1]=tptnVector[indexcs,1]+1
            else:
                tptnVector[indexcs,2]=tptnVector[indexcs,2]+1
            pnVector[indexcs]=pnVector[indexcs]+1
    
    tptnVectorR=tptnVector/pnVector
    #print(tptnVector,pnVector,tptnVectorR)
    #raw_input()
    return tptnVectorR, tptnVector, pnVector

def mln_tptn_multi(fMat,Rmat,Th):
    # reture true posivite true negative for connection strength [0,1]
    # Rmat is the ground truth
    import bisect
    csV=np.arange(0,1.01,0.1)
    
    tptnVector=np.zeros(shape=(len(csV),len(Th)+1))
    pnVector=np.zeros(shape=(len(csV),1))
    nchan,nchan=np.shape(fMat)
    indices=np.arange(nchan).astype(int)
    for i in indices:
        for j in indices:
            if i==j:
                continue
            indexall=np.where(csV<=Rmat[i][j])
            indexcs=indexall[0]
            indexcs=indexcs[-1]
            if csV[indexcs]<0.001 and Rmat[i][j]>0:
                continue
            #print Rmat[i][j],indexcs
            indc=bisect.bisect(Th, fMat[i][j])
            tptnVector[indexcs,indc]=tptnVector[indexcs,indc]+1

            pnVector[indexcs]=pnVector[indexcs]+1
    
    tptnVectorR=tptnVector/pnVector
    #print(tptnVector,pnVector,tptnVectorR)
    #raw_input()
    return tptnVectorR, tptnVector, pnVector


if __name__ == '__main__':

    Gs=[(0,0,1),(1,0,0)]
    th1s=50
    th2s=60
    th1f=60
    th2f=80
    basicdir='/Volumes/Lioncub/Users/huifangwang/MULAN2/SPM/nmmN10L15Vcs/'
    
    pltflag=0
    istru=10

    mmRule=[(2,[0.85],[0.85]),(1,[0.3],[0.3])] #(mlnValue,medianScore, meanScore)
    while len(sys.argv) > 1:
        option = sys.argv[1];                             del sys.argv[1]
        if   option == '-dir':   basicdir = sys.argv[1];  del sys.argv[1]
        if   option == '-nom':   basicnom = sys.argv[1];  del sys.argv[1]
        #elif option == '-filename':   data_filename = sys.argv[1];  del sys.argv[1]
        elif option == '-is':  istru = int(sys.argv[1]);                    del sys.argv[1]
        
        elif option == '-th1s': th1s=int(sys.argv[1]); del sys.argv[1]
        elif option == '-th1f': th1f=int(sys.argv[1]); del sys.argv[1]
        elif option == '-th2s': th2s=int(sys.argv[1]); del sys.argv[1]
        elif option == '-th2f': th2f=int(sys.argv[1]); del sys.argv[1]
        elif option == '-pltflag': pltflag=bool(sys.argv[1]); del sys.argv[1]
        
        #elif option == '-ts': ts=int(sys.argv[1]); del sys.argv[1]
        else:
            print('Options invalides :',option,'->',sys.argv[0])

    basicnom='nmmN10L15CS10100S'+str(istru)+'N6100'
    savedir=basicdir+'/MlnResults/'
    
    if not os.path.exists(savedir):
        os.makedirs(savedir)


    mdMat, mnMat, fMat,Rmat=mln_cal_mat_mltp(basicdir,basicnom,Gs,th1s,th1f,istru,th2s,th2f,mmRule,pltflag)
    sfilename='S_'+basicnom+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)+'_Rmat_'+str(mmRule[1][1])
    savefig=savedir+'Graph_'+sfilename+'.eps'
#md.mln_plt_mat_graph(savefig,fMat,mmRule,Rmat)

    savefig=savedir+'TPTN_'+sfilename+'.pdf'
    tptnVectorR, tptnVector, pnVector=mln_tptn(fMat,Rmat)
    md.mln_polar_bar_hf(savefig,tptnVectorR)