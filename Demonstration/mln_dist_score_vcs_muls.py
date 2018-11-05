
# coding: utf-8

"""For varying CS, which correponds to the distribution of scores of MULAN results
This version for multiple structures
"""
# In[1]:

import scipy.io
import os.path
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import sys

import sys, os.path
#sys.path.append( os.path.expanduser('/Users/huifangwang/BHBS/BHBS/Programs') )
#import ns_tools as nt


Gs=[(0,0,1),(1,0,0)]
th1s=50
th2s=60
th1f=60
th2f=80


while len(sys.argv) > 1:
    option = sys.argv[1];                                              del sys.argv[1]
    if   option == '-dir':   basicdir = sys.argv[1];  del sys.argv[1]
    #elif option == '-filename':   data_filename = sys.argv[1];  del sys.argv[1]
    #elif option == '-is':  istru = int(sys.argv[1]);                    del sys.argv[1]
    
    elif option == '-th1s': th1s=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th1f': th1f=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th2s': th2s=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th2f': th2f=int(sys.argv[1]); del sys.argv[1]
    
    #elif option == '-ts': ts=int(sys.argv[1]); del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

#basicdir='/Volumes/HWMac1/MULAN2/SPM/VaryingCS/nmmN10L15Vcs/'
#basicdir='/Volumes/MAC_HF/MULAN/SPM/nmmN10L15Vcs/'
basicdir='/Volumes/Lioncub/Users/huifangwang/MULAN2/SPM/nmmN10L15Vcs/'

savedir=basicdir+'/disScoreFigs/'


if not os.path.exists(savedir):
    os.makedirs(savedir)



def mln_computer_idistScore(idistScore,iGMat,csV,Rmat):
    nchan,nchan,nbt=np.shape(iGMat)
    indices=np.arange(nchan)
    mlnMat=np.empty((nchan,nchan))
    for i, ichan in enumerate(indices):
        for j, jchan in enumerate(indices):
            if j==i:
               continue 
            pij=iGMat[ichan,jchan,:]
            csR=Rmat[ichan,jchan]
            if csR<0.1:
                idistScore[0][1]=np.append(idistScore[0][1],pij)
            else: 
                idistScore[csV.index(csR)][1]=np.append(idistScore[csV.index(csR)][1],pij)
    return idistScore




isV=np.arange(1,11,1)

dirname=basicdir+'ToutResults/'
csV=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
idistScore=[[i,[]] for i in csV]

vth1=np.arange(th1s,th1f+1,5)
vth2=np.arange(th2s,th2f+1,5)

for istru in isV:

    basicnom='nmmN10L15CS10100S'+str(istru)+'N6100'
    toutprenom='Tout_'+basicnom+'.mat'
    tout=scipy.io.loadmat(dirname+toutprenom)
    Rmat=tout['Connectivity']

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
                idistScore=mln_computer_idistScore(idistScore,iGMat,csV,Rmat)


# In[33]:

fig=plt.figure(figsize=[10,4])
gs1 = gridspec.GridSpec(2, 7)
gs1.update(left=0.05, right=0.95, top=0.9,bottom=0.05,wspace=0.2,hspace=0.3)
widths=1./50
colorhist='b'
lineWidth=2.5
coloredge='red'
nbins=20



for indcsV,ics in enumerate(csV):
    pij=idistScore[indcsV][1]
    hist, bins = np.histogram(pij, bins=nbins, range=(0.,1.), density=False)
    ax0=plt.subplot(gs1[indcsV])
    ax0.bar(bins[:-1], hist, widths,
            color=colorhist, edgecolor=coloredge,linewidth=lineWidth)
    ax0.set_xlim(0-widths,1.+widths)
    ax0.set_yticklabels([])
    ax0.set_title('CS='+str(ics),fontsize=14)
    if indcsV==0:
        allhist=hist
        perhist=hist/float(np.sum(hist))
    else:
        allhist=np.r_['0,2',allhist,hist]
        perhist=np.r_['0,2',perhist,hist/float(np.sum(hist))]
        


# In[85]:


myCcmap='PiYG'
spe_cmap=plt.get_cmap(myCcmap)

idScorelist=[idistScore[i][1] for i in xrange(len(idistScore))]
idScoremean=[idistScore[i][1].mean() for i in xrange(len(idistScore))]
idScoremedian=[np.median(idistScore[i][1]) for i in xrange(len(idistScore))]

# In[181]:

fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])
sbins=bins[:-1]
bottomhist=np.zeros(np.size(perhist,axis=0))
width=0.05
for indbins,ibins in enumerate(sbins):
    iperhist=perhist[:,indbins]
    plt.bar(csV,iperhist,width=width,
            bottom=bottomhist,color=spe_cmap(ibins),edgecolor=spe_cmap(ibins),align='center')
    bottomhist=bottomhist+iperhist
for indcsV,icsV in enumerate(csV): 
    plt.plot([icsV-width/2,icsV+width/2],[idScoremedian[indcsV],idScoremedian[indcsV]],
             'red',linewidth=3)
plt.plot(csV,idScoremean,'b*',markersize=10)
ax.set_ylabel('percent of distributed scores;\n scores')
ax.set_xlim(0-width,1+width)
ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
cmap = cm.get_cmap(myCcmap, nbins)
ax1 = fig.add_axes([0.85, 0.1, 0.1, 0.8])
ax1.set_axis_off()
sm = cm.ScalarMappable(cmap=cmap)
sm._A = []
#cbar = fig.colorbar(sm, aspect=42)
ax.set_xlabel('CS')
ax.set_ylim([0-0.01,1+0.01])
sfilename='S_all'+str(istru)+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)+'.eps'
ax.set_title('S_all'+str(istru)+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)+'\n'+str(Gs), 
             fontsize=14, color='Blue')
plt.savefig(savedir+sfilename)