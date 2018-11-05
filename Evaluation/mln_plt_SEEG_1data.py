import scipy.io
import os.path
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
sys.path.append(os.path.expanduser('/Users/huifangwang/MULANII/programs/Inference/'))
import mln_cal_mat_mltp as mi

sys.path.append(os.path.expanduser('/Users/huifangwang/MULANII/programs/Demonstrate/'))
import mln_plt_demonstration as md

sys.path.append(os.path.expanduser('/Users/huifangwang/MULANII/programs/basicTools/'))
import mln_tools as mt

import itertools

Gs=[(0,0,1),(1,0,0)]
th1s=50
th2s=60
th1f=60
th2f=70
basicdir='/Users/huifangwang/MULAN2/SEEG/MotorNetwork/GUERFI/'

pltflag=0
istru=1
forceRMat=0
while len(sys.argv) > 1:
    option = sys.argv[1];                             del sys.argv[1]
    if   option == '-dir':   basicdir = sys.argv[1];  del sys.argv[1]
    if   option == '-nom':   basicnom = sys.argv[1];  del sys.argv[1]
    #elif option == '-filename':   data_filename = sys.argv[1];  del sys.argv[1]
    elif option == '-ts':  ts = int(sys.argv[1]);                    del sys.argv[1]
    
    elif option == '-th1s': th1s=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th1f': th1f=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th2s': th2s=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th2f': th2f=int(sys.argv[1]); del sys.argv[1]
    elif option == '-pltflag': pltflag=bool(sys.argv[1]); del sys.argv[1]
    
    #elif option == '-ts': ts=int(sys.argv[1]); del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]



mmRule=[(2,[0.95],[0.95]),(1,[0.92],[0.92])]
savedir=basicdir+'/MlnResults/'+'S'+str(istru)+'/'+'Th1_'+str(th1s)+str(th1f)+'_Th2_'+str(th2s)+str(th2f)+'mmRule'+str(mmRule[1][1])+'/'

if not os.path.exists(savedir):
    os.makedirs(savedir)

# for N30
 #(mlnValue,medianScore, meanScore)
data_path=basicdir

ts=600
postfix='ds5'
Pname='GUERFI'
data_filename='130830B-CEX_0000'
#basicnom='GUERFI130830B-CEX_0000SACCOR35ts'+str(ts)+'ds5'
prenomfile=Pname+data_filename

basicnom=Pname+data_filename+'SACCOR35ts'+str(ts)+postfix
mdMat, mnMat, fMat,Rmat=mi.mln_cal_mat_mltp(basicdir,basicnom,Gs,th1s,th1f,istru,th2s,th2f,mmRule,pltflag)
#print(fMat)
sfilename='S_'+basicnom+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)
savefig=savedir+'Graph_'+sfilename+'.eps'
area_list=[('SA',9,[10,800],[600,800]),('CC',12,[10,700],[700,700]),('OR',14,[10,500],[800,700])]
md.mln_plt_mat_graph_SEEG(savefig,fMat,mmRule,area_list,Rmat=[],figtitle=[])

# for N20*3
Chosen_area=[('SA',10,1),('CC', 13,2),('OR',15,3)]
mmRule=[(2,[0.95],[0.95]),(1,[0.9],[0.9])]

signalgroupN = [iarea_list[1] for iarea_list in area_list]
inodelabel=0
signalgroup=[]
for isingalN in signalgroupN:
    ig=(inodelabel+1,inodelabel+isingalN)
    inodelabel=inodelabel+isingalN
    signalgroup.append(ig)

import itertools
ilist=itertools.combinations([0,1,2],2)

seticom=set(ilist)
nchan,nchan=np.shape(fMat)
orginallabels=np.arange(1,nchan+1).astype(int)

N20fMat=[]
for ilistCa in seticom:
    iChosen_area=[Chosen_area[karea] for karea in ilistCa]
    iarea_list=[area_list[karea] for karea in ilistCa]
    nodelabels=np.concatenate([np.arange(signalgroup[ilistCa[0]][0],signalgroup[ilistCa[0]][1]+1).astype(int),np.arange(signalgroup[ilistCa[1]][0],signalgroup[ilistCa[1]][1]+1).astype(int)])
    
    ifilename=mt.mln_SEEG_extact_filename(ts,prenomfile,postfix,iChosen_area)
    imdMat, imnMat, ifMat,iRmat=mi.mln_cal_mat_mltp(basicdir,ifilename,Gs,th1s,th1f,istru,th2s,th2f,mmRule,pltflag=0)

#savefig=savedir+'Graph_'+ifilename+'.eps'
#md.mln_plt_mat_graph_SEEG(savefig,ifMat,mmRule,iarea_list,Rmat=[],figtitle=[])
    ifMatO=np.zeros([nchan,nchan])
    inchan=len(nodelabels)
    
    for iind, inode in enumerate(nodelabels):
        for jind,jnode in enumerate(nodelabels):
            ifMatO[inode-1,jnode-1]=ifMat[iind,jind]
    N20fMat.append(ifMatO)
N20fMat.append(fMat)

# Calcualte NfMat
NfMat=np.zeros([nchan,nchan])
sg=(np.arange(len(signalgroup))+1).astype(int)
sgnodelabels=[]
for isg in sg:
    sgnodelabels.append(np.arange(signalgroup[isg-1][0],signalgroup[isg-1][1]+1).astype(int))
print sgnodelabels

for iind, inode in enumerate(orginallabels):
    lind=mi.mln_find_subgroup(inode,sgnodelabels)
    for jind, jnode in enumerate(orginallabels):
        ac=sgnodelabels[lind]==jind
        ijfMatSum=[N20fMat[i][iind][jind] for i in np.arange(4)]
        if ac[ac.nonzero()]:
            NfMat[iind,jind]=np.sum(ijfMatSum)/3.
        
        else:
            
            NfMat[iind,jind]=np.sum(ijfMatSum)/2.

#
sfilename='S_'+basicnom+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)
savefig=savedir+'Graph_CV_'+sfilename+'.eps'
md.mln_plt_mat_graph_SEEG(savefig,NfMat,mmRule,area_list,Rmat=[],figtitle=[])