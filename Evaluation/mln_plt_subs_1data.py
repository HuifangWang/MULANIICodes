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
import itertools

Gs=[(0,0,1),(1,0,0)]
th1s=50
th2s=60
th1f=60
th2f=70
basicdir='/Volumes/Lioncub/Users/huifangwang/MULAN2/SPM/nmm16sN30Vcsds/'

pltflag=0
istru=1
forceRMat=0
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
        print 'Options invalides :',option,'->',sys.argv[0]



mmRule=[(2,[0.95],[0.95]),(1,[0.92],[0.92])]
savedir=basicdir+'/MlnResults/'+'S'+str(istru)+'/'+'Th1_'+str(th1s)+str(th1f)+'_Th2_'+str(th2s)+str(th2f)+'mmRule'+str(mmRule[1][1])+'/'

if not os.path.exists(savedir):
    os.makedirs(savedir)

# for N30
 #(mlnValue,medianScore, meanScore)

basicnom='nmmN30L15CS10100S'+str(istru)+'N40100ds5'
mi.mln_plt_mat_mltp(basicdir,savedir,basicnom,Gs,th1s,th1f,istru,th2s,th2f,mmRule,forceRMat=forceRMat)

# for N10*3
mmRule=[(2,[0.8],[0.8]),(1,[0.4],[0.4])] #(mlnValue,medianScore, meanScore)

signalgroup=[(1,10),(11,20),(21,30)];
sg=(np.arange(len(signalgroup))+1).astype(int)
for isg in sg:
    basicnomsg=basicnom+'sg'+str(isg)
    nodelabels=np.arange(signalgroup[isg-1][0],signalgroup[isg-1][1]+1).astype(int)
    
    mi.mln_plt_mat_mltp(basicdir,savedir,basicnomsg,Gs,th1s,th1f,istru,th2s,th2f,mmRule,nodelabels=nodelabels,forceRMat=forceRMat)

# for N20*3
mmRule=[(2,[0.95],[0.95]),(1,[0.9],[0.9])] #(mlnValue,medianScore,meanScore)
ilist=itertools.combinations(sg,2)
seticom=set(ilist)
for ilistCa in seticom:
    basicnomsg=basicnom+'b'+str(ilistCa[0])+str(ilistCa[1])
    nodelabels=np.concatenate([np.arange(signalgroup[ilistCa[0]-1][0],signalgroup[ilistCa[0]-1][1]+1).astype(int),np.arange(signalgroup[ilistCa[1]-1][0],signalgroup[ilistCa[1]-1][1]+1).astype(int)])
    mi.mln_plt_mat_mltp(basicdir,savedir,basicnomsg,Gs,th1s,th1f,istru,th2s,th2f,mmRule,nodelabels=nodelabels,forceRMat=forceRMat)



