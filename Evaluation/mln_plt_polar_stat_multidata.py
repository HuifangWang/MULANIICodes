import scipy.io
import os.path
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
sys.path.append(os.path.expanduser('/Users/huifangwang/MULANII/programs/Demonstrate/'))
import mln_plt_demonstration as md
sys.path.append(os.path.expanduser('/Users/huifangwang/MULANII/programs/Inference/'))
import mln_cal_mat_mltp as mi


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
    elif option == '-th1s': th1s=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th1f': th1f=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th2s': th2s=int(sys.argv[1]); del sys.argv[1]
    elif option == '-th2f': th2f=int(sys.argv[1]); del sys.argv[1]
    elif option == '-pltflag': pltflag=bool(sys.argv[1]); del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

savedir=basicdir+'/MlnResults/'

if not os.path.exists(savedir):
    os.makedirs(savedir)

stru=np.arange(1,10.5).astype(int)
for istru in stru:
    basicnom='nmmN10L15CS10100S'+str(istru)+'N6100'

    mdMat, mnMat, fMat,Rmat=mi.mln_cal_mat_mltp(basicdir,basicnom,Gs,th1s,th1f,istru,th2s,th2f,mmRule,pltflag)
    tptnVectorR, tptnVector, pnVector=mi.mln_tptn(fMat,Rmat)
    if istru==1:
        All_pnVector=pnVector
        All_tptnVector=tptnVector
    All_pnVector=All_pnVector+pnVector
    All_tptnVector=All_tptnVector+tptnVector

All_tptnVectorR=All_tptnVector/All_pnVector
sfilename='All_'+basicnom+'_th1_'+str(th1s)+'-'+str(th1f)+'_th2_'+str(th2s)+'-'+str(th2f)
savefig=savedir+'TPTN_'+sfilename+'.pdf'
md.mln_polar_bar(savefig,All_tptnVectorR)