# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 13:41:19 2014

@author: huifangwang
call mln_tools:
sys.path.append(os.path.expanduser('/Users/huifangwang/MULANIII/Codes/BasicTools/'))
import mln_tools as mt
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.colors as cm

def mln_SEEG_save2mat(ts,periods,fs,prenomfile,ch_names,bip_data,Chosen_area):
    """ save the normal SEEG data into mln format
    """
    Chosen_ch=[]
    AllChosen_area=''
    for iarea, valid_leed, areacolor in Chosen_area:
        Area_list=[iarea+str(i)+'-'+iarea+str(i+1) for i in np.arange(1,valid_leed)]
        Chosen_ch=Chosen_ch+Area_list
        AllChosen_area=AllChosen_area+iarea
    sub_EEG=np.array(nd.pick_channels(listregion=ch_names,chosen_ch=Chosen_ch,data=bip_data))
    Data=sub_EEG[:,int(ts*fs):int((ts+periods)*fs)]
    Params={'fs':fs,'str':Chosen_ch}
    filename=prenomfile+AllChosen_area+str(len(Chosen_ch))+'ts'+str(ts)+'.mat'
    a={}
    a['Data']=Data
    a['Params']=Params
    sio.savemat(filename,a)

def mln_save2mat_4cal(Data,fs,Chosen_ch,givenfilename):
    Params={'fs':fs,'str':Chosen_ch}
    filename=givenfilename+'.mat'
    a={}
    a['Data']=Data
    a['Params']=Params
    sio.savemat(filename,a)


def mln_SEEG_extact_filename(ts,prenomfile,postfix,Chosen_area):
    Chosen_ch=[]
    AllChosen_area=''
    for iarea, valid_leed, areacolor in Chosen_area:
        Area_list=[iarea+str(i)+'-'+iarea+str(i+1) for i in np.arange(1,valid_leed)]
        Chosen_ch=Chosen_ch+Area_list
        AllChosen_area=AllChosen_area+iarea
    filename=prenomfile+AllChosen_area+str(len(Chosen_ch))+'ts'+str(ts)+postfix
    return filename

def mln_read_wins_AUC(filename):
    matdata = sio.loadmat(filename)
    Meths=matdata['Meths']
    MSAUC=Meths['MSAUC'][0,0]
    bMs=Meths['Methodnames']
    sM=[bMs[0][0][i][0][0] for i in np.arange(np.size(bMs[0][0]))]

    return MSAUC, sM
    
    
def mln_read_wins_Sta(filename):
    matdata = sio.loadmat(filename)
    Eva=matdata['Eva']
    xinfo=matdata['xinfo']
    Nwins=xinfo['Nwins']
    EvaAUC=Eva['AUC']
    BMs=xinfo['BMG']
    bMG=BMs[0][0][0]
    
    calNWins=Nwins[0][0][0]
    Ncalwins=np.size(calNWins)
    EvaAUCr=[np.max(EvaAUC[0][i]) for i in np.arange(Ncalwins)]
        
    return EvaAUCr, calNWins, bMG
    
def mln_plot_3D2D_boxes(ax,data,myxtick):
    # data is 3D
    #fig.canvas.set_window_title('A Boxplot Example')
    #plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    nbars,nx,ndis = np.shape(data)
    data=np.reshape(data,[nx*nbars,ndis],order='F')
    data=np.swapaxes(data,0,1)
    box_widths=0.7
    ibarwidths=box_widths/nbars*0.9
    inerror=np.linspace(0,box_widths,nbars)-box_widths/2.
    error=np.tile(inerror,nx)
    bposition=np.arange(0,nx*nbars)/nbars+error
    bp = ax.boxplot(data, widths=ibarwidths,positions=bposition, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'],color='black')
    plt.setp(bp['whiskers'],color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
#
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.5)
#
#    # Hide these grid behind plot objects
    ax.set_axisbelow(True)
    
#    #ax1.set_xlabel('Distribution')
#    #ax1.set_ylabel('AUC Value')
    boxColors=[plt.cm.spectral(ibars/float(nbars)) for ibars in np.arange(nbars)]
#    # Now fill the boxes with desired colors
    #boxColors = ['ForestGreen','royalblue']
    #mColors=['DarkGreen','DarkBlue']
    mColors=boxColors
    numBoxes = nx*nbars
    medians = range(numBoxes)
    datamean=data.mean(axis=0)
    for i in range(numBoxes):
      box = bp['boxes'][i]
      boxX = []
      boxY = []
      for j in range(5):
          boxX.append(box.get_xdata()[j])
          boxY.append(box.get_ydata()[j])
      boxCoords = zip(boxX,boxY)
      # Alternate between Dark Khaki and Royal Blue
      k = i % nbars
      boxPolygon = plt.Polygon(boxCoords, facecolor=boxColors[k])
      ax.add_patch(boxPolygon)
      # Now draw the median lines back over what we just filled in
#     med = bp['medians'][i]
#      medianX = []
#      medianY = []
#      for j in range(2):
#          medianX.append(med.get_xdata()[j])
#          medianY.append(med.get_ydata()[j])
#          plt.plot(medianX, medianY, color=mColors[k] ,lw=3)
#          medians[i] = medianY[0]
#      # Finally, overplot the sample averages, with horizontal alignment
#      # in the center of each box
#      plt.plot([np.average(med.get_xdata())], [datamean[i]],
#               color=mColors[k], marker='o', markeredgecolor='k')
#
#                # Set the axes ranges and axes labels
    #ax.set_xlim(0, numBoxes+0.5)
    top = 1.01
    bottom = 0.6
    ax.set_ylim(bottom, top)
#    xtickNames = plt.setp(ax1, xticklabels=np.tile(['S1','S2'], nhn))
#    plt.setp(xtickNames, fontsize=14)
#    indxticks=np.arange(1.5,nhn*nrm+1,nrm)
    ax.set_xticks(np.arange(nx))
    ax.set_xticklabels(myxtick,fontsize=14)
    
def mln_legend(ax,colors,lables):
    dwidth=0.5
    dhight=0.1
    for ind,ilable in enumerate(lables):
        xi=[ind,ind,ind+dwidth,ind+dwidth]
        yi=[dhight,0,0,dhight]
        boxCoords = zip(xi,yi)
        #print xi
        boxpolygon=plt.Polygon(boxCoords, facecolor=colors[ind])
        ax.add_patch(boxpolygon)
        plt.text(ind+dwidth/2.,dhight*2,ilable,fontsize=14,ha='center')
    ax.set_xlim(-1,len(lables))
    ax.set_ylim(-0.2,0.6)
    
    ax.set_xticks([])
    ax.set_yticks([])

def mln_norm(CorrM,mode='diag'):
    CorrMd=CorrM-np.diag(np.diag(CorrM)) #
    CorrMd=CorrMd/float(np.max(CorrMd))
    CorrMd=CorrMd+np.diag(np.diag(CorrM))
    return CorrMd


def mln_calCorrM(BMwins,order=1):
    nchan,nchan,Nwins=np.shape(BMwins)
    Nwins=Nwins-order+1
    CorrM=np.zeros([Nwins,Nwins])
    for iBMwins in np.arange(Nwins):
       for jBMwins in np.arange(Nwins):
           if order==1:
              iBM=BMwins[:,:,iBMwins]
              jBM=BMwins[:,:,jBMwins]
           else:
              iBM=np.mean(BMwins[:,:,iBMwins:(iBMwins+order-1)],axis=2)
              jBM=np.mean(BMwins[:,:,jBMwins:(jBMwins+order-1)],axis=2)
           iBM=iBM.flatten()
           jBM=jBM.flatten()
           ijcorr=np.corrcoef(iBM,jBM)
           CorrM[iBMwins,jBMwins]=ijcorr[0,1]
    return CorrM

def mln_calCosM(BMwins,order=1):
    nchan,nchan,Nwins=np.shape(BMwins)
    Nwins=Nwins-order+1
    CosM=np.zeros([Nwins,Nwins])
    for iBMwins in np.arange(Nwins):
        for jBMwins in np.arange(Nwins):
            if order==1:
                iBM=BMwins[:,:,iBMwins]
                jBM=BMwins[:,:,jBMwins]
            else:
                iBM=np.mean(BMwins[:,:,iBMwins:(iBMwins+order-1)],axis=2)
                jBM=np.mean(BMwins[:,:,jBMwins:(jBMwins+order-1)],axis=2)
            ijcorr=mln_Sim_Cos(iBM,jBM)
            CosM[iBMwins,jBMwins]=ijcorr
    return CosM

def mln_Sim_Cos(M1,M2):
    vM1=M1.flatten()
    vM2=M2.flatten()
    normA=np.sqrt(sum(vM1**2))
    normB=np.sqrt(sum(vM2**2))
    return np.sum(vM1*vM2)/float(normA*normB)


def ns_channel_electro_color(area_color):
    electrodes_color=[]
    electrodes=[]
    for iarea,leeds,fcolor,lcolor in area_color:
        nleeds=(leeds[1]-leeds[0])/leeds[2]+1
        dicolor=(np.array(fcolor)-np.array(lcolor))/float(np.array(nleeds))
        #print fcolor, lcolor
        
        for indl,ileeds in enumerate(np.arange(leeds[0],leeds[1]+leeds[2],leeds[2])):
            icolor=dicolor*indl+lcolor
            #print icolor
            electrodes_color.append(icolor/255.)
            electrodes.append(iarea+str(ileeds))
    return electrodes_color, electrodes

def mln_nodeColors(levelmap,dV,lcolorname='white'):
    node_color=[]
    fcolor=hex_to_rgb(cm.cnames[levelmap])
    lcolor=hex_to_rgb(cm.cnames[lcolorname])
    nnode=len(dV)
    nleeds=20
    
    dicolor=(np.array(fcolor)-np.array(lcolor))/float(np.array(nleeds))
    for inode in np.arange(nnode):
        icolor=dicolor*dV[inode]/np.max(dV)*nleeds+lcolor
        node_color.append(icolor/255.)
    return node_color


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


