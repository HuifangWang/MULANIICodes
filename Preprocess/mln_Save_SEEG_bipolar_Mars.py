# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 00:26:44 2015

@author: huifangwang
"""


import sys, os.path
sys.path.append( os.path.expanduser('/Users/huifangwang/BHBS/NSeeg/Programs') )
import ns_tools as nt

import ns_read_eeg as nread
import mne
import ns_demonstration as ns

from array import array
import numpy as np

data_path = '/Users/huifangwang/MULAN3/mlnData/SEEG/MotorNetwork/GUERFI/SEEGdata/'
data_filename='130830B-CEX_0000'
raw_fname = data_filename+'_raw.fif'
Pname='GUERFI'

raw = mne.io.Raw(data_path+raw_fname) 
Seeg=raw

data, times = Seeg[:, :]

bip_data = []
bip_chan = []
 
#   16, 17, 'OF',   2,   3  results in 'OF2-3', _data[16] = _data[16] - _data[17]
for i1, i2, elec, ei1, ei2 in ns.find_bip_idx(Seeg):
    ch_info = Seeg.info['chs'][i1].copy()
    ch_info['ch_name'] = '%s%d-%s%d' % (elec, ei1, elec, ei2)
    bip_chan.append(ch_info)
    bip_data.append(data[i1] - data[i2])


chs = bip_chan
ch_names = [c['ch_name'] for c in chs]

### channel_area
Chosen_area=[('SA',10,1),('CC', 13,2),('OR',15,3)]
Chosen_ch=[]
Chosen_Color=[]
for iarea, valid_leed, areacolor in Chosen_area:
    Area_list=[iarea+str(i)+'-'+iarea+str(i+1) for i in np.arange(1,valid_leed)]
    Chosen_ch=Chosen_ch+Area_list
    vectorcolor=np.ones(valid_leed-1)*areacolor
    Chosen_Color=Chosen_Color+vectorcolor.tolist()
    
prenomfile=data_path+Pname+data_filename

periods=60
fs=raw.info['sfreq']

Chosen_area=[('SA',10,1),('CC', 13,2),('OR',15,3)]
vts=[1,400,600]
for ts in vts:
    nt.ns_mln_save2mat(ts,periods,fs,prenomfile,ch_names,bip_data,Chosen_area)
