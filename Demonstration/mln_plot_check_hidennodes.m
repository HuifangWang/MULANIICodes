function mln_plot_check_hidennodes()

%dir='/Volumes/HWMac1/MULAN2/SPM/hiddennodes/hERPN20l/hERPN20lwins600/data/';
%load([dir,'N20nmmCS80S5N4100N19hnl.mat']);
%% following for Paper fig.5 
% for high degree
%dir='/Volumes/HWMac1/MULAN2/SPM/hnNmmN30L10h/data/';
%load([dir,'nmmN30L10CS100S6N6100N25hnh.mat']);

% for low degree
dir='../Examples/forhn/';
load([dir,'nmmN30L10CS80S1N6100N28hnl.mat']);


subNnodes=5;
thresholdnw=0.01;
Params.fs=125;


Mat=Connectivity.pMat;
tMat=Connectivity.label;
hf=figure('Name','Check hidden nodes');
N=length(Mat);
Nrow=2;
Ncol=ceil(N/Nrow);
ha=subplot(Nrow,Ncol,[1,1+Ncol],'Parent',hf);
mln_showNetworkGraph_L(Mat{1},ha,subNnodes,thresholdnw,Params)
title(ha,tMat{1});
Np=[1:3,(4:5)+1];
for iN=2:N
    ha=subplot(Nrow,Ncol,Np(iN),'Parent',hf);
    mln_showNetworkGraph_hn(Mat{iN},ha,subNnodes,thresholdnw,Params,Connectivity.delnodes)
    title(ha,tMat{iN});
end