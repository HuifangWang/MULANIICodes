function mln_Cal_MULAN_fs_mlnvcs(dirname,prenom,th1,th2,flagTF,flagData)
% Calculate the connectivity from both basic methods and MULAN algorithm

% dirname: the folder to store all files
% prenom: spectial name for new datasets
% th1: the parameter for BM1-2
% th2: the parameter for BM3-4
% flagTF: We calculate Basic Methods flagTF=T; otherwise flagTF=F
% flagData: flagData='M' if the calculated datasets is simulated data,
% otherwise for the real data, the AUC and COS parts will be ignored
% 
% Examples:mln_Cal_MULAN_fs_mlnvcs Examples nmm 0.5 0.7 T M
% -------------------------------------------------------------------------
% Huifang, Sep 5, 2014,generate the data, calculated the basic methods, and MULAN algorithm and distribution 
% Huifang Nov 4, 2014 From this version, I updated icov for PcorrU
% This version the results of MULAN updated
% Huifang Wang, Marseille, Nov 27, 2014
% -------------------------------------------------------------------------


%VGroupMethlog={'TimeBasic','FreqBasic','Hsquare','Granger','FreqAH','TE'};
VGroupMethlog={'TimeBasic','TE','FreqAH'};
%VGroupMethlog={'TimeBasic','FreqAH'};



calParams.defwindow=300;
calParams.defoverlap=0.5;
calParams.defmodelorders=5;
calParams.minfreq=1;
calParams.maxfreq=40;
calParams.stepfreq=1;
calParams.defMaxDelay=10;
calParams.defbins=16;

kforMF=5;
Fthreshold=0.7;
nboot=1000;
nboot2=20;
%BM.G12={'BCorrU','PCorrU';'BTEU','PTEU'};
BM.G12={'BCorrU','PCorrU';'BTEU','PTEU'};
BM.G3={'BCorrD','BTED'};
BM.G4={'ffDTF','hmvar'};


datadir=[dirname,'/data'];
if ~exist(datadir,'dir')
    display('where did the data go?')
end

Resultsdir=[dirname,'/Results'];
if ~exist(Resultsdir,'dir')
    mkdir(Resultsdir);
end

Resultsdir=[dirname,'/ToutResults'];
if ~exist(Resultsdir,'dir')
    mkdir(Resultsdir);
end

indmat=strfind(prenom,'.mat');
if ~isempty(indmat)
 dataname=prenom(1:indmat-1);
else
dataname=prenom;
end

if flagTF=='T'

mln_Cal_BasicMethods(dirname,dataname,calParams,VGroupMethlog);
if flagData=='M'
mln_AUC_Tout(dirname,dataname);
end
else
outputname=mln_Stat_Valid_thv_mlnvcs(dirname,dataname,th1,th2,kforMF,Fthreshold,nboot,nboot2,BM);
if flagData=='M'
mln_Cos_mulan(dirname,outputname);
end
end
