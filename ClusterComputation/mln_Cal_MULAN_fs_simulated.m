function mln_Cal_MULAN_fs_simulated(dirname,prenom,threshold1,threshold2,flagTF,flagData)
% Huifang, Sep 5, 2014,generate the data, calculated the basic methods, and MULAN algorithm and distribution 
% Huifang Nov 4, 2014 From this version, I updated icov for PcorrU

%Examples:mln_Cal_MULAN_r_nboots tryhn1 Tp-try 0.5 0.7 T R
% Huifang Wang, Marseille, Nov 25, 2013,  Calculate all
% datasets, Get the AUC

%VGroupMethlog={'TimeBasic','FreqBasic','Hsquare','Granger','FreqAH','TE'};
VGroupMethlog={'TimeBasic','TE','FreqAH'};
%VGroupMethlog={'TimeBasic','FreqAH'};
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


calParams.defwindow=1000;
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

indmat=strfind(prenom,'.mat');
if ~isempty(indmat)
 dataname=prenom(1:indmat-1);
else
dataname=prenom;
end

if flagTF=='T'

MULAN_Cal_BasicMethods(dirname,dataname,calParams,VGroupMethlog);
if flagData=='M'
mln_AUC_Tout(dirname,dataname);
end
else
outputname=mln_Stat_Valid_thv(dirname,dataname,threshold1,threshold2,kforMF,Fthreshold,nboot,nboot2,BM);
if flagData=='M'
mln_Cos_mulan(dirname,outputname);
end
end
