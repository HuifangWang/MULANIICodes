function mln_Cal_LFP_swin(dirname,prenom,flagTF,wins,overlap,inwins,imlnwins)
% Calculate the connectivity from both basic methods and MULAN algorithm

% dirname: the folder to store all files
% prenom: spectial name for new datasets
% th1: the parameter for BM1-2
% th2: the parameter for BM3-4
% flagTF: We calculate Basic Methods flagTF=T; otherwise flagTF=F
% flagData: flagData='M' if the calculated datasets is simulated data,
% otherwise for the real data, the AUC and COS parts will be ignored
% inwins: MULAN takes how many windows to integrate
% imlnwins: the number order of MULAN

% Examples:mln_Cal_LFP_swin /Users/huifangwang/MULAN3/mlnData/LFP/SlowTheta/Win2500 20130516_CTRds2.mat F 
% 2500 0 3 3
% -------------------------------------------------------------------------
% Huifang, Sep 5, 2014,generate the data, calculated the basic methods, and MULAN algorithm and distribution 
% Huifang Nov 4, 2014 From this version, I updated icov for PcorrU
% This version the results of MULAN updated
% Huifang Wang, Marseille, Nov 27, 2014
% Huifang Wang, Marseille, Jan, 30, 2015: MULAN algorithm outputs the final
% results for given parameters and group methods
% -------------------------------------------------------------------------
% 

%VGroupMethlog={'TimeBasic','FreqBasic','Hsquare','Granger','FreqAH','TE'};
VGroupMethlog={'TimeBasic','TE','FreqAH','FreqBasic'};
%VGroupMethlog={'TimeBasic'};


inwins=str2double(inwins);
imlnwins=str2double(imlnwins);
calParams.defwindow=str2double(wins);
calParams.defoverlap=str2double(overlap)/100;
calParams.defmodelorders=5;
calParams.minfreq=1;
calParams.maxfreq=50;
calParams.stepfreq=1;
calParams.defMaxDelay=10;
calParams.defbins=16;
add=[wins,'_',overlap];
kforMF=5;
Fthreshold=0.7;
nboot=100;%%%
nboot2=20;%%%
th1=0.5:0.05:0.6;
th2=0.6:0.05:0.7;
%BM.G12={'BCorrU','PCorrU';'BTEU','PTEU'};
%BM.G12={'BCorrU','PCorrU';'BTEU','PTEU'};
%BM.G3={'BCorrD','BTED'};
%BM.G4={'ffDTF','hmvar'};
BMG={{'BCorrU','PCorrU','BCorrD','hmvar'},{'BTEU','PTEU','BCorrD','ffDTF'},{'COH1','BTEU','GGC','hmvar'}};


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

mln_Cal_BasicMethods_wins(dirname,dataname,calParams,VGroupMethlog,add);

else
dataname=[add,'_',dataname];
mln_alg_output_dy_inwins_swins(dirname,dataname,th1,th2,kforMF,Fthreshold,nboot,nboot2,BMG,inwins,add,imlnwins);
end


