function mln_adp_Cal(dirname,prenom,wins,inwins)
% Calculate the connectivity from MULAN algorithm based on the optimal
% parameters

% dirname: the folder to store all files
% prenom: spectial name for new datasets
% flagTF: We calculate Basic Methods flagTF=T; otherwise flagTF=F
% flagData: flagData='M' if the calculated datasets is simulated data,
% otherwise for the real data, the AUC and COS parts will be ignored
% 
% Examples:mln_Cal_aMULAN_LFP
% /Users/huifangwang/MULAN3/mlnData/LFP/SlowTheta/Win2500/ 20130516_CTRds2.mat F
% 2500 0 3
% -------------------------------------------------------------------------
% Huifang, Sep 5, 2014,generate the data, calculated the basic methods, and MULAN algorithm and distribution 
% Huifang Nov 4, 2014 From this version, I updated icov for PcorrU
% This version the results of MULAN updated
% Huifang Wang, Marseille, Nov 27, 2014
% Huifang Wang, Marseille, Jan, 30, 2015: MULAN algorithm outputs the final
% results for given parameters and group methods
% -------------------------------------------------------------------------


%VGroupMethlog={'TimeBasic','FreqBasic','Hsquare','Granger','FreqAH','TE'};

para.name = {'th11','th12','th21','th22','r1','r2','r3','r4','r5','thm'};
% this parameter is for N30L15
para.pE = [ 0.7410    0.5267    0.8203    0.8340    0.2089    0.1600    0.2274    0.1866    0.9417    0.6750];
para.var = [ 0.0221    0.0205    0.0009    0.0003    0.0066    0.0166    0.0119    0.0113    0.0010    0.0057];
para.methodlog = {'BCorrU','COH1','BCorrD','PDC'};



inwins=str2double(inwins);

add=[wins];
para.kforMF=5;
para.nboot=50;%%%
para.nboot2=20;%%%
%BM.G12={'BCorrU','PCorrU';'BTEU','PTEU'};
%BM.G12={'BCorrU','PCorrU';'BTEU','PTEU'};
%BM.G3={'BCorrD','BTED'};
%BM.G4={'ffDTF','hmvar'};
% BMG={{'BCorrU','PCorrU','BCorrD','hmvar'},{'BTEU','PTEU','BCorrD','ffDTF'},{'COH1','BTEU','GGC','hmvar'}};


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

dataname=[add,'_',dataname];
mln_adp_inwins(dirname,dataname,para,inwins,add);



