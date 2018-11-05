function mln_adp_Cal_pbg_main_SEEG(dirname,prenom,wins,inwins,parafile)
% Calculate the connectivity from MULAN algorithm based on the optimal
% parameters; only 1 way to take the optimal parameters
%dirname='/Users/huifangwang/MULAN3/mlnData/dyN20Hz125Mj12wins/';
%prenom = 'nmmN20L10CS80100S10Nt12Nstru4.mat';
%wins=num2str(500);
%inwpara.thetains=num2str(1);
%parafile=['/Users/huifangwang/MULANa/mlnData/N20CS8S50/opt_para_N20CS8L15.mat'];

wins=num2str(1500);
inwins=num2str(-1);
dirname='/Users/huifangwang/MULANa/SeegData/PS1001/';
parafile=[dirname,'opt_para_N30CS1L15.mat'];
prenom = 'WB_PS1001WI041222014ts0ds4.mat';

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

%parafile='opt_para_N30CS8L10.mat'

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
%methodgroup={{'BTEU','BCorrU','hmvar','ffDTF'},{'BCorrU','BTEU','BTED','ffDTF'},...,
  %  {'BCorrU','PCorrU','hmvar','ffDTF'},{'BCorrU','PCorrU','BTED','hmvar'},{'BCorrU','COH1','BCorrD','PDC'}};
load(parafile,'opt')
dataname=[add,'_',dataname];
dirname=[dirname,'/Win',wins,'/'];
[~,Nopt]=size(opt);
for indG = 1:Nopt
    %for indtype={'mean','media','hp','tak'}
        %indtype=char(indtype);
        indtype='tak';
        saveadd=['_G', num2str(opt(indG).iG),'_',indtype];
        %para.pE = pE(iG,:);
        %para.var = var(iG,:);%[ 0.0221    0.0205    0.0009    0.0003    0.0066    0.0166    0.0119    0.0113    0.0010    0.0057];
        para.methodlog = opt(indG).BMs;%{'BCorrU','COH1','BCorrD','PDC'};
        iopt=opt(indG);
        
        mln_adp_inwins_vcs(dirname,dataname,para,inwins,saveadd,indtype,iopt);
    
end