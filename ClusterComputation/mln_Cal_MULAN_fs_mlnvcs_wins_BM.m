function mln_Cal_MULAN_fs_mlnvcs_wins_BM(dirname,prenom,wins)
% Calculate the connectivity from both basic methods and MULAN algorithm

% dirname: the folder to store all files
% prenom: spectial name for new datasets
% th1: the parameter for BM1-2
% th2: the parameter for BM3-4
% flagTF: We calculate Basic Methods flagTF=T; otherwise flagTF=F
% flagData: flagData='M' if the calculated datasets is simulated data,
% otherwise for the real data, the AUC and COS parts will be ignored
% 

%wins=num2str(1500);
%dirname='/Users/huifangwang/MULANa/humanCon/HC30S50T/';
%prenom='erpN30S2N20100ds4';

%dirname='/Users/huifangwang/MULANa/SeegData/PS1001/';

%parafile=[dirname,'opt_para_N20CS8L15.mat'];
%prenom = 'WB_PS1001WI041222014ts0ds4.mat';
%% example for MULAN II paper  
%dirname='/Users/huifangwang/MULANa/SeegData/PS1001//hidNode/';
% is=3;
% prenom = ['WB_Hd',num2str(is),'PS1001WI041222014ts0ds4.mat'];
% flagTF='T';
% flagData='T';
%% for Cleveland data
%dirname='/Users/huifangwang/NormalSEEG_CL/PS1001/mln/';
%ind=3;
%prenom = ['ND',num2str(ind),'_PS1001WI041222014ts200ds4.mat'];
%flagTF='T';
%flagData='T';

% Examples:mln_Cal_MULAN_fs_mlnvcs_wins Examples nmmN20L10CS80CS100S1N8100 T M 800
% -------------------------------------------------------------------------
% Huifang, Sep 5, 2014,generate the data, calculated the basic methods, and MULAN algorithm and distribution 
% Huifang Nov 4, 2014 From this version, I updated icov for PcorrU
% This version the results of MULAN updated
% Huifang Wang, Marseille, Nov 27, 2014
% Huifang Wang, Marseille, Jan, 30, 2015: MULAN algorithm outputs the final
% results for given parameters and group methods
% -------------------------------------------------------------------------


%VGroupMethlog={'TimeBasic','FreqBasic','Hsquare','Granger','FreqAH','TE'};
VGroupMethlog={'TimeBasic','TE','FreqAH'};
%VGroupMethlog={'TimeBasic','FreqAH'};



calParams.defwindow=str2double(wins);
calParams.defoverlap=0.5;
calParams.defmodelorders=5;
calParams.minfreq=1;
calParams.maxfreq=40;
calParams.stepfreq=1;
calParams.defMaxDelay=10;
calParams.defbins=16;


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


mln_Cal_BasicMethods_wins(dirname,dataname,calParams,VGroupMethlog,wins);
