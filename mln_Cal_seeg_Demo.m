%function mln_Cal_seeg_Demo()

dirname = './Examples/SeegData/';
prenom = 'SEEG_PS1001ts0ds4.mat';
wins = '1500';

%% Calculate the Basic Methods
VGroupMethlog={'TimeBasic','TE','FreqAH'};

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
%% Calculate the MULAN results

parafile=[dirname,'opt_seeg_HHS30.mat'];
inwins='-1';
para.name = {'th11','th12','th21','th22','r1','r2','r3','r4','r5','thm'};

inwins=str2double(inwins);

add=[wins];
para.kforMF=5;
para.nboot=50;%%%
para.nboot2=20;


load(parafile,'opt')
dataname=[wins,'_',dataname];

indG = 2;

indtype='tak';
saveadd=['FRN_G', num2str(opt(indG).iG),'_',indtype];
para.methodlog = opt(indG).BMs;
iopt=opt(indG);

mln_adp_inwins(dirname,dataname,para,inwins,saveadd,indtype,iopt);
