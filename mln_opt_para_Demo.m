function mln_opt_para_Demo()

dirname='./Examples/AdapData/';
filedata = 'erpN30S10N20100ds4.mat';
wins='1500';
prenom = filedata;
%%Calculate Basic Methods
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







wins=str2double(wins);
winssub=[dirname,'Win', num2str(wins)];
if exist(winssub,'dir')
    
traindatafile = [dirname,'Win', num2str(wins), '/ToutResults/Tout_', num2str(wins),'_',filedata];
else
 traindatafile = [dirname, '/ToutResults/Tout_', num2str(wins),'_',filedata];
   
end
para.name = {'th11','th12','th21','th22','r1','r2','r3','r4','r5','thm'};
para.p0=[0.5,0.5,0.6,0.6,0,0,0,0,1,0.5];
para.rangeMin=[0.2,0.2,0.3,0.3,0,0,0,0,0.8,0.5];
para.rangeMax=[0.9,0.9,0.9,0.9,0.5,0.5,0.5,0.5,1,0.9];
fitness={'FPN'};

%methodlog={'BCorrU','COH1','BCorrD','PDC'};
methodgroup={{'BTEU','BCorrU','hmvar','ffDTF'},{'BCorrU','BTEU','BTED','ffDTF'},...,
    {'BCorrU','PCorrU','hmvar','ffDTF'},{'BCorrU','PCorrU','BTED','hmvar'},{'BCorrU','COH1','BCorrD','PDC'}};

for ibg=[2,5]%1:length(methodgroup)
para.fitness=fitness{1};
methodlog=methodgroup{ibg};
savefilefit = [dirname,'/ToutResults/Fit_',para.fitness, 'G',num2str(ibg),'_',num2str(wins),'_',filedata];
if ~exist(savefilefit,'file')
   mln_opt_para_vcs_fpn(traindatafile,para,methodlog,savefilefit);
end
end

