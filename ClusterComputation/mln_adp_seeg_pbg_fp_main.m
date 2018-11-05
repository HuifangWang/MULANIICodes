function mln_adp_seeg_pbg_main(dirname,prenom,wins,inwins,parafile)

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
%dirname=[dirname,'/Win',wins,'/'];
[~,Nopt]=size(opt);
for indG = 1:Nopt
    %for indtype={'mean','media','hp','tak'}
        %indtype=char(indtype);
        indtype='takfp';
        saveadd=['FRN_G', num2str(opt(indG).iG),'_',indtype];
        %para.pE = pE(iG,:);
        %para.var = var(iG,:);%[ 0.0221    0.0205    0.0009    0.0003    0.0066    0.0166    0.0119    0.0113    0.0010    0.0057];
        para.methodlog = opt(indG).BMs;%{'BCorrU','COH1','BCorrD','PDC'};
        iopt=opt(indG);
        
        mln_adp_inwins(dirname,dataname,para,inwins,saveadd,indtype,iopt);
    
end