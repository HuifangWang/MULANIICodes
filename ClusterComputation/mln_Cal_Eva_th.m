function mln_Cal_Eva_th(dirname,prenom,strfile,paramsfile,nc,is,npts,cs,models,Methodthresholdfile)
% Generate the data, calucate the connection methods by given methods, evalution the results by AUC, given the thresholds.
% dirname: the folder to store all files
% prenom: spectial name for new datasets
% strfile: the file stores the structures
% Paramsfile: the params used in the calucation
% nc: the channels of data
% is: which structure stored in the structure file
% npts: length of data
% cs: the connection strength
% models: 'nmm','fMRI','linear','rossler', 'henon'

%Examples: mln_Cal_Eva_th nmm nmm structureN5L5 nmmParams 5 20 200 1 nmm SetMethodsTh
% dirname is the direction name which

% Huifang Wang, Marseille, Nov 25, 2013,  Calculate all
% datasets, Get the AUC
VGroupMethlog={'TimeBasic','FreqBasic','Hsquare','Granger','FreqAH'};
%restoredefaultpath;
%path=cd;
%addpath(genpath(path));
%savepath;

if ~exist(dirname,'dir')
    mkdir(dirname);
end
datadir=[dirname,'/data'];
if ~exist(datadir,'dir')
    mkdir(datadir);
end
Resultsdir=[dirname,'/Results'];
if ~exist(Resultsdir,'dir')
    mkdir(Resultsdir);
end
Resultsdir=[dirname,'/ToutResults'];
if ~exist(Resultsdir,'dir')
    mkdir(Resultsdir);
end


%mln_generateParams(dirname);
%if exist([paramsfile,'.mat'],'file')
  % strfile1=which([paramsfile,'.mat']); 
    
%strfilethere=['./',dirname,'/',strfile,'.mat'];

copyfile ([paramsfile,'.mat'], ['./',dirname]);
%end

%if exist([Methodthresholdfile,'.mat'],'file')
  % strfile2=which([Methodthresholdfile,'.mat']); 
    

copyfile ([Methodthresholdfile,'.mat'], ['./',dirname]);
%end

is=str2double(is);
nc=str2double(nc);
npts=str2double(npts);
cs=str2double(cs);
switch models
    case 'nmm'
    dataname=mln_generate_nmm(is,nc,strfile,dirname,prenom,cs,npts); 
    case 'fMRI'
   dataname=mln_generate_fMRI(is,nc,strfile,dirname,prenom,cs,npts); 
end

MULANCalMUltiBP(dirname,dataname,paramsfile,VGroupMethlog);
mln_MULAN_th(nc,dirname,dataname,Methodthresholdfile);
%mln_rule_Validation(dirname,prenom)
mln_MethodStructuresAUC(dirname,dataname);
display('finished!');
error('er')

