function mln_hn_Gen_from_mch(basicdir)
% Divided the datasets with 30 nodes to 7 datasets with given subnetworks (3 N10 + 3 N20)

% basicdir: the folder name to store all files

% -----------------------------------------------------
% Examples:mln_hn_Gen_from_mch 
% -------------------------------------------------------------------------
% Huifang Wang, Marseille, Nov 25, 2013,  
%basicdir='/Users/huifangwang/MULANa/mlnData/FRN30CS1S30/';
basicdir='/Users/huifangwang/MULANa/humanCon/HH30S50T/'

olddatadir=[basicdir,'data/'];
newdatadir=[basicdir,'hidden/data/'];

list=dir(fullfile(olddatadir,'*.mat'));


if ~exist(newdatadir,'dir')
    mkdir(newdatadir);
end

for ilist=1:length(list)
    prenom=list(ilist).name;
    indmat=strfind(prenom,'.mat');
if ~isempty(indmat)
 dataname=prenom(1:indmat-1);
else
dataname=prenom;
end

ofile=load([olddatadir,dataname]);
signalgroup={[1:10],[11:20],[21:30]};
%% divided the datasets with every 10 nodes
%for isg=1:length(signalgroup)
 %   Data=ofile.Data(signalgroup{isg},:);
  %  Connectivity=ofile.Connectivity(signalgroup{isg},signalgroup{isg});
   % Params=ofile.Params;
   % Params.str=num2cell(signalgroup{isg});
   % save([newdatadir,dataname,'sg',num2str(isg),'.mat'],'Data','Connectivity','Params')
%end  

%% divided the datasets with every 20 nodes
 isgv=1:length(signalgroup);
 bgv=nchoosek(isgv,2);
% 
 for ibgv=1:length(bgv)
    clearvars Params Data Connectivity
    indexsg=[signalgroup{bgv(ibgv,1)},signalgroup{bgv(ibgv,2)}];
    Data=ofile.Data(indexsg,:);
    Connectivity=ofile.Connectivity(indexsg,indexsg);
    Params=ofile.Params;
    Params.str=num2cell(indexsg);
    sfilename=[newdatadir,dataname,'b',num2str(bgv(ibgv,1)),num2str(bgv(ibgv,2)),'.mat']
    save(sfilename,'Data','Connectivity','Params')
end
end