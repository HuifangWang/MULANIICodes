function allprenom=mln_gen_hidden_givendatasets(prenom,olddatadir,newdatadir,Nsubs,subNchan)

indmat=strfind(prenom,'.mat');
if ~isempty(indmat)
    dataname=prenom(1:indmat-1);
else
    dataname=prenom;
end

ofile=load([olddatadir,dataname]);
signalgroup=cell(1,Nsubs);
for isubs= 1:Nsubs
    signalgroup{isubs}=[1:subNchan]+(isubs-1)*subNchan;
end


%% divided the datasets with every 2 nodes
isgv=1:length(signalgroup);
bgv=nchoosek(isgv,2);
allprenom=cell(2,length(bgv));
for ibgv=1:length(bgv)
    clearvars Params Data Connectivity
    indexsg=[signalgroup{bgv(ibgv,1)},signalgroup{bgv(ibgv,2)}];
    Data=ofile.Data(indexsg,:);
    Connectivity=ofile.Connectivity(indexsg,indexsg);
    Params=ofile.Params;
    Params.str=num2cell(indexsg);
    sfilename=[dataname,'b',num2str(bgv(ibgv,1)),num2str(bgv(ibgv,2)),'.mat'];
    save([newdatadir,sfilename],'Data','Connectivity','Params')
    allprenom{1,ibgv}=sfilename;
    
    allprenom{2,ibgv}=mln_ifsamegroup(Nsubs,bgv(ibgv,:));
end

function flagSameGroup=mln_ifsamegroup(Nsubs,ibgv)
flagSameGroup = 'den';
subgroup{1} = 1:Nsubs/2;
subgroup{2} = Nsubs/2+1:Nsubs;
if ~isempty(find(subgroup{1}==ibgv(1))) && ~isempty(find(subgroup{2}==ibgv(2)))
    flagSameGroup='few';
end

