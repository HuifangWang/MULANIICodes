%function mln_gen_data_hn(filename,num_hn,subNode)
function mln_gen_data_hn_degree(odirname,ndirname,num_hn,filename,dg)
%Huifang Wang, April 16, 2014
% this function is to delete some channels to creat the cases with hidden
% nodes from filename.mat
% dg degree high/low
%odirname='./';
%ndirname='./';
%dg='l';
%filename='N20nmmCS80S3N4100';
TH_deg=3;
% Marseille, April 18, 2014 for degree
%num_hn=[1 2 5 10];
subNnode=5;
load([odirname,'/',filename]); % Data, Connectivity,Params
Nnode=length(Connectivity);
Nhn=length(num_hn);
%sC=sparse(Connectivity);
%% according to the degrees
[~,~,deg] = mln_degrees_dir(Connectivity);

%richnodes=mln_find_rn(sC,subNnode);
possible_nodes=1:Nnode;
%possiblenodes=mln_del_nodes(possible_nodes,richnodes);

cs=max(max(Connectivity));
allData=Data;
OldConnectivity=Connectivity;
%% correct the filename
indmat=strfind(filename,'.mat');
if ~isempty(indmat)
 filename=filename(1:indmat-1);
end

for idhn=1:Nhn
    Connectivity=OldConnectivity;
    Ndhn=num_hn(idhn);
    Nnode_hn=Nnode-Ndhn;
    switch dg
        case 'h'
            possiblenodes=find(deg>TH_deg);
            for i=1:TH_deg
                if length(possiblenodes)>=Ndhn
                    break;
                else
                    possiblenodes=find(deg>TH_deg-i);
                end
            end
        case 'l'
            possiblenodes=find(deg<TH_deg);
            for i=1:TH_deg
                if length(possiblenodes)>=Ndhn
                    break;
                else
                    possiblenodes=find(deg<TH_deg+i);
                end
            end
    end
    Nposs=length(possiblenodes);
    delnodes=mln_int_randsample(possiblenodes,Ndhn);
    %delnodes=[4 9 10 20 16];
    restnodes=mln_del_nodes(possible_nodes,delnodes);
    savedfilename=[filename,'N',num2str(Nnode_hn),'hn',dg];
    RM1=Connectivity(restnodes,restnodes); % RM1 is the structures which deleted the deleted nodes
    for indhn=1:Ndhn
        % add the common source
        RM2=Connectivity;
        Del_source=Connectivity(:,delnodes(indhn));% we only consider the first order common source
        ind_cs=find(Del_source>0);
        if ~isempty(ind_cs) && length(ind_cs)>=2
            indC = combnk(ind_cs,2);
            for i_indC=1:size(indC,1)
                RM2(indC(i_indC,1),indC(i_indC,2))=cs;
                RM2(indC(i_indC,2),indC(i_indC,1))=cs;
            end
        end
        % add the directly path
        RM3=Connectivity;
        source_Del=Connectivity(delnodes(indhn),:);
        ind_source=find(source_Del>0);
        ind_j=repmat(ind_source,length(ind_cs),1);
        ind_i=sort(repmat(ind_cs,length(ind_source),1));
        for i_i=1:size(ind_i,1)
            RM3(ind_i(i_i),ind_j(i_i))=cs;
        end
        
    end
    
    DRM1=Connectivity;
    DRM1(delnodes,:)=0;
    DRM1(:,delnodes)=0;
    DRM2=RM2;
    DRM2(delnodes,:)=0;
    DRM2(:,delnodes)=0;
    
    DRM3=RM3;
    DRM3(delnodes,:)=0;
    DRM3(:,delnodes)=0;
    
    DRM4=((DRM1+DRM2+DRM3)>0).*cs;
    DRM4(delnodes,:)=0;
    DRM4(:,delnodes)=0;
    
    
    RMmats.delnodes=delnodes;
    RM2=RM2(restnodes,restnodes);
    RM3=RM3(restnodes,restnodes);
    RM4=((RM1+RM2+RM3)>0).*cs;
    Data=allData(restnodes,:);
    
    RMmats.Mat={Connectivity,RM1,RM2,RM3,RM4};
    RMmats.pMat={Connectivity,DRM1,DRM2,DRM3,DRM4};
    RMmats.label={'Original Graph','Deleted nodes','Common Source','Indirect paths','Common Source+Indirect Paths'};
    Connectivity=RMmats;
    
    save([ndirname,'/',savedfilename],'Connectivity','Data','Params')
end


function possiblenodes=mln_del_nodes(possible_nodes,richnodes)
Nrh=length(richnodes);
possiblenodes=possible_nodes;
for i=1:Nrh
    possiblenodes=possiblenodes(possiblenodes~=richnodes(i));
end
