%% function produce the training data for the neurofuzzy controller
%% Huifang Wang, Aug 24 2013, Marseille
%% Huifang Wang, Aug 28, 2013, Marseille, generate Nwin-1 more data
%% Huifang Wang, Aug 30, Marseille consider two order indirect connection
%% and common source
%Huifang Wang, 29/10/13, Marseille, bootstrapping results

function [InputMM,InputBS,RMat]=mln_generate_BS_Input(dataprenom,List4Methods,nboot)

% this function will generate input for MULAN inference systems
%[x_1,x_2,x_3,x_4,C_1,C_2, Y], x_1 is the cdf of CS of each links, Y is the normalized reference link
%[0, 1],C_i is the two order indirect connection and common source (1/0),
%which is set to 0

NmethG=length(List4Methods);
methodlog=List4Methods{1};
Ninput=length(methodlog);

% read the files

nbins=100;
datafile=['Tout_',dataprenom,'.mat'];
%Ninput=length(methodlog);
MatA=load(datafile,methodlog{:});
[Nnode,~,Nwin]=size(MatA.(methodlog{1}));
Nlink=Nnode*(Nnode-1);

InputMM=zeros(Nlink,Ninput+3,NmethG);
InputBS=zeros(Nlink,Ninput+3,NmethG,nboot);


for imethG=1:NmethG
    
iTDM=zeros(Nlink,Ninput+3);
methodlog=List4Methods{imethG};
MatA=load(datafile,methodlog{:});
for imethod=1:Ninput
    clear a
    cMethod=methodlog{imethod};
    
    if exist('AllGbs')
        idx=find(strcmp(cMethod,AllGMeth));
        if ~isempty(idx)
             InputMM(:,imethod,imethG)=AllGbs{idx,1};
    InputBS(:,imethod,imethG,:)=AllGbs{idx,2};
    continue;
        end
    end
        
       
        Mat1=MatA.(methodlog{imethod});
    
    [Nnode,~,Nwin]=size(Mat1);
    isizeboot=Nlink;
    Ndata=isizeboot*(nboot+1);
    MatBoot=NaN(Ndata,1);
    iMethodBs=MatBoot;
    
    iMat=mean(abs(Mat1),3);
    iMat(1:size(iMat,1)+1:end) = NaN;
    a=reshape(iMat,Nnode*Nnode,1);
    MatBoot(1:isizeboot,1)=a(~isnan(a));
    
    for iboot=1:nboot
        iBoot_ind=randi(Nwin,[1,Nwin]);
        iBoot_Mat=Mat1(:,:,iBoot_ind);
        iMat=mean(abs(iBoot_Mat),3);
        iMat(1:size(iMat,1)+1:end) = NaN;
        ia=reshape(iMat,Nnode*Nnode,1);
        MatBoot((iboot)*isizeboot+1:(iboot+1)*isizeboot,1)=ia(~isnan(ia));
    end
    xi=linspace(min(MatBoot),max(MatBoot),nbins);
    n_elements = histc(MatBoot,xi);
    c_elements = cumsum(n_elements);
    c_elements=c_elements./max(c_elements);
    %% the means for all windows
    
    for ipairs=1:Ndata
        iMethodBs(ipairs,1)=c_elements(find(xi>=MatBoot(ipairs),1));
    end
    cInputM=iMethodBs(1:Nlink,1);
    cMethodBS=iMethodBs(Nlink+1:Ndata,1);
    cInputBS=reshape(cMethodBS,Nlink,nboot);
    
    if ~exist('AllGbs')
        AllGbs{1,1}=cInputM;
        AllGbs{1,2}=cInputBS;
        AllGMeth{1}=cMethod;
    else
        inew=size(AllGbs,1)+1;
        AllGbs{inew,1}=cInputM;
        AllGbs{inew,2}=cInputBS;
        AllGMeth{inew}=cMethod;
    end
        
    InputMM(:,imethod,imethG)=cInputM;
    InputBS(:,imethod,imethG,:)=cInputBS;
    
    
end
end
load(datafile,'Connectivity');
RMat=Connectivity;
%MatC=mln_removediag(MatC);

%iTDM(:,Ninput+3)=MatC./max(MatC);
       
               
        